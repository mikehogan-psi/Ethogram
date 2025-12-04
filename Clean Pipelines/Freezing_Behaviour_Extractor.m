%% Setup (input needed)
% !!!Specify which session is to be analysed!!!
% Acquisition, Extinction, or Renewal
session = 'Renewal';

% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

%!!!Specify stim set for each mouse (mouse number only)!!!
received_stim_set_1 = [1 2 3 6 7];
received_stim_set_2 = [4 5 8];

%% Get directory for SSM datafiles for each mouse for specified session
% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

% Initialise 
SSM_session_folders = cell(length(mouse_files), 1);

% Get SSM folders for specified session
for mouse = 1:length(mouse_files)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    SSM_session_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'SSM Fitted Data');
end

%% Sort into mouse order
% Extract mouse number from folder names
mouse_numbers = cellfun(@(x) sscanf(regexp(x, 'Mouse\s*(\d+)', 'match', 'once'), 'Mouse%d'), SSM_session_folders);

% Sort by mouse number and get index for sorted order
[~, sort_idx] = sort(mouse_numbers, 'ascend');

% Apply sorting to folder list
SSM_session_folders = SSM_session_folders(sort_idx);
mouse_files = mouse_files(sort_idx); 

% Optional: display sorted folder names
disp('Sorted Folder Names:');
disp(SSM_session_folders);
%% Extracting SSM data

num_mice = size(mouse_files, 1);

% Number of parts in video for each mouse (always 2)
num_parts = 2;

SSM_file_paths = cell(num_mice, num_parts);

for mouse = 1:length(SSM_session_folders)
    current_SSM_folder = SSM_session_folders{mouse};
    SSM_file_list = dir(fullfile(current_SSM_folder, 'mouse*_p*_SSM_fitted.mat'));
    SSM_file_list = SSM_file_list(~[SSM_file_list.isdir]);    

    [~, idx] = sort({SSM_file_list.name});
    SSM_file_list = SSM_file_list(idx);

    % Extract file paths from list of SSM files
    for part = 1:num_parts
        SSM_file_paths{mouse, part} =  fullfile(SSM_file_list(part).folder, SSM_file_list(part).name);
    end
    
end

%% Extract central body landmark from Xfit and concatenate parts of videos
% Initialise storage for Xfit body anterior data (num mice x sesion parts)
Xfit_matrix_central = cell(num_mice, num_parts);

for file = 1:length(SSM_file_paths)
    for part = 1:2
    load(SSM_file_paths{file, part}, 'Xfit');  
    % Extract x and y for 8th landmark (body anterior)
    Xfit_matrix_central{file, part} = squeeze(Xfit(8,1:2,:));
    end
end

% Putting parts of videos together for each mouse
% Initialising new cell array half the size of original (part-wise)
Xfit_paired_central = cell(num_mice, num_parts/2);

for mouse = 1:num_mice
    Xfit_paired_central{mouse} = cat(2, Xfit_matrix_central{mouse, 1}, Xfit_matrix_central{mouse, 2});
end
%% Getting xy positional data and calculating Euclidian distance between frames

num_frames = 502; % 502 frames per trial
num_trials = 40; % 40 trials
all_velocity_data = zeros(num_trials, num_frames, num_mice);

for mouse_idx = 1:num_mice
    x_pos = Xfit_paired_central{mouse_idx}(1,:);
    y_pos = Xfit_paired_central{mouse_idx}(2,:);
    x_pos = fillmissing(x_pos, 'linear');
    y_pos = fillmissing(y_pos, 'linear');


    dist = zeros(num_trials, num_frames);
    
    for trial_ix = 1:num_trials

        % Calculating delta x and delta y for all frames for each trial (e.g., 1:502,
        % 503:1004 etc. and adds a 0 at start to maintain correct trial length)

        dx = [0, diff(x_pos(((trial_ix-1)*num_frames + 1):(trial_ix*num_frames)))];
        dy = [0, diff(y_pos(((trial_ix-1)*num_frames + 1):(trial_ix*num_frames)))];
      
        % Calcuting Euclidian distance and adding values to columns of each
        % trial row
        dist(trial_ix, :) = sqrt(dx.^2 + dy.^2);
        
    end

all_velocity_data(:, :, mouse_idx) = dist;

end

nan_idx = isnan(all_velocity_data);
all_velocity_data(nan_idx) = 0;

%% Defining freezing behaviour

fps = 15; % fps is always 15

% cm per frame (use this for xfit data because triangulation transforms to cm)
freezing_velocity_threshold = 0.0528; 

% Convert seconds to frames (1 second × 15 FPS)
freezing_duration_threshold = 1*fps; 

freeze_counter_matrix = all_velocity_data < freezing_velocity_threshold;

validated_freeze_matrix = zeros(size(freeze_counter_matrix));

% Iterating over each video
for mouse_ix = 1:num_mice
    % Iterating over each trial    
    for trials = 1:num_trials 
        % Labeled_freezing is a vector same size as freeze_counter_matrix, but
        % has a unique value for each 'blob' of consecutive trials.
        % num_segments is as scalar representing how many 'blobs' there are
        [labeled_freezing, num_segments] = bwlabel(freeze_counter_matrix(trials, :, mouse_ix));
        % Extracts length of connected regions in labeled_freezing and stores
        % it in segment_props.Area
        segment_props = regionprops(labeled_freezing, 'Area');
        % Iterating over number of consecutive segments
        for segment = 1:num_segments
            if segment_props(segment).Area >=freezing_duration_threshold
                validated_freeze_matrix(trials, labeled_freezing == segment, mouse_ix) = 1;
            end
        end
    end
end

%% Splitting in trial type and saving to respective folders

% Stim sets: 0 = flash, 1 = loom
stim_set_1 = [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1,...
    0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1];

stim_set_2 = [0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0,...
    1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0 ,0 ,1 ,0];

% Setting indexes for looms and flashes
stim_set_1_looms_idx = logical(stim_set_1);
stim_set_1_flashes_idx = ~stim_set_1_looms_idx;
stim_set_2_looms_idx = logical(stim_set_2);
stim_set_2_flashes_idx = ~stim_set_2_looms_idx;

% Loop through each filename and extract base name for creating save paths
base_names = cell(length(SSM_session_folders), 1);
for i = 1:length(SSM_session_folders)
% Match pattern: mouseX_session_pY
    tokens = regexp(SSM_file_paths{i, 1}, ['mouse\d+_' session] , 'match',...
        'ignorecase');
    if ~isempty(tokens)
      base_names{i} = tokens{1};
    end
end

% Remove double entries 
base_names = unique(base_names(~cellfun('isempty', base_names)));

for mouse = 1:num_mice
    current_freeze_data = validated_freeze_matrix(:, :, mouse);

    % Select trial type using index provided at start of script
    if ismember(mouse, received_stim_set_1)
        loom_freezing = current_freeze_data(stim_set_1_looms_idx, :);
        flash_freezing = current_freeze_data(stim_set_1_flashes_idx, :);
    elseif ismember(mouse, received_stim_set_2)
        loom_freezing = current_freeze_data(stim_set_2_looms_idx, :);
        flash_freezing = current_freeze_data(stim_set_2_flashes_idx, :);
    end
    
    % Extract folder for saving freeze data
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    freeze_data_folder = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Freezing');
    
    % Create savepath for each file
    freeze_data_save_path_flash = fullfile(freeze_data_folder,...
            [base_names{mouse} '_flashes_freezing.mat']);

    freeze_data_save_path_loom = fullfile(freeze_data_folder,...
            [base_names{mouse} '_looms_freezing.mat']);
    
    % % Save (or skip if already created)
    % if exist(freeze_data_save_path_flash, 'file')
    %     warning([base_names{mouse} '_flashes_freezing.mat already exists, skipping save.']);
    % else
        disp(['Saving ', base_names{mouse}, '_flashes_freezing.mat'])
        save(freeze_data_save_path_flash, 'flash_freezing')
    % end

    % if exist(freeze_data_save_path_loom, 'file')
    %     warning([base_names{mouse} '_looms_freezing.mat already exists, skipping save'])
    % else
        disp(['Saving ', base_names{mouse}, '_looms_freezing.mat'])
        save(freeze_data_save_path_loom, 'loom_freezing')
    % end

end

%% Save velocity data

for mouse = 1:num_mice
    current_velocity_data = all_velocity_data(:, :, mouse);

    % Select trial type using index provided at start of script
    if ismember(mouse, received_stim_set_1)
        loom_velocity = current_velocity_data(stim_set_1_looms_idx, :);
        flash_velocity = current_velocity_data(stim_set_1_flashes_idx, :);
    elseif ismember(mouse, received_stim_set_2)
        loom_velocity = current_velocity_data(stim_set_2_looms_idx, :);
        flash_velocity = current_velocity_data(stim_set_2_flashes_idx, :);
    end
    
    % Extract folder for saving velocity data
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    velocity_data_folder = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Velocity');
    
    % Create savepath for each file
    velocity_data_save_path_flash = fullfile(velocity_data_folder,...
            [base_names{mouse} '_flashes_velocity.mat']);

    velocity_data_save_path_loom = fullfile(velocity_data_folder,...
            [base_names{mouse} '_looms_velocity.mat']);
    
    % % Save (or skip if already created)
    % if exist(velocity_data_save_path_flash, 'file')
    %     warning([base_names{mouse} '_flashes_velocity.mat already exists, skipping save.']);
    % else
        disp(['Saving ', base_names{mouse}, '_flashes_velocity.mat'])
        save(velocity_data_save_path_flash, 'flash_velocity')
    % end

    % if exist(velocity_data_save_path_loom, 'file')
    %     warning([base_names{mouse} '_looms_velocity.mat already exists, skipping save'])
    % else
        disp(['Saving ', base_names{mouse}, '_looms_velocity.mat'])
        save(velocity_data_save_path_loom, 'loom_velocity')
    % end

end



%% Extracting freezing data for habituation

%% Extracting SSM data

num_mice = size(mouse_files, 1);

SSM_file_paths = cell(num_mice, 1);

for mouse = 1:length(SSM_session_folders)
    current_SSM_folder = SSM_session_folders{mouse};
    SSM_file_list = dir(fullfile(current_SSM_folder, 'mouse*_habituation_SSM_fitted.mat'));
    SSM_file_list = SSM_file_list(~[SSM_file_list.isdir]);    

    [~, idx] = sort({SSM_file_list.name});
    SSM_file_list = SSM_file_list(idx);

    % Extract file paths from list of SSM files

    SSM_file_paths{mouse} =  fullfile(SSM_file_list.folder, SSM_file_list.name);
    
end

for file = 1:length(SSM_file_paths)
    load(SSM_file_paths{file}, 'Xfit');  
    % Extract x and y for 8th landmark (body anterior)
    Xfit_matrix_central{file} = squeeze(Xfit(8,1:2,:));
end

hab_velocity_data = cell(num_mice, 1);

for mouse_idx = 1:num_mice
    x_pos = Xfit_matrix_central{mouse_idx}(1,:);
    y_pos = Xfit_matrix_central{mouse_idx}(2,:);
    x_pos = fillmissing(x_pos, 'linear');
    y_pos = fillmissing(y_pos, 'linear');

    dx = diff(x_pos);
    dy = diff(y_pos);

    hab_velocity_data{mouse_idx} = sqrt(dx.^2 + dy.^2);
    nan_idx = isnan(hab_velocity_data{mouse_idx});
    hab_velocity_data{mouse_idx}(nan_idx) = 0;
    
end


fps = 15; % fps is always 15

% cm per frame (use this for xfit data because triangulation transforms to cm)
freezing_velocity_threshold = 0.0528; 

% Convert seconds to frames (1 second × 15 FPS)
freezing_duration_threshold = 1*fps; 

%% Thresholding & validating freezing (cell-based)

freeze_counter_matrix   = cell(num_mice, 1);   % logical per-frame freezing
validated_freeze_matrix = cell(num_mice, 1);   % logical per-frame validated freezing

for mouse_idx = 1:num_mice
    % Logical vector: frames below velocity threshold
    freeze_counter_matrix{mouse_idx} = hab_velocity_data{mouse_idx} < freezing_velocity_threshold;

    % Label consecutive freezing segments for this mouse
    [labeled_freezing, num_segments] = bwlabel(freeze_counter_matrix{mouse_idx});

    % Preallocate validated vector (same length as this session)
    validated = false(size(freeze_counter_matrix{mouse_idx}));

    if num_segments > 0
        % Get length of each connected segment
        segment_props = regionprops(labeled_freezing, 'Area');

        % Keep only segments at least as long as freezing_duration_threshold
        for segment = 1:num_segments
            if segment_props(segment).Area >= freezing_duration_threshold
                validated(labeled_freezing == segment) = true;
            end
        end
    end

    % Store per-mouse validated freezing vector
    validated_freeze_matrix{mouse_idx} = validated;
end

%% Saving

% Loop through each filename and extract base name for creating save paths
base_names = cell(num_mice, 1);
for i = 1:num_mice
    % Match pattern: mouseX_session
    tokens = regexp(SSM_file_paths{i}, 'mouse\d+_' , 'match', ...
        'ignorecase');
    if ~isempty(tokens)
        base_names{i} = tokens{1};
    end
end

for mouse = 1:num_mice
    % Now a logical vector for this mouse, variable length OK
    hab_freeze_data = validated_freeze_matrix{mouse};

    % Extract folder for saving freeze data
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    freeze_data_folder = fullfile(mouse_path, mouse_name, ...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Freezing');

    save_name = [base_names{mouse}, lower(session) '_habituation_freezing.mat'];

    if ~exist(fullfile(freeze_data_folder, save_name), 'file')
        save(fullfile(freeze_data_folder, save_name), "hab_freeze_data");
        fprintf('%s has been saved to %s\n', save_name, freeze_data_folder);
    else
        warning('%s already exists, skipping save', save_name);
    end
end
