%% Setup

% !!!Specify which session is to be analysed!!!
% Acquisition, Extinction, or Renewal
session = 'Renewal';

% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!! Define likelihood threshold (only DLC data which exceed this
% confidence value will be used for triangulation) !!!
TH = 0.9; 

% !!! Provide filepath for SSM for estimating data
SSM_model_path = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\SSMs\SSM_all_body_points_v1.mat';

%% Get directory for DLC datafiles for each mouse for specified session
% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

% Initialise 
dlc_session_folders = cell(length(mouse_files), 1);

% Get DLC folders for specified session
for mouse = 1:length(mouse_files)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    dlc_session_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Raw DLC Data');
end

%% Sort into mouse order
% Extract mouse number from folder names
mouse_numbers = cellfun(@(x) sscanf(regexp(x, 'Mouse\s*(\d+)', 'match', 'once'), 'Mouse%d'), dlc_session_folders);

% Sort by mouse number and get index for sorted order
[~, sort_idx] = sort(mouse_numbers, 'ascend');

% Apply sorting to folder list
dlc_session_folders = dlc_session_folders(sort_idx);

% Optional: display sorted folder names
disp('Sorted Folder Names:');
disp(dlc_session_folders);

%% Triangulate DLC data and save in appropriate subfolder

% Load camera projection matrix (P)
load('C:\Users\G71044MH\OneDrive - The University of Manchester\Documents\GitHub\3D_camera_calibration\p_matrices\Pcal_hogan_9cameras.mat', 'P')

% Master loop goes over each mouse and accesses its DLC data
for mouse = 1:length(dlc_session_folders)
    current_dlc_folder_path = dlc_session_folders{mouse};
    dlc_file_list = dir(fullfile(current_dlc_folder_path, 'camera*_mouse*.csv'));
    dlc_file_list = dlc_file_list(~[dlc_file_list.isdir]);
    base_names = cell(length(dlc_file_list), 1);
    dlc_file_names = cell(length(dlc_file_list), 1);
    
    % Extract file names from list of DLC files
    for n = 1:length(dlc_file_list)
        dlc_file_names{n} =  dlc_file_list(n).name;
    end
    
    % Loop through each filename and extract base name using regex
    for i = 1:length(dlc_file_names)
    % Match pattern: mouseX_session_pY
        tokens = regexp(dlc_file_names{i}, ['mouse\d+_' session '_p\d+'] , 'match',...
            'ignorecase');
        if ~isempty(tokens)
          base_names{i} = tokens{1};
        end
    end

    % Remove double entries 
    base_names = unique(base_names(~cellfun('isempty', base_names)));

    % Loop through session parts (1 and 2)
    for part = 1 : length(base_names)
    
        % Extract filenames for each camera for each part
        camera_files_struct = dir(fullfile(current_dlc_folder_path, ['camera*' base_names{part} '*.csv'])); 
        % Ensure cameras are in the correct order (matching order of P matrix)
        
        % Load camera file paths into cell array
        camera_files = cell(1, length(camera_files_struct)); % initiate cell array with length corresponding to number of cameras used
    
        for camera = 1:length(camera_files_struct)
            camera_files{camera} = fullfile(current_dlc_folder_path, camera_files_struct(camera).name);
            % -> each cell will contain file paths for CSV files recorded by one camera   
        end 

        % Run load_data function 
        disp(['Loading DLC data for ' base_names{part} '...']);
        [x,L] = load_data_og(camera_files); % Extacts 2D coordinates and their likelihood values and loads them into x and L variables
        
        disp(['Triangulating data for ' base_names{part} '...'])
        % Run triangulation function 
        [X, W] = triangulate_simple_og(x, L, P, TH); 
          

        % Save data in appropriate subfolder and with appropriate name
        % Get mouse names for file creation
        mouse_name = regexp(current_dlc_folder_path, 'Mouse\s*\d+', 'match', 'once');

        triangulated_data_folder = fullfile(master_directory, mouse_name,...
        session, 'Behavioural Data', 'Triangulated Data'); 
        % Skip if save directory cannot be found and display name of
        % directory that needs creating
        if ~exist(triangulated_data_folder, 'dir')
            warning(['Save directory: ' triangulated_data_folder ' does not exist, skipping this file'])
        continue
        end

        triangulated_data_save_path = fullfile(triangulated_data_folder,...
            [base_names{part} '_triangulated.mat']);
        % Check if triangulated data already exists, skip this file if it
        % does
        if exist(triangulated_data_save_path, 'file')
        warning([triangulated_data_save_path ' already exists, skipping save.']);
        else
        save(triangulated_data_save_path,"W","X","x","L","TH","camera_files"); 
        disp([base_names{part} '_triangulated.mat' ' saved to ' triangulated_data_folder])
        end
        
    end

end

% call 'plot3d_video(X,false)' if you want to check triangulation has worked
%% Fit SSM to triangulated data

% Initialise 
triangulated_session_folders = cell(length(mouse_files), 1);

% Get triangulated data folders for specified session
for mouse = 1:length(mouse_files)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    triangulated_session_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Triangulated Data');
end

% Loop through each triangulated data file and fit SSM
for mouse = 1:length(triangulated_session_folders)
    current_triangulated_folder_path = dlc_session_folders{mouse};
    triangulated_file_list = dir(fullfile(current_triangulated_folder_path, 'mouse*_p*_triangulated.mat'));
    SSM_data_save_folder = fullfile(master_directory, mouse_files(mouse).name,...
            session, 'Behavioural Data', 'SSM Fitted Data');
    for part = 1:length(triangulated_file_list)
        SSM_save_name = fullfile([SSM_data_save_folder triangulated_file_list{part},...
            '_SSM_fitted']);
        Fit_SSM_3D_new(triangulated_file_list{part}, SSM_save_name, SSM_model_name_path)
    end
end





%%
% Initialise a cell arrays to store data information
T_matrices = cell(1, length(fileList)); % for each mouse (2 cells - part 1 and part 2) x and y coordinates (row) of body anterior body marker for each frame (colums)
file_names = cell(length(fileList), 1); % for each mouse (2 cells - part 1 and part 2)  filename



%% Extracting data (using SSM fitted data)

% iterate over each file, load it, and extract the T variable

for i = 1:length(fileList)
    filePath = fullfile(folderPath, fileList(i).name);

     % load T(Translation Vector) 
    load(filePath, 'Xfit'); 

    % extract x and y coordinates of anteriour body landmark
    T_matrices{i} = squeeze(Xfit(8,1:2,:));  

    % stores body_anterior x and y positional data in T_matrices array
    file_names{i} = fileList(i).name;
    % stores file names for sorting later
end

%% Sorting velocity matrix into mouse/part order

mouse_numbers = cellfun(@(x) sscanf(regexp(x, 'mouse(\d+)', 'match', 'once'), 'mouse%d'), file_names);
% extracting the mouse number from the file names

[~, sort_idx] = sort(mouse_numbers);
% sort by mouse number amd get index for sorted order

T_matrices = T_matrices(sort_idx); % sort matrix

disp('Sorted File Names:');
disp(file_names(sort_idx));

file_names = file_names(sort_idx); % sort filenames

%% Putting parts of videos together for each mouse
T_matrices_paired = cell(1, length(T_matrices)/2);
%initialising new cell array half the size of original

for video_ix = 1:2:length(T_matrices)
    %iterating over pairs of cells
    T_matrices_paired{(video_ix+1)/2} = cat(2, T_matrices{video_ix}, T_matrices{video_ix+1});
    %stores concatenated matrices in correct indices of new matrix
end

T_matrices = T_matrices_paired;
%% Getting xy positional data and calculating Euclidian distance between frames

num_frames = 502; % 502 frames per trial
num_trials = 40; % 40 trials
num_mice = length(T_matrices); 
all_velocity_data = zeros(num_trials, num_frames, num_mice);

for mouse_idx = 1:num_mice
    x_pos = T_matrices{mouse_idx}(1,:);
    y_pos = T_matrices{mouse_idx}(2,:);

    dist = zeros(num_trials, num_frames);
    
    for trial_ix = 1:num_trials

        dx = [0, diff(x_pos(((trial_ix-1)*num_frames + 1):(trial_ix*num_frames)))];
        dy = [0, diff(y_pos(((trial_ix-1)*num_frames + 1):(trial_ix*num_frames)))];

        % calculating dx and dy for all frames for each trial (e.g., 1:502,
        % 503:1004 etc. and adds a 0 at start to maintain correct trial length)


        dist(trial_ix, :) = sqrt(dx.^2 + dy.^2);
        % calcuting Euclidian distance and adding values to columns of each
        % trial row
    end

all_velocity_data(:, :, mouse_idx) = dist;

end

%% Defining freezing behaviour

%freezing_velocity_threshold = 1.1; % pixels per frame, change as needed (*0.048*15 to get cm/s value)
freezing_velocity_threshold = 0.0528; % cm per frame (use this for xfit data because triangulation transforms to cm)

freezing_duration_threshold = 1*15; % convert seconds to frames (1 second × 15 FPS)


freeze_counter_matrix = all_velocity_data < freezing_velocity_threshold;

validated_freeze_matrix = zeros(size(freeze_counter_matrix));

for mouse_ix = 1:num_mice
    %iterating over each video
    for trials = 1:num_trials 
        % iterating over each trial
        [labeled_freezing, num_segments] = bwlabel(freeze_counter_matrix(trials, :, mouse_ix));
        % labeled_freezing is a vector same size as freeze_counter_matrix, but
        % has a unique value for each 'blob' of consecutive trials
        %num_segments is as scalar representing how many 'blobs' there are
        segment_props = regionprops(labeled_freezing, 'Area');
        % extracts length of connected regions in labeled_freezing and stores
        % it in segment_props.Area
        for segment = 1:num_segments
            %iterates over number of consecutive segments
            if segment_props(segment).Area >=freezing_duration_threshold
                validated_freeze_matrix(trials, labeled_freezing == segment, mouse_ix) = 1;
            end
        end
    end
end
%% Bin freezing data to match neural data
fps = 15; % 15 frames per second
bin_size_frames = 15; %15 frames per bin
bin_size_s = 15/fps;
num_bins = floor(size(all_velocity_data,2)/bin_size_frames);
num_frames = 502;

% Initialise new freeze matrix
binned_freezing_data = zeros(num_trials, num_bins, num_mice);

for mouse = 1:num_mice
    for trial = 1:num_trials
    current_freezing = validated_freeze_matrix(trial, :, mouse);
        for bin = 1:num_bins
            start_idx = (bin-1) * bin_size_frames+1;
            end_idx = bin_size_frames*bin;
            binned_freezing_data(trial, bin, mouse) = mean(current_freezing(start_idx:end_idx));
        end
    end
end

% Classify bins with 70% freezing as freezing bins
is_freezing = double(binned_freezing_data > 0.7); % trials x bins x mice

% Calculate each bin's start time in s
bin_start_s = 0:bin_size_s:(num_bins*bin_size_s)-1; 

disp(['Number of bins: ' num2str(num_bins)])

%% END OF PROCESSING FREEZE DATA %%
%% Neural data directory setup
% Define folder path that contains raw neuropixel data (open ephyis output) of this session
root_ephys_dir = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\'; 

% Define folder path for kilosorted data 
kilosort_dir = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\kilosort4\'; % 

%% Load neural data

dt = 1; % Define neural bin size in seconds
tpre= 0; % Time before event (not applicable here, sampling entire trial)
tpost = num_frames/fps; % Time after (trial end time)
fs = 30000; % Neural sampling rate in Hz

% Generate aligned trigger files
sesh_parts = ['p1' 'p2'];

 for p = 1:length(sesh_parts) % loop through relevant session parts

    % get filepaths to open_ephys data   
    file_path = dir(fullfile(ephys_dir, ['*' sesh_parts{p} '*']));
    TTL_path  = [ephys_dir '\' file_path(1).name '\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\'];
    cont_path = [ephys_dir '\' file_path(1).name '\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\'];
    
    % convert python into matlab files 
    evt    = readNPY([TTL_path   'sample_numbers.npy']);  % extracts trigger counts (one at each video frame) 
    states = readNPY([TTL_path   'states.npy']);          % 1 or -1 (TTL signal on/off)
    cont   = readNPY([cont_path  'sample_numbers.npy']);  % extracts neuronal data sample counts
    
    % adjust trigger events so that align with neuronal data samples
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);

    % save aligned triggers
    save([triggers_path mouse '_' sesh '_' sesh_parts{p}],'evt');
 
  end

load() % just load extinction for now

% evt = double(evt)/fs; % convert sampling number into sec (through dividing by sampling rate)
evt = evt(2:end);     % remove first event (artefact)

samples_per_frame = fs/fps;
counter = 0;

trial_start_fs = zeros(1, num_trials);
for trial = 1:num_trials
    trial_start_fs(trial) = evt(1) + counter * (num_frames * samples_per_frame);
    counter = counter + 1;
end

trial_start_s = trial_start_fs/30000;


% converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

spk = double(spk)/fs;   % converts sample number (by dividing by sampling rate) into seconds
clu_val = unique(clu);  % gets list of all individual clusters, i.e., potential neuronal cells (starts counting at 0)
Ncell = numel(clu_val); % get number of detected clusters 



mfr_all_trials = [];
sfr_all_trials = [];
t_trial =[];
fr_trial = [];

for neuron = 1:Ncell
    tsp = spk(clu == clu_val(neuron));  % extract spike times for current neuron
    [mfr_all_trials(:, neuron), ...
     sfr_all_trials(:, neuron), ...
     t_trial(:, :, neuron), ...
     fr_trial(:, :, neuron), ...
     ~] = raster_NM(tsp, trial_start_s, tpre, tpost, bin_size, false);
end

