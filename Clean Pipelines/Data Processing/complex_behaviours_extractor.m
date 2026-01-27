%%
% !!! Provide directory with all data in !!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!! Provide session being processed !!!
session = 'Renewal';

% !!! Provide filepaths to previously generated RFMs for each behaviour !!!
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\RFM Training Files\Darting\darting_model_1.mat")
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\RFM Training Files\Grooming\grooming_model_1.mat")
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\RFM Training Files\Rearing\rearing_model_1.mat")

received_stim_set_1 = [1 2 3 6 7];
received_stim_set_2 = [4 5 8 9];

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
    mouse_files = mouse_files(sort_idx); 

    % Extract file paths from list of SSM files
    for part = 1:num_parts
        SSM_file_paths{mouse, part} =  fullfile(SSM_file_list(part).folder, SSM_file_list(part).name);
    end
    
end
%% Extract info needed from SSM fitted data

% Initialise storage for each variable
Xfit_matrix_whole = cell(num_mice, num_parts);
b_matrix = cell(num_mice, num_parts);
angle_matrix = cell(num_mice, num_parts);
T_matrix = cell(num_mice, num_parts);
missing_matrix = cell(num_mice, num_parts);

% Extract variables needed from SSM files
for file = 1:num_mice
    for part = 1:num_parts
    load(SSM_file_paths{file, part}, 'Xfit', 'b', 'R', 'T', 'missing');
    Xfit_matrix_whole{file, part} = Xfit;
    b_matrix{file, part} = b;
    angle_matrix{file, part} = rotm2eul(R, 'ZYX')';
    T_matrix{file, part} = T;
    missing_matrix{file, part} = double(missing);
    end
end

% Initialise concantenated storage for each variable
Xfit_matrix_whole_paired = cell(num_mice, num_parts/2);
b_matrix_paired = cell(num_mice, num_parts/2);
R_matrix_paired = cell(num_mice, num_parts/2);
T_matrix_paired = cell(num_mice, num_parts/2);
missing_matrix_paired = cell(num_mice, num_parts/2);

%  Concantenated parts for each variable
for mouse = 1:num_mice
    Xfit_matrix_whole_paired{mouse} = cat(3, Xfit_matrix_whole{mouse, 1}, Xfit_matrix_whole{mouse, 2});
    b_matrix_paired{mouse} = cat(2, b_matrix{mouse, 1}, b_matrix{mouse, 2});
    R_matrix_paired{mouse} = cat(2, angle_matrix{mouse, 1}, angle_matrix{mouse, 2});
    T_matrix_paired{mouse} = cat(2, T_matrix{mouse, 1}, T_matrix{mouse, 2});
    missing_matrix_paired{mouse} = cat(2, missing_matrix{mouse, 1}, missing_matrix{mouse, 2});
end
%% Extract features and make predictions

behaviours = {'Grooming', 'Rearing', 'Darting'};

models = {grooming_model, rearing_model, darting_model};

behaviour_thresholds  = [15, 7, 7]; % Minimum time in frames the behaviour lasts for

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

for behaviour_idx = 1:length(behaviours)
    current_behaviour = behaviours{behaviour_idx};
    current_model = models{behaviour_idx};
    current_threshold = behaviour_thresholds(behaviour_idx);
    extracted_features_matrix = cell(num_mice, 1); 

    for mouse = 1:num_mice
    
        Xfit = Xfit_matrix_whole_paired{mouse};
        b = b_matrix_paired{mouse};
        R = R_matrix_paired{mouse};
        T = T_matrix_paired{mouse};
        missing = missing_matrix_paired{mouse};
        
        [extracted_features] = extract_behavioural_features(Xfit, b, R, T, missing, current_behaviour);
    
        extracted_features_matrix{mouse} = extracted_features;
    
        % Initialise 
        behaviour_session_folders = cell(length(mouse_files), 1);
        
        % Get extracted behaviour folders for specified session
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        behaviour_session_folders{mouse} = fullfile(mouse_path, mouse_name,...
            session, 'Behavioural Data', 'Extracted Behaviours', current_behaviour);
        
         % Dynamically construct save path
         save_name_looms = lower(['mouse', num2str(mouse), '_', session, '_', ...
            current_behaviour, '_looms']);
        save_name_scores_looms = lower(['mouse', num2str(mouse), '_', session, '_', ...
            current_behaviour, '_probabilities_looms']);
        save_path_looms = fullfile(behaviour_session_folders{mouse}, [save_name_looms, '.mat']);
        save_path_scores_looms = fullfile(behaviour_session_folders{mouse}, [save_name_scores_looms, '.mat']);
        
        save_name_flashes = lower(['mouse', num2str(mouse), '_', session, '_', ...
            current_behaviour, '_flashes']);
        save_name_scores_flashes = lower(['mouse', num2str(mouse), '_', session, '_', ...
            current_behaviour, '_probabilities_flashes']);
        save_path_flashes = fullfile(behaviour_session_folders{mouse}, [save_name_flashes, '.mat']);
        save_path_scores_flashes = fullfile(behaviour_session_folders{mouse}, [save_name_scores_flashes, '.mat']);
       
       
        % Predict frames with behaviour
        [predictions, scores] = predict(current_model, extracted_features_matrix{mouse});
        predicted_labels = str2double(predictions);
        scores = scores(:, 2);

        % Post-processing: remove false positives
        % Identify continuous segments of behaviour
        [labeled_behaviour, num_segments] = bwlabel(predicted_labels);
        % Extract length of each segment
        segment_props = regionprops(labeled_behaviour, 'Area');
   
        % Iterate over detected segments
        for segment = 1:num_segments
            % Replace short sequences with 0s
            if segment_props(segment).Area < current_threshold    
                predicted_labels(labeled_behaviour == segment) = 0; 
            end
        end

        predicted_labels = reshape(predicted_labels, 502, 40)'; 
        scores = reshape(scores, 502, 40)'; 
        
        % Select trial type using index provided at start of section
        if ismember(mouse, received_stim_set_1)         
            behaviour_looms = predicted_labels(stim_set_1_looms_idx, :);
            behaviour_flashes = predicted_labels(stim_set_1_flashes_idx, :);
            looms_scores = scores(stim_set_1_looms_idx, :);
            flashes_scores = scores(stim_set_1_flashes_idx, :);
        elseif ismember(mouse, received_stim_set_2)            
            behaviour_looms = predicted_labels(stim_set_2_looms_idx, :);            
            behaviour_flashes = predicted_labels(stim_set_2_flashes_idx, :);         
            looms_scores = scores(stim_set_2_looms_idx, :);         
            flashes_scores = scores(stim_set_2_flashes_idx, :);
        end
                
        % Save predictions
        % if exist(save_path_looms, 'file')
        %     warning('%s already exists, skipping extraction', save_name_looms)
        % else
            disp(['Making predictions for ', save_name_looms, '...'])
            save(save_path_looms, 'behaviour_looms');
            disp([save_name_looms, ' has been saved!'])       
        % end
        
        % if exist(save_path_flashes, 'file')
        %     warning('%s already exists, skipping extraction', save_name_flashes)
        % else
            disp(['Making predictions for ', save_name_flashes, '...'])
            save(save_path_flashes, 'behaviour_flashes');
            disp([save_name_flashes, ' has been saved!'])       
        % end
        
        % Save behaviour probabilities for GLM covariate input
        % if exist(save_path_scores_looms, 'file')
        %     warning('%s already exists, skipping probability extraction', save_name_scores_looms)
        % else
            save(save_path_scores_looms, 'looms_scores');
            disp([save_name_scores_looms, ' has been saved!'])
        % end
        
        % if exist(save_path_scores_flashes, 'file')
        %     warning('%s already exists, skipping probability extraction', save_name_scores_flashes)
        % else
            save(save_path_scores_flashes, 'flashes_scores');
            disp([save_name_scores_flashes, ' has been saved!'])
        % end

    
    end

end

