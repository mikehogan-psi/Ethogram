%%
% !!! Provide directory with all data in !!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!! Provide session being processed !!!
session = 'Renewal';

% !!! Provide filepaths to previously generated RFMs for each behaviour !!!
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\RFM Training Files\Darting\darting_model_1.mat")
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\RFM Training Files\Grooming\grooming_model_1.mat")
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\RFM Training Files\Rearing\rearing_model_1.mat")



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

SSM_file_paths = cell(num_mice, 1);



for mouse = 1:length(SSM_session_folders)
    current_SSM_folder = SSM_session_folders{mouse};
    SSM_file_list = dir(fullfile(current_SSM_folder, 'mouse*_*habituation_SSM_fitted.mat'));
    SSM_file_list = SSM_file_list(~[SSM_file_list.isdir]);    
    
    [~, idx] = sort({SSM_file_list.name});
    SSM_file_list = SSM_file_list(idx);
    

    % Extract file paths from list of SSM files
    SSM_file_paths{mouse} =  fullfile(SSM_file_list.folder, SSM_file_list.name);
    
end


%% Extract info needed from SSM fitted data

% Initialise storage for each variable
Xfit_matrix_whole = cell(num_mice, 1);
b_matrix = cell(num_mice, 1);
angle_matrix = cell(num_mice, 1);
T_matrix = cell(num_mice, 1);
missing_matrix = cell(num_mice, 1);

% Extract variables needed from SSM files
for file = 1:num_mice
    load(SSM_file_paths{file}, 'Xfit', 'b', 'R', 'T', 'missing');
    Xfit_matrix_whole{file} = Xfit;
    b_matrix{file} = b;
    angle_matrix{file} = rotm2eul(R, 'ZYX')';
    T_matrix{file} = T;
    missing_matrix{file} = double(missing);
end

%% Extract features and make predictions

behaviours = {'Grooming', 'Rearing', 'Darting'};

models = {grooming_model, rearing_model, darting_model};

behaviour_thresholds  = [15, 7, 7]; % Minimum time in frames the behaviour lasts for


for behaviour_idx = 1:length(behaviours)
    current_behaviour = behaviours{behaviour_idx};
    current_model = models{behaviour_idx};
    current_threshold = behaviour_thresholds(behaviour_idx);
    extracted_features_matrix = cell(num_mice, 1); 

    for mouse = 1:num_mice
    
        Xfit = Xfit_matrix_whole{mouse};
        b = b_matrix{mouse};
        R = angle_matrix{mouse};
        T = T_matrix{mouse};
        missing = missing_matrix{mouse};
        
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
         save_name_hab = lower(['mouse', num2str(mouse), '_', session, '_', ...
            current_behaviour, '_habituation']);
        save_name_scores_hab = lower(['mouse', num2str(mouse), '_', session, '_', ...
            current_behaviour, '_habituation_probabilities']);
        save_path_hab = fullfile(behaviour_session_folders{mouse}, [save_name_hab, '.mat']);
        save_path_scores_hab = fullfile(behaviour_session_folders{mouse}, [save_name_scores_hab, '.mat']);
        
      
       
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
        
        behaviour_habituation = predicted_labels;
        habituation_scores = scores;
                
        % Save predictions
        % if exist(save_path_looms, 'file')
        %     warning('%s already exists, skipping extraction', save_name_looms)
        % else
            disp(['Making predictions for ', save_name_hab, '...'])
            save(save_path_hab, 'behaviour_habituation');
            disp([save_name_hab, ' has been saved!'])       
        % end
        
        % Save behaviour probabilities for GLM covariate input
        % if exist(save_path_scores_looms, 'file')
        %     warning('%s already exists, skipping probability extraction', save_name_scores_looms)
        % else
            save(save_path_scores_hab, 'habituation_scores');
            disp([save_name_scores_hab, ' has been saved!'])
        % end
        
      
    end

end

