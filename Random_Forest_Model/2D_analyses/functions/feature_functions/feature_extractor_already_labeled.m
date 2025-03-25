function feature_extractor_already_labeled(all_frames_labels, SSM_file_path, manual_labels_path, video_file_path, frame_features_strings, window_size)

% labeled_data_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\extinction_analysis\rearing_data\corrected_labels_data\mouse13_extinction_p2_2024-10-16-131736-0000_verified_labels.mat';
% mat_file_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\data_all_mice\SSM_fitted_data\SSM_fitted_data_extinction\mouse13_extinction_p2_2024-10-16-131736-0000DLC_resnet50_Fear Extinction No ImplantOct15shuffle1_500000_body_fit.mat';
% output_folder = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\extinction_analysis\rearing_data\labeled_frames_data\';
% video_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\data_all_mice\video_data\video_data_extinction\mouse13_extinction_p2_2024-10-16-131736-0000.avi';

    % Load the labeles
    labels = all_frames_labels;

    % Load data
    data = load(SSM_file_path);
    b = data.b;
    A = data.A;
    T = data.T;
    missing = data.missing; 

    % Initialise feature matrix 
    num_frames = length(data.A);
    features = [];  % Store features for all frames as rows


    % balance dataset and extract corresponding features
        behaviour_indices = find(labels == 1);
        non_behaviour_indices = find(labels == 0);

        % select non-behaviour frames (want eg. 5x more than behaviour frames)
        num_behaviour = numel(behaviour_indices);
        num_non_behaviour_to_sample = min(5 * num_behaviour, numel(non_behaviour_indices)); % to have a 5-to-1 non-behaviour/behaviour ratio
                                                                                            % OR if not enough non-behaviour frames to reach ratio -> just takes all non-behaviour frames
        selected_non_behaviour_indices = datasample(non_behaviour_indices, num_non_behaviour_to_sample, 'Replace', false);

        % Combine indices
        balanced_indices = [behaviour_indices; selected_non_behaviour_indices];

        % Extract features for the balanced dataset
        balanced_features = arrayfun(@calculate_features, balanced_indices, 'UniformOutput', false);
        features = cell2mat(balanced_features');
        features = features';
        labels = labels(balanced_indices);

        % Remove NaN values
        valid_rows = all(~isnan(features), 2);
        features = features(valid_rows, :);
        labels = labels(valid_rows);

        disp('Balanced feature set generated.');

    % save balanced labels and corresponding features    
        % Ensure the output folder path ends with a separator
        if ~endsWith(manual_labels_path, filesep)
            manual_labels_path = [manual_labels_path, filesep];
        end

        % Create the full path for saving the file
        [~, video_name, ~] = fileparts(video_file_path);
        save_path = [manual_labels_path, video_name, '_labels.mat']; 

        try     
            save(save_path, 'labels', 'features', 'all_frames_labels');
            disp(['Labels and features saved to ', save_path]);
        catch ME
            disp(['Error saving labels: ', ME.message]);
        end
   

    % Feature extraction function 
        function frame_features = calculate_features(frame_idx)
        
        % Define window size that seems relevant to behaviour being analysed
        start_idx = max(1, round(frame_idx - (window_size / 2)));
        end_idx = min(num_frames, round(frame_idx + (window_size / 2)));

        % Extract data for the window
        b_window = b(:, start_idx:end_idx);
        T_window = T(:, start_idx:end_idx);
        A_window = A(start_idx:end_idx);

        % find how many frames are missing datapoints
        miss_window = missing(start_idx:end_idx);
        % Count number of frames qith missing datapoints
        miss_num = sum(miss_window);

        % If 1-6 frames with missing datapoints, exclude these  
        % (If 7 or more columns have NaNs, we proceed without excluding them)
        if miss_num >= 1 && miss_num <= 6
            valid_cols = ~any(isnan([b_window; T_window; A_window]), 1);
            b_window = b_window(:, valid_cols);
            T_window = T_window(:, valid_cols);
            A_window = A_window(valid_cols);
        end   
        
        % compute all features for this frame  
        frame_features = zeros(length(frame_features_strings), 1); % Preallocate for efficiency

        for i = 1:length(frame_features_strings)
            frame_features(i) = eval(frame_features_strings{i});
        end


    end

end
