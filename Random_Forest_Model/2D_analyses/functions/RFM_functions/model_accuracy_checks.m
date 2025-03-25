%% Checking accuracy of model

%load the required data
% labels predicted by model (predicted_labels)
    load('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\extinction_analysis\rearing_data\predicted_labels_data\predicted_rearing_labels_mouse13_extinction_p2_2024-10-16-131736-0000DLC_resnet50_Fear Extinction No ImplantOct15shuffle1_500000_body_fit.mat');
% manually corrected labeles (verified_labels)
    load('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\extinction_analysis\rearing_data\corrected_labels_data\mouse13_extinction_p2_2024-10-16-131736-0000_verified_labels.mat');

% identify false positives and false negatives   
    false_neg = (predicted_labels == 0) & (verified_labels == 1);
    nr_false_neg = sum(false_neg);
    false_pos = (predicted_labels == 1) & (verified_labels == 0);
    nr_false_pos = sum(false_pos);

% total model prediction accuracy 
    total_errors = nr_false_pos + nr_false_neg;
    model_accuracy = 100 - (total_errors / length(predicted_labels));


%% finding missing data from extracting features

% load SSM fitted data (A,b,C,T,X,Ffit and missing)
load('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\data_all_mice\SSM_fitted_data\SSM_fitted_data_extinction\mouse13_extinction_p2_2024-10-16-131736-0000DLC_resnet50_Fear Extinction No ImplantOct15shuffle1_500000_body_fit.mat');

num_frames = length(verified_labels);
all_features = [];


    for k = 1: num_frames 
        
    window_size = 30;  % Define window size (/fps to get time over which features are being computed)
    start_idx = max(1, k - floor(window_size / 2));
    end_idx = min(num_frames, k  + floor(window_size / 2));

    % Extract data for the window
    b_window = b(:, start_idx:end_idx);
    T_window = T(:, start_idx:end_idx);
    A_window = A(start_idx:end_idx);

    % Calculate velocity and angular velocity, these can be changed to
    % capture the desired behaviour      
    angular_velocity = [0, abs(diff(A_window))];  % Prepend 0 for alignment
    velocity = [0, sqrt(sum(diff(T_window, 1, 2).^2))];  % Euclidean velocity

    % Compute features for the current window
    frame_features = [
        mean(b_window(3, :));  % Mean of eigenpose 3 (looking up/down)
        var(b_window(3, :));   % Variance of eigenpose 3
        mean(diff(b_window(3, :)));  % Temporal derivative of eigenpose 3
        mean(b_window(2, :));  % Mean of eigenpose 2 (elongation/hunching)
        mean(velocity);        % Mean velocity
        mean(angular_velocity) % Mean angular velocity
    ];

    all_features = [all_features, frame_features];
    end


% identify frames that features could not be computed for (hence have NAN)
 missing_features = any(isnan(all_features), 1);
 num_missing_features = sum(missing_features);
 indx_missing_features  = find(missing_features);                   %indexes of frames missing features
 missing_features_false_pos = sum(false_pos(indx_missing_features));%number of false postives arrising from missing features  
 missing_features_false_neg = sum(false_neg(indx_missing_features));%number of false negatives arrising from missing features  


%% Checking influence of missing datapoints

% load SSM fitted data (A,b,C,T,X,Ffit and missing)
load('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\data_all_mice\SSM_fitted_data\SSM_fitted_data_extinction\mouse13_extinction_p2_2024-10-16-131736-0000DLC_resnet50_Fear Extinction No ImplantOct15shuffle1_500000_body_fit.mat');

%extract details of missing data points
total_missing = sum(missing);                    %total number of missing datapoints 
indx_missing  = find(missing);                   %indexes of missing datapoints
missing_false_pos = sum(false_pos(indx_missing));%number of false postives arrising from missing datapoints  
missing_false_neg = sum(false_neg(indx_missing));%number of false negatives arrising from missing datapoints  

%% Checking amount of missing data in general

 folder_path = ('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\data_all_mice\SSM_fitted_data\SSM_fitted_data_extinction');
 all_files = dir(fullfile(folder_path, '*mouse*'));
 num_mice_files = length(all_files);

 for i = 1:num_mice_files
    
     filename = all_files(i).name;
     load(filename);

     total_missing = sum(missing);  

     mouse_id = regexp(filename, 'mouse(\d+)', 'tokens', 'once');
     mouse_id = str2double(mouse_id);

    fprintf('Mouse %d had %d missing datapoints.\n', mouse_id, total_missing);

end





