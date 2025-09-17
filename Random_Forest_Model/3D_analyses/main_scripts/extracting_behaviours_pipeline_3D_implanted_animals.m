% current work in progress %%%%

%% Behavioural Analysis Pipeline to define an predict specific behaviours (for 3D data)

% This analysis pipeline allows for the identification (a.k.a. extraction) of desired mouse behaviours ("target behaviour")

% input: 3D coordinates of triangulated data (use triangulation3D_DLC_data to generate 3D data from 2D DLC (Deep lab cut) labeled data)
% end-outout: frames in which desired behaviour is exhibited

% Individual analysis steps:
%            Step 0: General SETUP                                  - define directories, model names and general variables
%            STEP 1: Generate Statistical Shape model (SSM)         - using 3D-upper function to correct outliers and generate mean pose & eigenposes
%            STEP 2: Generate Statistical Shape model fitted data   - uses SSM to calculate mouse body shape parameter measures 
%            STEP 3: SETUP calculate features function              - define which features (derived from  body shape parameters) and time window shall be used to train model to predict behaviour
%            STEP 4: Manually label frames of 3-4 videos            - label frames that exibit desried behaviour (to later train model with)
%            STEP 5: Train Random Forest Model                      - associates the previously defined features with behaviour
%            STEP 6: predict labels in full dataset                 - trained model to predict behaviour in each frame based on exhibited features
%            STEP 7: Post-processing                                - Remove small groups of false positive labels
%            STEP 8: Manual check predicted label accuracy and correct labels  
%            STEP 9: Asessment of model senstivity                  - manual check of predcited labels for estimation of model accuracy and sensitivity     



% Behaviour definitions (for concistency the behaviours were consequently defined if at least one of the points was fullfilled)
    % rearing and exploratory behaviour (non-anxiety):
    %   - standing on hind legs
    %   - looking up
    %   - front paws touching side of arena AND looking up
    %   - average timewindow = 30s

    % darting behaviour (anxiety):
    %   - short, quick bursts of running and/or abrupt changes in direction
    %   - freezing before and after darting
    %   - running into courner of box 
    %   - average timewindow = 30s

    % freezing behaviour (anxiety):
    %   - mody completly immmobile for at least a few seconds

%% Step 0: General SETUP 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % !!! CHANGE/ ADJUST THESE AS NECCECSARY !!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Define variables that remain concistent throughout entire analyses

     % A: mouse part that shall be analysied (depends on behaviour that will be investigated)
            which_part = 'body';    % for behaviours that entail changes in body posture e.g. rearing
            % which_part = 'tail';   % for behaviours that entail changes in tail movements e.g. tail rattling 

     % B: define target behaviour that shall be analysed
            behaviour = 'tail_rattling'; 
    
     % C: define session 
          sesh = 'extinction';
       %    sesh = 'renewal';


     % D: define which name generated models shall be given
            SSM_model_name = 'SSM_all_body_points_v1.mat';      % Can be left out if SSM fitted data already generated
            RFM_model_name = ['RFM_' behaviour '_1'];          % this model will be used to predict behaviour -> folders in which behavioral data is stored is named after the model
                                                               % !!! make a note of all the configurations used for this model (see Random forest model configurations document)                                                                
                             
                                                                
% 2. Define directories to where specific data can be accessed/saved 
     % A: folder containing 'raw' data from all mice (data before any behaviours predicted)  
            % common_raw_dir = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\data_all_mice\Extinction\data_from_DLC_iteration_3'; % Laptop
            %common_raw_dir = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\data_all_mice\Renewal'; % Laptop
            common_raw_dir = 'C:\Cohort 4 Temp Data Storage\data_all_mice';
     % subdirectories containing video and DLC data (must be already defined) + folder where SSM fitted data shall be saved 
            % ! make sure these folders already exist and contain correct data prior to starting analysis !
            video_path    = [common_raw_dir '\video_data'];                  % raw video data (videos of mice) 
            DLC_data_path = [common_raw_dir '\DLC_data'];                    % DLC data(2D coordinates of body label markers recorded from different camera angles - generated by DeepLabCut 
            triangulated_data_path = [common_raw_dir '\triangulated_data'];  % contains 3D coordinates of body label markers - generated trough triangulation
            SSM_data_path = [common_raw_dir '\SSM_fitted_data'];             % SSM fitted data (3D-Upper generated data - corrected body marker positions and mean/eigenpose data)
          
      
     % B: folder where data labeled with desired behaviour shall be saved   
              common_behav_dir = 'C:\Cohort 4 Temp Data Storage\data_all_mice'; % Laptop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % ! will be automatically created for each new analyses (based on the previous definitions !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % 1. subdirectories containing behaviour labeled data 
    
        manual_labels_path    = [common_behav_dir '\' behaviour '_analyses\' RFM_model_name '\manually_labeled_data\'];   % data with frames manually labeled as behaviour (to train model)
                                 mkdir(manual_labels_path);
        predicted_labels_path = [common_behav_dir '\' behaviour '_analyses\' RFM_model_name '\predicted_labels_data\'];   % data with behaviour labels predicted by model
                                 mkdir(predicted_labels_path);
        corrected_labels_path = [common_behav_dir '\' behaviour '_analyses\' RFM_model_name '\corrected_labeled_data\'];   % predcited_labvelds data which was manually corrected        
                                 mkdir(corrected_labels_path);

                                 addpath(genpath([common_behav_dir '\' behaviour '_analyses\' RFM_model_name ])); % adds new model folder and subfolders to matlab path

 % 2. Define directories to where models shall be saved (both the SSM and RFM will be saved here)
        model_path =  [common_behav_dir '\models\']; 
                       mkdir(model_path);
            SSM_model_name_path = [model_path SSM_model_name]; % automatically defines path to find SSM model
            RFM_model_name_path = [model_path RFM_model_name]; % automatically defines path to find SSM model


 % 3. Extract DLC data files and count them 
        all_files = dir(fullfile(triangulated_data_path, 'mouse*')); % makes a struct containing details of all DLC files 
        num_files = length(all_files);                   % identifies number of files (for each mouse two files so 2x nr. of mice)

%% STEP 1: Generate Statistical Shape model (SSM) - skip if SSM fitted data already generated
% here the triangulated 3D data is used to generate a SSM of the mouse using the 3D_upper analysis pipeline
% performs following operations: 1.removes outliers
%                                2.performs data alignment (computes mean pose i.e. average shape of mouse)
%                                3.applies PCA to model shape variation (computes eigenposes which represent the primary modes of variation from the mean pose)

% run function that will estimate SSM for moude body or tail 
        SSM_est_mike_abi_implanted(triangulated_data_path, SSM_model_name_path) 


% the SSM model will be saved with the following variables:
% - template: contains mean pose i.e. aligned and centered mean 3D coordinates of each body marker [Np,3] 
% - lambda: Eigenvalues, describe amount of variance explained by one eigenvector 
% - eignV3D: 3D reshaped version of the eigenvectors [Np × 3 × Eigenposes]
% - sigma2: tell you how much of the variation in the data is not explained by the eigenvectors (if small most of the variation is explained by the SSM

% results of PCA (the tree main eigenposes that together explain 80% of the data variation from mean pose):
%      SSM_3D_implant_mouse1_head_fixed (body): -> model made from data of implanted mouse1
                    % Eigenpose 1: "bending side to side" of body
                    % Eigenpose 2: "up/down movement" of head with some body bending 
                    % Eigenpose 3: "elongation/hunching" of body 


%% STEP 2: Generate Statistical Shape model fitted data (using 2D-Upper trained model) - skip if SSM fitted data already generated

% the SSM is used to generate mouse body shape parameter measures for every single frame :
%     X(Input Coordinates)          - The original input coordinates after applying the necessary filters (body or tail). 
%     Xtf/Xfit(Fitted Coordinates)  - The reconstructed 3D coordinates of the shape after fitting the Statistical Shape Model (after smoothing data)
%     Ttf/T(Translation Vector)     - The 3D translation vector that aligns the fitted model to the observed data, describes the offset in x, y and z directions.
%     Rtf/R(Rotation Angle)         - rotational alignment of the fitted Statistical Shape Model (SSM) with the observed data for each frame (via 3x3 rotation matrix)
%     btf/b(Shape Coefficients)     - These are the weights for the eigenvectors in the shape model, describe how much each principal component contributes to the reconstructed shape.
% (the ...tf variable forms have data split into infividual trials)


% run SSM on all 3D data files to generate body shape parameters
 cd(triangulated_data_path) 

for f = 1 : length(all_files) -1

    current_3Ddata_file = all_files(f).name;

    Fit_SSM_3D_new(current_3Ddata_file, SSM_data_path, SSM_model_name_path)
    % Fit_SSM_tail_3D(current_3Ddata_file, SSM_data_path, SSM_model_name_path)
end

% can manually check how fitted data looks like (just load a specific SSM_fit.mat)
  plot3d_video(Xfit,false)


%% STEP 3: SETUP calculate features function
 
%%%%%%%%%%%!!! DEFINE VARIABLES TO CALCULATE RELEVANT FEATURES !!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. out of the following list of all possible to calculate features decide which onces are relevant for the target behaviour
%    ! for detailed explanation of features see mouse_feature_selectiontable !
%    # can also modify this and add new features #

 frame_features_all_strings = {
 
  'mean(b_window(3, :))';                          % Feature 1:  Average "elongation/hunching" (General horizontal body deformation )
  'var(b_window(3, :))';                           % Feature 2:  Variability in "elongation/hunching"  ( horizontal deformation dynamism / general fluctuation)
  'mean(diff(b_window(3, :)))';                    % Feature 3:  Trend of change "elongation/hunching"  posture (rate and  direction of horizontal deformation dynamism)
  'mean(b_window(1, :))';                          % Feature 4:  Average "bending side to side" (general lateral deformation) 
  'var(b_window(1, :))';                           % Feature 5:  Variability in "bending side to side" (lateral motion dynamism)
  'mean(diff(b_window(1, :)))';                    % Feature 6:  Trend of change in "bending side to side" (rate and  direction of lateral motion dynamism)
  'mean(b_window(2, :))';                          % Feature 7:  Average "looking up/down" (General vertical posture trend)
  'var(b_window(2, :))';                           % Feature 8:  Variability in "looking up/down" posture (vertical motion dynamism / general fluctuation)
  'mean(diff(b_window(2, :)))';                    % Feature 9:  Trend of change in "looking up/down" posture (rate and  direction of vertical motion dynamism)

  'mean([0, sqrt(sum(diff(T_window(1:2,:), 1, 2).^2, 1))])';                           % Feature 10: Average 2D velocity - ONLY USING X AND Y COODINATES !!!
  'mean([0, abs(diff(R_window(1,:)))])'                                                % Feature 11: Average Angular Velocity / Average yaw (Describes rotational movement around z-axis)
  'var([0, diff([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])])'                 % Feature 12: Variability in acceleration of movement / How much acceleration fluctuates (2D)
  'sum(abs([0, diff([0, diff([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])])]))' % Feature 13: Cumulative change in acceleration (2D)  
  'max([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])'                            % Feature 14: Peak velocity within the window (2D)
  'sum(abs(diff([0, diff(b_window(1, :))])))'                                          % Feature 15: Cumulative change in rate of bending
  'sum([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])'                            % Feature 16: Total distance travelled (2D)
  'mean([0, abs(diff(R_window(2,:))) + abs(diff(R_window(3,:)))])';                    % Feature 17: average angular velocity (magnitude only) in pitch and roll directions combined (Describes rotational movement around x- and y-axis)


  'Xfit(4,3,frame_idx)-Xfit(11,3,frame_idx)'                                         % Feature 18: z-difference between nose (body-marker 4) and tail-anterior (body-marker 11)
  'mean([0, sqrt(sum(diff(T_window(1:2,:), 1, 2).^2, 1))])';                         % Feature 19: Average 3D velocity (overall movement intensity) - ONLY USING X AND Y COODINATES
  'sum([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])'                          % Feature 20: Total distance travelled (3D)
    
  'sum(abs(diff(T_window(3, :))))'                	     % Feature 21: cumulative vertical movement 
  'Xfit(7,3,frame_idx)-Xfit(10,3,frame_idx)'             % Feature 22: z-difference between neck base (body-marker 7) and tail-base (body-marker 10)
  'norm(Xfit(7,1:2,frame_idx) - Xfit(10,1:2,frame_idx))' % Feature 23: euclidian distance in x-y plane bewteen neck base and tail base 
  'Xfit(4,3,frame_idx)'                                  % Feature 24: z-value of nose
  'max(T_window(3, :))'                                  % Feature 25: max vertical translation
  'Xfit(1,3,frame_idx)'                                  % Feature 26: z-value of cable tip
  'Xfit(4,3,frame_idx)-Xfit(7,3,frame_idx)'              % Feature 27: z-difference between nose (body-marker 4) and tail-anterior (body-marker 11) - implanted animals

 };

% 2. define which features are relevant for the target behaviour
    
    % feature_selectection_indx = [5,7,8,9,17,18,19,20,21,22,23,24,25,26];  % feature selection RFM_rearing_1
    % feature_selectection_indx = [21, 24, 25, 26, 27];  % feature selection RFM_rearing_5 and 6

    % feature_selectection_indx = [2,5,6,11,12,13,14];  % feature selection RFM_darting_1
    % feature_selectection_indx = [2, 5, 6, 11, 12, 13, 14, 23 ];  % feature selection RFM_darting_4
    % feature_selectection_indx = [11, 12, 14, 15, 16];  % feature selection RFM_darting_9

    %feature_selectection_indx = [2,5,8,10,12,16]; % freezing 1;
    %feature_selectection_indx = [2,3,5,6,8,9,10,12,16]; % freezing 2;
    % feature_selectection_indx = [2,5,6,10,11,12,13,14,16,21]; % freezing 3;

    % feature_selectection_indx = [2 5 9 24 21 27 20]; - grooming

    frame_features_strings = frame_features_all_strings(feature_selectection_indx); % automatically excludes all features that are not relevant for behaviour


% 3. define the timewindow which is relevant to the target behaviour in frames (sampling rate = 15fps - e.g. 30 frames = 2s )
       % window_size = 30; % time window for RFM_rearing_1 & 2
       %  window_size = 4; % time window for RFM_rearing_5 and 6
       %  window_size = 26; % darting_9
        window_size = 20; % freezing_2

%% STEP 4(A): Manually label frames that exibit desried behaviour (to later train model with)

% define file name of raw video and mat file containing results of SSM (2D_UPPER output)
base_name = ['mouse3_' sesh '_p1']; % !!! CHANGE THIS !!!

% find specific camera filepaths
all_camera_files = dir(fullfile(video_path, ['camera*' base_name, '*.avi']));
video_files_paths{1} = [video_path '\' all_camera_files(1).name];
video_files_paths{2} = [video_path '\' all_camera_files(4).name];
video_files_paths{3} = [video_path '\' all_camera_files(6).name];
video_files_paths{4} = [video_path '\' all_camera_files(9).name];

% Find the correct SSM(.mat) file
SSM_file_path = dir(fullfile(SSM_data_path, [base_name, '*.mat']));
SSM_file_path = [SSM_data_path '\' SSM_file_path.name];

% run feature_extractor function
feature_extractor_3D(base_name, video_files_paths, SSM_file_path , manual_labels_path, frame_features_strings, window_size); 
% open GUI to scrub through video and lable all frames that mouse is showing desired behaviour (eg. rearing)

% Saved output files contain a lables and features variable for each mouse
    % lables: vector containing 1s for all labled frames and 0s for balanced number of non-labled frames (randomly selects x-times (e.g. 5) non-labeled frames) 
    % features: contains x columns (depending on how many features analysed) each with values representing de degree to which a cetrain behavioural feature is expresseed in a certain frame (per row)
         

%% STEP 4(B): extract balanced features from already labeled videos (to later train model with) 
% in the case that you have aleady previousely labelled some videos with target behaviour and now want to extact new features 
% (eg. added or removed a feature to see if it improves the model) then run this section instead of STEP 4(A) 

% go into to folder that contains .mat files with behavior labels for all frames 
% (either manually_labeled_data or corrected_labeled_data) 

files = dir(fullfile(cd, '*mouse*')); % get list of all files within the folder to loop through

for f = 1 : length(files) 

    % load file
    load(files(f).name);
    
        % incase corrected_labels data loaded variable has to be renamed
        if exist('verified_labels', 'var')
            all_frames_labels = verified_labels;
            clear verified_labels; % Optional: remove the original variable
        end
    
    % extract basename
    base_name = regexp(files(f).name, ['(mouse\d+_' sesh '_p\d+)'], 'match', 'once');  
    
    % Find the correct video and SSM(.mat) file
    SSM_file_path = dir(fullfile(SSM_data_path, [base_name, '*.mat']));
    SSM_file_path = [SSM_data_path '\' SSM_file_path.name];
    
    % run feature_extractor_already_labeled function
    feature_extractor_already_labeled_3D(base_name, all_frames_labels, SSM_file_path, manual_labels_path, frame_features_strings, window_size)

end 

% Clear unnecessary variables
  clearvars all_frames_labels base_name features labels 


%% STEP 5: Train Random Forest Model using 5-Fold Cross-Validation

% Load all labeled data
labelled_files = dir(fullfile(manual_labels_path, '*_labels.mat'));
all_features = [];
all_labels = [];

for i = 1:length(labelled_files)
    data = load(fullfile(manual_labels_path, labelled_files(i).name));
    all_features = [all_features; data.features(~isnan(data.labels), :)]; % Extract valid features
    all_labels = [all_labels; data.labels(~isnan(data.labels))];
end

disp(['Total usable frames labelled: ', num2str(length(all_labels))]);

% Define number of folds
k = 5;
cv = cvpartition(all_labels, 'KFold', k, 'Stratify', true);

% Initialize performance metrics
accuracy_scores = zeros(k, 1);
d_prime_scores = zeros(k, 1);
numTrees = 100; % Number of trees for Random Forest

for fold = 1:k
    % Get train/test indices for this fold
    train_idx = training(cv, fold);
    test_idx = test(cv, fold);
    
    % Split data
    X_train = all_features(train_idx, :);
    y_train = all_labels(train_idx);
    X_test = all_features(test_idx, :);
    y_test = all_labels(test_idx);
    
    % Train Random Forest Model
    model = TreeBagger(numTrees, X_train, y_train, 'Method', 'classification');
    
    % Predict on test set
    y_pred = predict(model, X_test);
    y_pred = str2double(y_pred);
    
    % Calculate accuracy
    accuracy_scores(fold) = sum(y_pred == y_test) / length(y_test);
    
    % Calculate d-prime
    hit_rate = sum(y_pred == 1 & y_test == 1) / sum(y_test == 1);
    false_alarm_rate = sum(y_pred == 1 & y_test == 0) / sum(y_test == 0);
    
    % Prevent extreme hit/false alarm rates (avoid infinity in z-score)
    hit_rate = max(min(hit_rate, 0.99), 0.01);
    false_alarm_rate = max(min(false_alarm_rate, 0.99), 0.01);
    
    % Compute d-prime
    d_prime_scores(fold) = (norminv(hit_rate) - norminv(false_alarm_rate)) / sqrt(2);
    
    disp(['Fold ', num2str(fold), ' Accuracy: ', num2str(accuracy_scores(fold) * 100), '%']);
    disp(['Fold ', num2str(fold), ' d-prime: ', num2str(d_prime_scores(fold))]);
end

% Calculate final averaged results
final_accuracy = mean(accuracy_scores);
final_d_prime = mean(d_prime_scores);

disp(['Final Cross-Validation Accuracy: ', num2str(final_accuracy * 100), '%']);
disp(['Final Cross-Validation d-prime: ', num2str(final_d_prime)]);

% Save the final trained model using all data
final_model = TreeBagger(numTrees, all_features, all_labels, 'Method', 'classification');
save([RFM_model_name_path '.mat'], 'final_model');


%% STEP 6: Use Model on to predict labels in full dataset

% Load SSM files for all mice
file_list_SSM = dir(fullfile(SSM_data_path, '*_fit.mat')); % change if body or tail data needed

% Load model that shall be used to predict labels
loaded_model = load(RFM_model_name_path);
model = loaded_model.final_model; % load saved model, make sure to delete old one if retraining

for i = 1:length(file_list_SSM)

    % get the mouseID to name predicted_label file at the end
    filename = file_list_SSM(i).name;
    new_file_name = erase(filename, '_3D_SSM_fit.mat');

    % load the SSM data
    mat_path = fullfile(SSM_data_path, file_list_SSM(i).name);
    data = load(mat_path);
    
    b = data.b;
    R = data.R;
    T = data.T;
    missing = data.missing;
    Xfit = data.Xfit; 

    % calculate Euler angles from R (rotational matrix) -> [yaw, pitch, roll] / [yaw, roll, pitch] for each frame
    R = rotm2eul(R, 'ZYX');
    R = R';

  
    % define nr of frames and features
    num_frames = size(b, 2);
    all_frames_features = zeros(num_frames, length(feature_selectection_indx));    

    
    % Extract behaviour-relevant features for every frame based on specified time-window 
    % (identical to calculation of features in feature_extractor function)
    for frame_idx = 1:num_frames

     % Define window size that seems relevant to behaviour being analysed
        start_idx = max(1, frame_idx - (window_size / 2));
        end_idx = min(num_frames, frame_idx + (window_size / 2));

     % Extract data for the window
        b_window = b(:, start_idx:end_idx);
        T_window = T(:, start_idx:end_idx);
        R_window = R(:, start_idx:end_idx);

    % find how many frames are missing datapoints
        miss_window = missing(start_idx:end_idx);      
        miss_num = sum(miss_window); % Count number of frames qith missing datapoints

        % If 1-6 frames with missing datapoints, exclude these  
        % (If 7 or more columns have NaNs, we proceed without excluding them)
        if miss_num >= 1 && miss_num <= 6
            valid_cols = ~any(isnan([b_window; T_window; R_window]), 1);
            b_window = b_window(:, valid_cols);
            T_window = T_window(:, valid_cols);
            R_window = R_window(:, valid_cols);
        end   
       
        % compute all features for this frame  
        for s = 1:length(frame_features_strings)
            all_frames_features(frame_idx,s) = eval(frame_features_strings{s});
        end

    end

    % Predict Rearing Frames
    predictions = predict(model, all_frames_features);
    predicted_labels = str2double(predictions);
    % rearing_frames = find(predicted_labels == 1);
    
    % Save predicted behaviour frame indices 
    save([predicted_labels_path 'predicted_' behaviour '_labels_' new_file_name '.mat'], 'predicted_labels'); 
    disp('Predicted rearing frames saved.');
end

%% STEP 7: Post-processing - Remove small groups of false positive labels
% false postives often occur singluar or small groups of random labelled frames or frames labelled in close proximity to behavioural event
%   -> to correct for these all behavioural labels which occur in less than a certain amount of consectuve frames are removed 

% Load predicted labels files for all mice
file_list_pred_labels = dir(fullfile(predicted_labels_path, ['predicted_' behaviour '_labels_*'])); % change if body or tail data needed

for i = 1:length(file_list_pred_labels) % loop through all files

    % load predicted labels variable
    load([predicted_labels_path file_list_pred_labels(i).name]);

    % Identify continuous segments of 1s (behavior)
      [labeled_behaviour, num_segments] = bwlabel(predicted_labels);
    % Extract properties of each segment (length of consecutive 1s)
      segment_props = regionprops(labeled_behaviour, 'Area');
    
    % Define minimum duration threshold !!! DEFINE THE CUTOFF BASED ON TYPICAL DURATION OF BEHAVIOUR !!! 
            %duration_threshold = 10;  % for darting_8
            %duration_threshold = 7;  % for rearing
            duration_threshold = 15;   % for freezing 1

    % Iterate over detected darting sequences
    for segment = 1:num_segments
        if segment_props(segment).Area < duration_threshold    

            predicted_labels(labeled_behaviour == segment) = 0; % Replace short darting sequences with 0s
        end
    end

    % diplay progress
    disp(['label groups <' num2str(duration_threshold) 'filtered for' file_list_pred_labels(i).name ])

    % save the filtered predicted labels variable to replace original file
      save([predicted_labels_path file_list_pred_labels(i).name], 'predicted_labels'); 

end


%% STEP 8: Manual check predicted label accuracy and correct labels  

% choose a file you want to double-check the mouse behaviour of the frames predicted by the model to show the desired behaviour 
base_name = ['mouse1_' sesh '_p1'];

    % Find the correct video and predicted labels file
    video_file_path = dir(fullfile(video_path, ['camera4_' base_name, '*.avi']));
    video_file_path = [video_path '/' video_file_path.name];
    % predicted_file_path = dir(fullfile(predicted_labels_path, ['predicted_' behaviour '_labels_', base_name, '*.mat']));
    % predicted_file_path = [predicted_labels_path '/' predicted_file_path.name];
    % 
 predicted_file_path =  'C:\Cohort 4 Temp Data Storage\data_all_mice\grooming_analyses\RFM_grooming_1\predicted_labels_data\predicted_grooming_labels_mouse1_extinction_p1.mat';

    % fun validate predcitions function
         validate_predictions(base_name, video_file_path, predicted_file_path, corrected_labels_path);
         % should open video in which can manually scrub through frames and see if the predicted labels match actual behaviour and if necessary correct it
    
    % repeat this for a couple of videos to generate enough verified_labels data to update model (to increase its accuracy)     

%% STEP 9: Assess model accuracy and sensitivity 

% choose a video which you have the manual or verified labels and predicted labels for 
base_name = ['*mouse2_' sesh '_p2'];

    % load the correct label variable (either from manually labelled or corrected label files)
    matFiles = dir(fullfile(corrected_labels_path, [base_name, '*.mat']));

            if isempty(matFiles)
               matFiles = dir(fullfile(manual_labels_path, [base_name, '*.mat']));
               load([manual_labels_path, matFiles.name]);
               correct_labels = all_frames_labels;  
    
            else     
                load([corrected_labels_path, matFiles.name]);  
                correct_labels = verified_labels; 
            end
    
   % load the by model predicted label variable      
   matFiles = dir(fullfile(predicted_labels_path, ['predicted_' behaviour '_labels_', base_name, '*.mat']));
              load([predicted_labels_path matFiles.name]);  

    % Calculate accuracy
    accuracy_score = sum(predicted_labels == correct_labels) / length(predicted_labels);
    
    % Calculate d-prime
    hit_rate = sum(predicted_labels == 1 & correct_labels == 1) / sum(correct_labels == 1);
    false_alarm_rate = sum(predicted_labels == 1 & correct_labels == 0) / sum(correct_labels == 0);
    
    % Prevent extreme hit/false alarm rates (avoid infinity in z-score)
    hit_rate = max(min(hit_rate, 0.99), 0.01);
    false_alarm_rate = max(min(false_alarm_rate, 0.99), 0.01);
    
    % Compute d-prime
    d_prime_score = (norminv(hit_rate) - norminv(false_alarm_rate)) / sqrt(2);
    
    disp([base_name ' ' behaviour ' predictions:']);
    disp(['Accuracy = ', num2str(accuracy_score * 100), '%']);
    disp(['d-prime = ', num2str(d_prime_score)]);



