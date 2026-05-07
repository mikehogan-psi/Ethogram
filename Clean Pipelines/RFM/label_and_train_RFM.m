%% Setup
% !!! Provide top directory with all data in !!!
% master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';
master_directory = 'W:\Mike\Neuropixels_Fear_Conditioning\Data';

% !!! Provide path to where RFM files are !!!
RFM_folder = 'W:\Mike\Neuropixels_Fear_Conditioning\Models\RFMs';

% !!! Provide session !!!
session = 'Extinction'; % Extinction or Renewal

% !!! Provide behaviour being labelled !!!
behaviour = 'Freezing'; % Rearing, Grooming or Darting (or Freezing for feature extraction only)

% Provide number of mouse to label
mouse_to_label = 7; 

% Provide video part to label
part = 1;

% Provide model number
model_number = 3;

do_save = true;
%% Extract data and setup filepaths
mouse_files = dir(fullfile(master_directory, 'mouse*'));

num_mice = size(mouse_files, 1);

video_data_path = cell(num_mice, 1);

for mouse_idx = 1:num_mice
    mouse_name = mouse_files(mouse_idx).name;
    mouse_path = mouse_files(mouse_idx).folder;
    current_mouse_path = fullfile(mouse_path, mouse_name, session, 'Behavioural Data',...
        'Video Data');
    video_data_path{mouse_idx} = current_mouse_path; 
end


mouse_files = dir(fullfile(master_directory, 'mouse*'));

SSM_data_path = cell(num_mice, 1);

for mouse_idx = 1:num_mice
    mouse_name = mouse_files(mouse_idx).name;
    mouse_path = mouse_files(mouse_idx).folder;
    current_mouse_path = fullfile(mouse_path, mouse_name, session, 'Behavioural Data',...
        'SSM Fitted Data');
    SSM_data_path{mouse_idx} = current_mouse_path; 
end

labels_output_folder = fullfile(RFM_folder, behaviour, 'Labelled Data');
models_output_folder = fullfile(RFM_folder, behaviour, 'Models');
%% Load GUI and label data
video_paths = cell(4, 1);

video_list = dir(fullfile(video_data_path{mouse_to_label},...
    ['camera*_mouse' num2str(mouse_to_label) '*_p' num2str(part) '*']));

video_paths{1} = fullfile(video_data_path{mouse_to_label}, video_list(1).name);
video_paths{2} = fullfile(video_data_path{mouse_to_label},video_list(4).name);
video_paths{3} = fullfile(video_data_path{mouse_to_label},video_list(6).name);
video_paths{4} = fullfile(video_data_path{mouse_to_label}, video_list(9).name);

base_name = regexp(video_paths{1}, 'mouse\s*(\d+).*?_p(\d+)', 'match');
base_name = base_name{1};
base_name = lower(strrep(base_name, ' ', ''));

SSM_path = dir(fullfile(SSM_data_path{mouse_to_label}, ['*p' num2str(part) '*']));

SSM_data_path = fullfile(SSM_path.folder, SSM_path.name);

load(SSM_data_path, 'Xfit', 'b', 'R', 'T', 'missing')

[labels, features] = label_behaviour_gui(video_paths, Xfit, b, R, T, missing, behaviour, labels_output_folder, base_name);
%% Train model with extracted features/labels

% Load all labeled data and put into single variables
labelled_files = dir(fullfile(labels_output_folder, '*_labels.mat'));

num_videos = length(labelled_files);

filtered_data = cell(num_videos, 1);

% Remove any row that has a feature = NaN
for i = 1:length(labelled_files)
    data = load(fullfile(labels_output_folder, labelled_files(i).name));
    good_idx  = ~isnan(data.extracted_features);
    good_idx = all(good_idx, 2);
    data.extracted_features = data.extracted_features(good_idx, :); % Extract valid features
    data.labels = data.labels(good_idx);
    filtered_data{i}.extracted_features = data.extracted_features;
    filtered_data{i}.labels = data.labels;
    disp(['Total usable frames labelled from video ' num2str(i),...
        ': ', num2str(length(data.labels))]);
end

% Define number of folds
k = num_videos;
numTrees = 500; % Number of trees for Random Forest

% Initialize performance metrics
accuracy_scores = zeros(k, 1);
d_prime_scores = zeros(k, 1);
precision_scores = zeros(k,1);
recall_scores    = zeros(k,1);
f1_scores        = zeros(k,1);
balacc_scores    = zeros(k,1);
mcc_scores       = zeros(k,1);
prauc_scores     = zeros(k,1);
rocauc_scores    = zeros(k,1);
brier_scores     = zeros(k,1);

prev_scores      = zeros(k,1);
baseline_brier_scores = zeros(k,1);
brier_skill_scores = zeros(k,1);

p1_all_test = [];
y_all_test = [];

% Test model using k-folds
for fold = 1:k
    % Get train/test indices for this fold
    k_range = 1:k;
    test_idx = fold;  % Extract sequential video per fold for testing  
    train_idx = k_range ~= test_idx; % Others are train data

    train_data = filtered_data(train_idx); % Apply these indexes to filtered data
    test_data = filtered_data{test_idx};
    
    all_train_features = [];
    all_train_labels = [];
    
    
    for i = 1:length(train_data)
        features = train_data{i}.extracted_features;
        labels = train_data{i}.labels;
        all_train_features = [all_train_features; features];
        all_train_labels = [all_train_labels; labels];
    end

    all_test_features = test_data.extracted_features;
    all_test_labels = test_data.labels;

    % Naming variables for readability

    % Prepare training and testing data for model training
    X_train = all_train_features;
    y_train = all_train_labels;
    X_test = all_test_features;
    y_test = all_test_labels;

    % Train Random Forest Model
    model = TreeBagger(numTrees, X_train, y_train, 'Method', 'classification');
    

    % Predict on test set
    [y_pred, scores] = predict(model, X_test);
    y_pred = str2double(y_pred);

    % get probability of class "1"
    positive_col = find(strcmp(model.ClassNames,'1'));
    p1 = scores(:,  positive_col);

    p1_all_test = [p1_all_test; p1];
    y_all_test = [y_all_test; y_test(:)];
        
    % Calulate model performance metrics per fold
    [accuracy_score, d_prime_score, precision_score, recall_score, f1_score,...
    balacc_score, mcc_score, prauc_score, rocauc_score, brier_score] = ...
    test_RFM(p1, y_pred, y_test);
    
    accuracy_scores(fold, 1) = accuracy_score;
    d_prime_scores(fold, 1) = d_prime_score;
    precision_scores(fold,1) = precision_score;
    recall_scores(fold,1) = recall_score;
    f1_scores(fold,1) = f1_score;
    balacc_scores(fold,1) = balacc_score;
    mcc_scores(fold,1) = mcc_score;
    prauc_scores(fold,1) = prauc_score;
    rocauc_scores(fold,1) = rocauc_score;
    brier_scores(fold, 1) = brier_score;
    
    prev_scores(fold,1) = mean(y_test == 1);
    baseline_brier_scores(fold,1) = prev_scores(fold) * (1 - prev_scores(fold));
    brier_skill_scores(fold,1) = 1 - (brier_scores(fold) / baseline_brier_scores(fold));

end

% ===== Display summary metrics (mean ± SEM across folds) =====
k_eff = numel(accuracy_scores);  % in case k changes

mean_sem = @(x) deal(mean(x,'omitnan'), std(x,'omitnan')/sqrt(sum(isfinite(x))));

[mAcc, seAcc] = mean_sem(accuracy_scores);
[mDp,  seDp ] = mean_sem(d_prime_scores);
[mPre, sePre] = mean_sem(precision_scores);
[mRec, seRec] = mean_sem(recall_scores);
[mF1,  seF1 ] = mean_sem(f1_scores);
[mBal, seBal] = mean_sem(balacc_scores);
[mMCC, seMCC] = mean_sem(mcc_scores);
[mPRA, sePRA] = mean_sem(prauc_scores);
[mROC, seROC] = mean_sem(rocauc_scores);
[mBri, seBri] = mean_sem(brier_scores);

[mPrev, sePrev] = mean_sem(prev_scores);
[mBSS,  seBSS ] = mean_sem(brier_skill_scores);

disp('==================== CV PERFORMANCE SUMMARY ====================');
disp(['Folds used: ' num2str(k_eff)]);
disp(['Prevalence    : ' num2str(100*mPrev,'%.2f') ' ± ' num2str(100*sePrev,'%.2f') ' %']);
disp(['Accuracy      : ' num2str(100*mAcc,'%.2f') ' ± ' num2str(100*seAcc,'%.2f') ' %']);
disp(['d-prime       : ' num2str(mDp,'%.3f')      ' ± ' num2str(seDp,'%.3f')]);
disp(['Precision (1) : ' num2str(mPre,'%.3f')     ' ± ' num2str(sePre,'%.3f')]);
disp(['Recall (1)    : ' num2str(mRec,'%.3f')     ' ± ' num2str(seRec,'%.3f')]);
disp(['F1 (1)        : ' num2str(mF1,'%.3f')      ' ± ' num2str(seF1,'%.3f')]);
disp(['Balanced Acc  : ' num2str(mBal,'%.3f')     ' ± ' num2str(seBal,'%.3f')]);
disp(['MCC           : ' num2str(mMCC,'%.3f')     ' ± ' num2str(seMCC,'%.3f')]);
disp(['PR-AUC        : ' num2str(mPRA,'%.3f')     ' ± ' num2str(sePRA,'%.3f')]);
disp(['ROC-AUC       : ' num2str(mROC,'%.3f')     ' ± ' num2str(seROC,'%.3f')]);
disp(['Brier         : ' num2str(mBri,'%.3f')     ' ± ' num2str(seBri,'%.3f')]);
disp(['Brier skill   : ' num2str(mBSS,'%.3f')      ' ± ' num2str(seBSS,'%.3f') ' (vs baseline p(1-p))']);
disp('===============================================================');

% Collate all data to train final model
all_labels = [];
all_features = [];

for i = 1:num_videos
    current_data = filtered_data{i};
    all_features = [all_features; current_data.extracted_features]; 
    all_labels = [all_labels; current_data.labels]; 
end

% Calculate final averaged results
final_accuracy = mean(accuracy_scores);
final_d_prime = mean(d_prime_scores);

disp(['Final Cross-Validation Accuracy: ', num2str(final_accuracy * 100), '%']);
disp(['Final Cross-Validation d-prime: ', num2str(final_d_prime)]);

model_name = [lower(behaviour), '_model'];
trained_model = TreeBagger(numTrees, all_features, all_labels, 'Method', 'classification');

model_struct = struct();
model_struct.(model_name) = trained_model;


if do_save
    save_name = [lower(behaviour), '_model_' num2str(model_number) '.mat'];
    save(fullfile(models_output_folder, save_name), '-struct', 'model_struct', model_name);
    disp([save_name, ' saved to ', models_output_folder]);
    
    % Save metrics variables
    save(fullfile(models_output_folder, save_name), ...
        'accuracy_scores','d_prime_scores','precision_scores','recall_scores','f1_scores', ...
        'balacc_scores','mcc_scores','prauc_scores','rocauc_scores', ...
        'prev_scores','brier_scores','baseline_brier_scores','brier_skill_scores', 'p1_all_test', ...
        'y_all_test','-append');
end