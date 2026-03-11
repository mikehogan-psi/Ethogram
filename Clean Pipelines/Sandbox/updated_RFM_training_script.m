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
k = 4;

% Initialize performance metrics
accuracy_scores = zeros(k, 1);
d_prime_scores = zeros(k, 1);
numTrees = 100; % Number of trees for Random Forest

for fold = 1:k
    % Get train/test indices for this fold
    k_range = 1:k;
    test_idx = fold;  % Extract one video for testing  
    train_idx = k_range ~= test_idx;

    train_data = filtered_data(train_idx);
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
    y_pred = predict(model, X_test);
    y_pred = str2double(y_pred);
    
    % Calculate accuracy
    accuracy_scores(fold) = sum(y_pred == y_test) / length(y_test);
    
    % Confusion counts
    hits  = sum(y_pred == 1 & y_test == 1);
    miss  = sum(y_pred == 0 & y_test == 1);
    fa    = sum(y_pred == 1 & y_test == 0);
    cr    = sum(y_pred == 0 & y_test == 0);

    nSignal = hits + miss;   % number of positives in test
    nNoise  = fa + cr;       % number of negatives in test

    % Compute hit/false alarm rates with log-linear (Hautus) correction
    hit_rate = (hits + 0.5) / (nSignal + 1);
    false_alarm_rate = (fa + 0.5) / (nNoise + 1);
       
    % Compute d-prime
    d_prime_scores(fold) = (norminv(hit_rate) - norminv(false_alarm_rate));
    
    disp(['Fold ', num2str(fold), ' Accuracy: ', num2str(accuracy_scores(fold) * 100), '%']);
    disp(['Fold ', num2str(fold), ' d-prime: ', num2str(d_prime_scores(fold))]);
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

save_name = [lower(behaviour), '_model_' num2str(model_number) '.mat'];
save(fullfile(models_output_folder, save_name), '-struct', 'model_struct', model_name);
disp([save_name, ' saved to ', models_output_folder]);