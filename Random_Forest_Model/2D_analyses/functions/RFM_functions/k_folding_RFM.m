%% STEP 3: Train Random Forest Model using 5-Fold Cross-Validation

% Define model name
model_name = 'rearingA_5mice_5kfolding_without_normalization_model';

% Load all labeled data
labelled_files = dir(fullfile(output_folder, '*_labels.mat'));
all_features = [];
all_labels = [];

for i = 1:length(labelled_files)
    data = load(fullfile(output_folder, labelled_files(i).name));
    all_features = [all_features; data.features(~isnan(data.labels), :)]; % Extract valid features
    all_labels = [all_labels; data.labels(~isnan(data.labels))];
end

disp(['Total usable frames labelled: ', num2str(length(all_labels))]);

% Define number of folds
k = 5;
cv = cvpartition(length(all_labels), 'KFold', k);

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
save(['C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\extinction_analysis\rearing_data\' model_name '.mat'], 'final_model');
