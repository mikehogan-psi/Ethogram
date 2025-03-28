load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_extinction_cohort1.mat')
load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_extinction_cohort1.mat')

% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_renewal_cohort1.mat')
% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_renewal_cohort1.mat')

% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_extinction_cohort2.mat')
% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_extinction_cohort2.mat')

decay_period = 233:502;

psi_loom_data_freeze_decay = psi_loom_data(:, decay_period, :);
veh_loom_data_freeze_decay = veh_loom_data(:, decay_period, :);

% Define parameters
num_trials = 20;
num_frames = 240;
num_mice = 10;
time_bin_size = 30; % 2 seconds = 30 frames at 15 FPS
num_bins = floor(num_frames / time_bin_size);

% Preallocate arrays for binned data
freezing_psi_binned = zeros(num_bins * num_trials * num_mice, 1);
freezing_veh_binned = zeros(num_bins * num_trials * num_mice, 1);
time_bins = repmat((1:num_bins)', num_trials * num_mice, 1); % Time bin indices
trials_binned = repmat(repelem((1:num_trials)', num_bins), num_mice, 1); % Trial numbers
mice_binned = repelem((1:num_mice)', num_trials * num_bins); % Mouse IDs

% Group labels
group_psi = ones(num_bins * num_trials * num_mice, 1); % Psilocybin
group_veh = zeros(num_bins * num_trials * num_mice, 1); % Vehicle

% Process psilocybin data
idx = 1;
for mouse = 1:num_mice
    for trial = 1:num_trials
        current_trial = squeeze(psi_loom_data_freeze_decay(trial, :, mouse));
        for bin = 1:num_bins
            start_idx = (bin - 1) * time_bin_size + 1;
            end_idx = start_idx + time_bin_size - 1;
            freezing_psi_binned(idx) = mean(current_trial(start_idx:end_idx));
            idx = idx + 1;
        end
    end
end

% Process vehicle data
idx = 1;
for mouse = 1:num_mice
    for trial = 1:num_trials
        current_trial = squeeze(veh_loom_data_freeze_decay(trial, :, mouse));
        for bin = 1:num_bins
            start_idx = (bin - 1) * time_bin_size + 1;
            end_idx = start_idx + time_bin_size - 1;
            freezing_veh_binned(idx) = mean(current_trial(start_idx:end_idx));
            idx = idx + 1;
        end
    end
end

% Combine data into tables
data_psi_binned = table(freezing_psi_binned, time_bins, trials_binned, group_psi, mice_binned, ...
    'VariableNames', {'Freezing', 'TimeBin', 'Trial', 'Group', 'Mouse'});
data_veh_binned = table(freezing_veh_binned, time_bins, trials_binned, group_veh, mice_binned, ...
    'VariableNames', {'Freezing', 'TimeBin', 'Trial', 'Group', 'Mouse'});

full_data_binned = [data_psi_binned; data_veh_binned];

% Set threshold to rebinarise data (should this be here or does the glmfit
% do this again later?)
% full_data_binned.Freezing = full_data_binned.Freezing >= 0.5;

% Prepare predictors for the GLM
time_bin = full_data_binned.TimeBin;
trial = full_data_binned.Trial;
group = full_data_binned.Group;

% Compute mean and std for inverse transformation
mu_time_bin = mean(time_bin);
sigma_time_bin = std(time_bin);
mu_trial = mean(trial);
sigma_trial = std(trial);

% Interaction terms
time_bin_group = time_bin .* group;
trial_time_bin = trial .* time_bin;
trial_group = trial .* group;

% Response variable
y = full_data_binned.Freezing;

% Add intercept and predictors to matrix X
X = [ones(size(y)), time_bin, trial, group, time_bin_group, trial_time_bin, trial_group];

% Standardise predictors (Z-transform), excluding the intercept
X(:, 2:end) = zscore(X(:, 2:end));

% Fit logistic regression
[b, dev, stats] = glmfit(X(:, 2:end), y, 'binomial', 'link', 'logit');

% Display results
disp('Coefficients:');
disp(b);
disp('P-values:');
disp(stats.p);



%% Inverse transformation of coefficients
b_orig = b; % Copy original
b_orig(2) = b(2) / sigma_time_bin; % TimeBin coefficient
b_orig(3) = b(3) / sigma_trial; % Trial coefficient
b_orig(5) = b(5) / sigma_time_bin; % TimeBin*Group coefficient
b_orig(6) = b(6) / (sigma_trial * sigma_time_bin); % Trial*TimeBin coefficient
b_orig(7) = b(7) / sigma_trial; % Trial*Group coefficient

% Adjust the intercept
b_orig(1) = b(1) - (b(2) * mu_time_bin / sigma_time_bin) - (b(3) * mu_trial / sigma_trial);

% Display transformed coefficients
disp('Transformed Coefficients (Original Scale):');
disp(b_orig);

% Extract standard errors from original model
se_orig = stats.se;  

% Transform standard errors using the same scaling factors
se_transformed = se_orig;  
se_transformed(2) = se_orig(2) / sigma_time_bin; % TimeBin standard error
se_transformed(3) = se_orig(3) / sigma_trial; % Trial standard error
se_transformed(5) = se_orig(5) / sigma_time_bin; % TimeBin*Group standard error
se_transformed(6) = se_orig(6) / (sigma_trial * sigma_time_bin); % Trial*TimeBin standard error
se_transformed(7) = se_orig(7) / sigma_trial; % Trial*Group standard error

% Compute z-scores for transformed coefficients
z_transformed = b_orig ./ se_transformed;  

% Compute p-values for transformed coefficients
p_transformed = 2 * (1 - normcdf(abs(z_transformed)));  

% Display transformed p-values
disp('Transformed P-values:');
disp(p_transformed);

%% Make predictions using un-normalised data
time_vals = linspace(min(time_bin), max(time_bin), 100);
trial_vals = linspace(min(trial), max(trial), 100);

% Compute probability for time bins for both groups
X_time_veh = [ones(length(time_vals), 1), time_vals', repmat(mean(trial), length(time_vals), 1), ...
          zeros(length(time_vals), 1), time_vals' * 0, time_vals' * mean(trial), repmat(mean(trial), length(time_vals), 1) * 0];

X_time_psi = [ones(length(time_vals), 1), time_vals', repmat(mean(trial), length(time_vals), 1), ...
          ones(length(time_vals), 1), time_vals' * 1, time_vals' * mean(trial), repmat(mean(trial), length(time_vals), 1) * 1];

logit_time_veh = X_time_veh * b_orig;
logit_time_psi = X_time_psi * b_orig;

prob_time_veh = 1 ./ (1 + exp(-logit_time_veh));
prob_time_psi = 1 ./ (1 + exp(-logit_time_psi));

% Compute probability for trials for both groups
X_trial_veh = [ones(length(trial_vals), 1), repmat(mean(time_bin), length(trial_vals), 1), trial_vals', ...
           zeros(length(trial_vals), 1), repmat(mean(time_bin), length(trial_vals), 1) * 0, ...
           trial_vals' * mean(time_bin), trial_vals' * 0];

X_trial_psi = [ones(length(trial_vals), 1), repmat(mean(time_bin), length(trial_vals), 1), trial_vals', ...
           ones(length(trial_vals), 1), repmat(mean(time_bin), length(trial_vals), 1) * 1, ...
           trial_vals' * mean(time_bin), trial_vals' * 1];

logit_trial_veh = X_trial_veh * b_orig;
logit_trial_psi = X_trial_psi * b_orig;

prob_trial_veh = 1 ./ (1 + exp(-logit_trial_veh));
prob_trial_psi = 1 ./ (1 + exp(-logit_trial_psi));

%% Plot results
figure;

% Freezing probability vs. Time Bin
plot(time_vals, prob_time_veh, '-', 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Vehicle');
hold on;
plot(time_vals, prob_time_psi, '-', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Psilocybin');
xlabel('Time Bin (2s)');
ylabel('Probability of Freezing');
title('Freezing Probability vs. Time Bin (All Trials)');
legend;
xlim([1 8])
ylim([0 1])
grid on;
hold off;
print(gcf, 'freezing_probability1.png', '-dpng', '-r300');


% Freezing probability vs. Trial
figure;
plot(trial_vals, prob_trial_veh, '-', 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Vehicle');
hold on;
plot(trial_vals, prob_trial_psi, '-', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Psilocybin');
xlabel('Trial');
ylabel('Probability of Freezing');
title('Freezing Probability vs. Trial');
legend;
xlim([1 20])
ylim([0 1])
grid on;
hold off;

print(gcf, 'freezing_probability.png', '-dpng', '-r300');

%% Plot real data over models predictions: across trials
% show average probability for two groups across trials
yhat = glmval(b,X(:,2:end),'logit');
trial_val = unique(trial);
Nval = numel(trial_val);
for n = 1:Nval
    ind = find((trial == trial_val(n))&(group == 0));
    ym0(n) = mean(y(ind));
    ym_hat0(n) = mean(yhat(ind));
    ind = find((trial == trial_val(n))&(group == 1));
    ym1(n) = mean(y(ind));
    ym_hat1(n) = mean(yhat(ind));
end

figure;
hold on;

% Plot actual freezing for vehicle group
plot(trial_val, ym0, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1);

% Plot predicted freezing for vehicle group
plot(trial_val, ym_hat0, 'k--', 'LineWidth', 1.5);

% Plot actual freezing for psilocybin group
plot(trial_val, ym1, 'ro-', 'MarkerFaceColor', 'r', 'LineWidth', 1);

% Plot predicted freezing for psilocybin group
plot(trial_val, ym_hat1, 'r--', 'LineWidth', 1.5);

% Labels and title
xlabel('Trial Number', 'FontWeight', 'bold');
ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
% % title('Actual vs Predicted Freezing Across Trials');
% legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
%         'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
%        'Location', 'Best');
ylim([0 1])
xlim([1 20])
grid on;
hold off;

% Save figure
print(gcf, 'freezing_across_trials.png', '-dpng', '-r300');

%% Plot real data over models predictions: across timebins 
% show average probability for two groups across trials
yhat = glmval(b,X(:,2:end),'logit');
bin_val = unique(time_bins);
Nval = num_bins;
for n = 1:Nval
    ind = find((time_bin == bin_val(n))&(group == 0));
    ym0(n) = mean(y(ind));
    ym_hat0(n) = mean(yhat(ind));
    ind = find((time_bin == bin_val(n))&(group == 1));
    ym1(n) = mean(y(ind));
    ym_hat1(n) = mean(yhat(ind));
end

figure;
hold on;

% Plot actual freezing for vehicle group
plot(bin_val, ym0, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1);

% Plot predicted freezing for vehicle group
plot(bin_val, ym_hat0, 'k--', 'LineWidth', 1.5);

% Plot actual freezing for psilocybin group
plot(bin_val, ym1, 'ro-', 'MarkerFaceColor', 'r', 'LineWidth', 1);

% Plot predicted freezing for psilocybin group
plot(bin_val, ym_hat1, 'r--', 'LineWidth', 1.5);

% Labels and title
xlabel('Time Bin Number (2s)', 'FontWeight', 'bold');
ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
% title('Actual vs Predicted Freezing Within Trials');
legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
        'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
       'Location', 'Best');
ylim([0 1])
xlim([1 8])
grid on;
hold off;

% Save figure
print(gcf, 'freezing_across_timebins.png', '-dpng', '-r300');