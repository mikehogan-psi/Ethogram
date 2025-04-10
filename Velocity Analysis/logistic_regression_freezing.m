% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_extinction_cohort1.mat')
% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_extinction_cohort1.mat')

load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_renewal_cohort1.mat')
load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_renewal_cohort1.mat')

% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_extinction_cohort2.mat')
% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_extinction_cohort2.mat')

% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_renewal_cohort2.mat')
% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_renewal_cohort2.mat')

% Define decay period

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

% Apply binarisation for logistic regression
full_data_binned.Freezing = full_data_binned.Freezing >= 0.5;

% Prepare predictors for the GLM
time_bin = full_data_binned.TimeBin;
trial = full_data_binned.Trial;
group = full_data_binned.Group;

% Standardise continuous predictors only
mu_time_bin = mean(time_bin);
sigma_time_bin = std(time_bin);
mu_trial = mean(trial);
sigma_trial = std(trial);
time_bin_standardised = (time_bin - mu_time_bin) / sigma_time_bin;
trial_standardised = (trial - mu_trial) / sigma_trial;

% Interaction terms
time_bin_group = time_bin .* group;
trial_group = trial .* group;

% Response variable
y = full_data_binned.Freezing;

% Add intercept and predictors to matrix X
X = [ones(size(y)), time_bin_standardised, trial_standardised, group, time_bin_group, trial_group];

% Fit logistic regression
[b, dev, stats] = glmfit(X(:, 2:end), y, 'binomial', 'link', 'logit');

% Define labels for each coefficient
labels = { 'Intercept', ...
           'Time Bin', ...
           'Trial', ...
           'Group (Psi vs Veh)', ...
           'Time Bin × Group', ...
           'Trial × Group' };

% Display results with labels
disp('**COEFFICIENTS**:');
for i = 1:length(b)
    fprintf('%s: %.6f\n', labels{i}, b(i));
end

disp('**P-VALUES**:');
for i = 1:length(stats.p)
    fprintf('%s: %.6f\n', labels{i}, stats.p(i));
end


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
% title('Actual vs Predicted Freezing Across Trials');
legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
        'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
       'Location', 'Best');
ylim([0 1])
xlim([1 20])
% Set axis properties
ax = gca;
ax.FontWeight = 'bold';  % Make all text in axes bold
ax.FontSize = 12;         % Increase font size for better visibility
ax.LineWidth = 2;         % Make the axes lines bold
ax.GridLineWidth = 0.5;   % Keep grid lines thin
grid on;
hold off;

% Save figure
print(gcf, 'freezing_across_trials_Acq_ren.png', '-dpng', '-r300');

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
% Set axis properties
ax = gca;
ax.FontWeight = 'bold';  % Make all text in axes bold
ax.FontSize = 12;         % Increase font size for better visibility
ax.LineWidth = 2;         % Make the axes lines bold
ax.GridLineWidth = 0.5;   % Keep grid lines thin
grid on;
hold off;

% save_folder = 'D:\PhD 2nd Year\Figures'; 
% save_path = fullfile(save_folder, 'freezing_across_timebins_preext.png');
% % Save figure
% print(gcf, '-dpng', '-r300');

 print(gcf, 'freezing_across_timebins_acq_ren.png', '-dpng', '-r300');