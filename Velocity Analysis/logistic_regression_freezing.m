% Load data
% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_renewal_cohort1.mat')
% load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_renewal_cohort1.mat')
%load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_extinction_cohort1.mat')
%load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_extinction_cohort1.mat')
load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_extinction_cohort2.mat')
load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_extinction_cohort2.mat')

% Define decay period
decay_period = 233:502;
psi_loom_data_freeze_decay = psi_loom_data(:, decay_period, :);
veh_loom_data_freeze_decay = veh_loom_data(:, decay_period, :);

% Define parameters
num_trials = 20;
num_frames = 270;
num_mice = 10; % mice per group (change as needed)
time_bin_size = 30; % 2 seconds = 30 frames at 15 FPS
num_bins = floor(num_frames / time_bin_size);

% Preallocate arrays for binned data
freezing_psi_binned = zeros(num_bins * num_trials * num_mice, 1);
freezing_veh_binned = zeros(num_bins * num_trials * num_mice, 1);
time_bins = repmat((1:num_bins)', num_trials * num_mice, 1); % Time bin indices
trials_binned = repmat(repelem((1:num_trials)', num_bins), num_mice, 1); % Trial numbers
psi_mouse_ids = 2:2:20;  % Original even-numbered mice (change total as needed)
veh_mouse_ids = 1:2:19;  % Original odd-numbered mice (change total as needed)
mice_binned_psi = repelem(psi_mouse_ids', num_trials * num_bins);
mice_binned_veh = repelem(veh_mouse_ids', num_trials * num_bins);
group_psi = ones(num_bins * num_trials * num_mice, 1);
group_veh = zeros(num_bins * num_trials * num_mice, 1);

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

% Combine data into a single table
data_psi_binned = table(freezing_psi_binned, time_bins, trials_binned, group_psi, mice_binned_psi, ...
    'VariableNames', {'Freezing', 'TimeBin', 'Trial', 'Group', 'Mouse'});
data_veh_binned = table(freezing_veh_binned, time_bins, trials_binned, group_veh, mice_binned_veh, ...
    'VariableNames', {'Freezing', 'TimeBin', 'Trial', 'Group', 'Mouse'});

full_data_binned = [data_psi_binned; data_veh_binned];

% Binarise freezing (>= 0.5 frozen = frozen)
full_data_binned.Freezing = full_data_binned.Freezing >= 0.5;

% Standardise continuous predictors
full_data_binned.TimeBinStd = zscore(full_data_binned.TimeBin);
full_data_binned.TrialStd = zscore(full_data_binned.Trial);

% Create a categorical variable for Group (0 = veh, 1 = psi)
full_data_binned.Group = categorical(full_data_binned.Group);

% Fit Generalised Linear Mixed Model (GLMM)
glme = fitglme(full_data_binned, ...
    'Freezing ~ TimeBinStd + TrialStd + Group + TimeBinStd:Group + TrialStd:Group + (1|Mouse)', ...
    'Distribution', 'Binomial', ...
    'Link', 'Logit');

% Display model summary
disp(glme);

% % Calculate raw probability values from log odds values (if you want)
% log_odds = dataset2table(glme.Coefficients(:,2));
% log_odds = table2array(log_odds);
% convert_log_odds(log_odds);
% disp(raw_probability);

%% Plot: Actual vs Predicted Freezing Across Trials

trial_val = unique(full_data_binned.Trial);
group_val = unique(full_data_binned.Group);
Nval = numel(trial_val);

% Recalculate predicted probabilities
yhat = fitted(glme);

% Preallocate arrays
ym_actual = zeros(Nval, numel(group_val));
ym_pred = zeros(Nval, numel(group_val));
sem_actual = zeros(Nval, numel(group_val));

for g = 1:numel(group_val)
    for n = 1:Nval
        % Logical indices
        group_idx = full_data_binned.Group == group_val(g);
        trial_idx = full_data_binned.Trial == trial_val(n);
        idx = group_idx & trial_idx;

        % Get mouse IDs in this group
        mice = unique(full_data_binned.Mouse(idx));
        n_mice = numel(mice);
        per_mouse_mean = zeros(n_mice, 1);

        for m = 1:n_mice
            mouse_idx = full_data_binned.Mouse == mice(m);
            mouse_trial_idx = idx & mouse_idx;
            per_mouse_mean(m) = mean(full_data_binned.Freezing(mouse_trial_idx));
        end

        % Mean & SEM across mice
        ym_actual(n, g) = mean(per_mouse_mean);
        sem_actual(n, g) = std(per_mouse_mean) / sqrt(n_mice);

        % Predicted freezing (mean across observations)
        ym_pred(n, g) = mean(yhat(idx));
    end
end

% COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM ERROR BARS
figure; hold on;

errorbar(trial_val, ym_actual(:,1), sem_actual(:,1), 'ko-', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1);
plot(trial_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);

errorbar(trial_val, ym_actual(:,2), sem_actual(:,2), 'ro-', ...
    'MarkerFaceColor', 'r', 'LineWidth', 1);
plot(trial_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);

% (Rest of your plotting code stays the same)

xlabel('Trial Number', 'FontWeight', 'bold');
ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
        'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
       'Location', 'Best');
ylim([0 1]);
xlim([1 20]);
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
ax.LineWidth = 2;
ax.GridLineWidth = 0.5;
grid on;
hold off;

% print(gcf, 'freezing_across_trials_Acq_ren_GLMEM.png', '-dpng', '-r300');

% % COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM SHADED AREAS
% figure; hold on;
% 
% % Define colours
% vehicle_colour = [0 0 0];        % black
% psilo_colour   = [0.8 0 0];      % dark red
% 
% % Plot shaded SEM area: Vehicle
% x = trial_val;
% y1 = ym_actual(:,1) - sem_actual(:,1);
% y2 = ym_actual(:,1) + sem_actual(:,1);
% fill([x; flipud(x)], [y1; flipud(y2)], vehicle_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Plot shaded SEM area: Psilocybin
% y1 = ym_actual(:,2) - sem_actual(:,2);
% y2 = ym_actual(:,2) + sem_actual(:,2);
% fill([x; flipud(x)], [y1; flipud(y2)], psilo_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Now plot the means and predicted lines
% plot(trial_val, ym_actual(:,1), 'k-', 'LineWidth', 1.5);
% plot(trial_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);
% 
% plot(trial_val, ym_actual(:,2), 'r-', 'LineWidth', 1.5);
% plot(trial_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);
% 
% xlabel('Trial Number', 'FontWeight', 'bold');
% ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
% legend({'Vehicle SEM', 'Psilocybin SEM', ...
%         'Vehicle (Actual)', 'Vehicle (Predicted)', ...
%         'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
%        'Location', 'Best');
% ylim([0 1]);
% xlim([1 20]);
% ax = gca;
% ax.FontWeight = 'bold';
% ax.FontSize = 12;
% ax.LineWidth = 2;
% ax.GridLineWidth = 0.5;
% grid on;
% hold off;

% % print(gcf, 'freezing_across_trials_Acq_ren_GLMEM.png', '-dpng', '-r300');

%% Plot: Actual vs Predicted Freezing Across Time Bins
bin_val = unique(full_data_binned.TimeBin);
Nval = numel(bin_val);

ym_actual = zeros(Nval, numel(group_val));
ym_pred = zeros(Nval, numel(group_val));
sem_actual = zeros(Nval, numel(group_val));

for g = 1:numel(group_val)
    for n = 1:Nval
        group_idx = full_data_binned.Group == group_val(g);
        time_idx = full_data_binned.TimeBin == bin_val(n);
        idx = group_idx & time_idx;

        mice = unique(full_data_binned.Mouse(idx));
        n_mice = numel(mice);
        per_mouse_mean = zeros(n_mice, 1);

        for m = 1:n_mice
            mouse_idx = full_data_binned.Mouse == mice(m);
            mouse_time_idx = idx & mouse_idx;
            per_mouse_mean(m) = mean(full_data_binned.Freezing(mouse_time_idx));
        end

        ym_actual(n, g) = mean(per_mouse_mean);
        sem_actual(n, g) = std(per_mouse_mean) / sqrt(n_mice);
        ym_pred(n, g) = mean(yhat(idx));
    end
end

% COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM ERROR BARS
figure; hold on;

errorbar(bin_val, ym_actual(:,1), sem_actual(:,1), 'ko-', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1);
plot(bin_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);

errorbar(bin_val, ym_actual(:,2), sem_actual(:,2), 'ro-', ...
    'MarkerFaceColor', 'r', 'LineWidth', 1);
plot(bin_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);




xlabel('Time Bin Number (2s)', 'FontWeight', 'bold');
ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
        'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
       'Location', 'Best');
ylim([0 1]);
xlim([1 9]);
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
ax.LineWidth = 2;
ax.GridLineWidth = 0.5;
grid on;
hold off;

% print(gcf, 'freezing_across_timebins_acq_ren_GLMEM.png', '-dpng', '-r300');

% % COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM SHADED AREAS
% figure; hold on;
% 
% % Define colours
% vehicle_colour = [0 0 0];        % black
% psilo_colour   = [0.8 0 0];      % dark red
% 
% % Plot shaded SEM area: Vehicle
% x = bin_val;
% y1 = ym_actual(:,1) - sem_actual(:,1);
% y2 = ym_actual(:,1) + sem_actual(:,1);
% fill([x; flipud(x)], [y1; flipud(y2)], vehicle_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Plot shaded SEM area: Psilocybin
% y1 = ym_actual(:,2) - sem_actual(:,2);
% y2 = ym_actual(:,2) + sem_actual(:,2);
% fill([x; flipud(x)], [y1; flipud(y2)], psilo_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Now plot the means and predicted lines
% plot(bin_val, ym_actual(:,1), 'k-', 'LineWidth', 1.5);
% plot(bin_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);
% 
% plot(bin_val, ym_actual(:,2), 'r-', 'LineWidth', 1.5);
% plot(bin_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);
% 
% xlabel('Bin Number', 'FontWeight', 'bold');
% ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
% legend({'Vehicle SEM', 'Psilocybin SEM', ...
%         'Vehicle (Actual)', 'Vehicle (Predicted)', ...
%         'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
%        'Location', 'Best');
% ylim([0 1]);
% xlim([1 8]);
% ax = gca;
% ax.FontWeight = 'bold';
% ax.FontSize = 12;
% ax.LineWidth = 2;
% ax.GridLineWidth = 0.5;
% grid on;
% hold off;

% % print(gcf, 'freezing_across_timebins_acq_ren_GLMEM.png', '-dpng', '-r300');