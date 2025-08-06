% Get indices of loom trials (stim_set_1 == 1)
%0 = flash, 1 = loom
stim_set_1 = [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1,...
    0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1];

post_stim_period = 202:502;

stim_set_1_looms_idx = logical(stim_set_1);

% Extract only loom trials from the freezing data
loom_trials = validated_freeze_matrix(stim_set_1_looms_idx, :);  % size: n_looms x 502

% Compute mean freezing across loom trials at each frame
mean_freezing_loom = mean(loom_trials, 1);  % 1 x 502

% Plot
figure;
plot(mean_freezing_loom, 'LineWidth', 2);
xlabel('Frame');
ylabel('Proportion Freezing');
title('Mean Freezing Across Loom Trials');
xlim([1 502]);
ylim([0 1]);
grid on;

%%
% Logical index for loom trials
stim_set_1_looms_idx = logical(stim_set_1);

% Extract loom trials from binary freezing matrix (40 x 502)
loom_trials = validated_freeze_matrix(stim_set_1_looms_idx, post_stim_period);  % (n_looms x 502)

% Compute mean freezing per trial (mean over frames)
mean_freezing_per_trial = mean(loom_trials, 2);  % (n_looms x 1)

% Plot
figure;
plot(mean_freezing_per_trial, 'o-', 'LineWidth', 2);
xlabel('Loom Trial #');
ylabel('Mean Freezing (Proportion of 502 frames)');
title('Mean Freezing Per Loom Trial');
ylim([0 1]);
xlim([0 20])
grid on;
