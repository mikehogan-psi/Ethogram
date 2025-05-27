%% Load data
load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\psi_loom_data_extinction_cohort1.mat')
load('D:\PhD 2nd Year\MATLAB\Processed_freezing_data\veh_loom_data_extinction_cohort1.mat')

%% Calculate percentage change in freezing from trial 1 to 20

num_trials = 20;
num_mice = 15;

% Collapse frame dimension
mean_freeze_trials_psi = squeeze(mean(psi_loom_data, 2));
mean_freeze_trials_veh = squeeze(mean(veh_loom_data, 2));

% Preallocate arrays for percentage change data
percentage_change_psi = zeros(num_mice, 1);
percentage_change_veh = zeros(num_mice, 1);

% Define early v late trials here
early_trials_idx = 1:5;
late_trials_idx = num_trials-4:num_trials;

% Calculating percentage change between early and late trials for treatment
% groups
for mouse_idx = 1:num_mice
    temp_change = (sum(mean_freeze_trials_psi(early_trials_idx, mouse_idx))...
        -sum(mean_freeze_trials_psi(late_trials_idx, mouse_idx)))...
        /sum(mean_freeze_trials_psi(early_trials_idx, mouse_idx));
    percentage_change_psi(mouse_idx) = temp_change*100;
end


for mouse_idx = 1:num_mice
    temp_change = (sum(mean_freeze_trials_veh(early_trials_idx, mouse_idx))...
        -sum(mean_freeze_trials_veh(late_trials_idx, mouse_idx)))...
        /sum(mean_freeze_trials_veh(early_trials_idx, mouse_idx));
    percentage_change_veh(mouse_idx) = temp_change*100;
end

%%
figure;
scatter(ones(1, num_mice), percentage_change_psi, 60, 'ro', 'filled', 'XJitter', 'rand'); hold on;
scatter(2*ones(1, num_mice), percentage_change_veh, 60, 'ko', 'filled', 'XJitter', 'rand');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Psilocybin','Vehicle'});
ylabel('% reduction in freezing');
yline(0, lineStyle="--");
ax=gca;
ax.FontWeight = 'bold';
title('Extinction Rate By Mouse');

%% Get mouse IDs for extinction responders

% Create tables
psi_table = table((1:num_mice)', percentage_change_psi, ...
    'VariableNames', {'Mouse', 'ExtinctionRate'});

veh_table = table((1:num_mice)', percentage_change_veh, ...
    'VariableNames', {'Mouse', 'ExtinctionRate'});

extinction_responders_idx_psi = psi_table.ExtinctionRate>0;
extinction_responders_idx_veh = veh_table.ExtinctionRate>0;
extinction_non_responders_idx_psi = ~extinction_responders_idx_psi;
extinction_non_responders_idx_veh = ~extinction_responders_idx_veh;

extinction_responders_psi = psi_table.Mouse(extinction_responders_idx_psi);
extinction_responders_veh = veh_table.Mouse(extinction_responders_idx_veh);
extinction_non_responders_psi = psi_table.Mouse(extinction_non_responders_idx_psi);
extinction_non_responders_veh = veh_table.Mouse(extinction_non_responders_idx_veh);

%% Create table with all data and perform logistic regression

% Create tables
psi_table_temp = table((2:2:num_mice*2)', percentage_change_psi, ...
    'VariableNames', {'Mouse', 'ExtinctionRate'});

veh_table_temp = table((1:2:num_mice*2)', percentage_change_veh, ...
    'VariableNames', {'Mouse', 'ExtinctionRate'});

all_data_table = vertcat(psi_table_temp, veh_table_temp);
treatment_type = [zeros(15, 1); ones(15, 1)];
all_data_table.Treatment = treatment_type;

extinction_responders_idx = all_data_table.ExtinctionRate>0;
extinction_non_responders_idx = ~extinction_responders_idx;
extinction_responders = all_data_table.Mouse(extinction_responders_idx);
extinction_non_responders = all_data_table.Mouse(extinction_non_responders_idx);
responders_double = double(extinction_responders_idx);
all_data_table.ResponseType = responders_double;

load("D:\PhD 2nd Year\full_data_binned.mat")

num_trials = 20;
num_bins = 9;

responders_double_binned = repelem(responders_double, num_trials*num_bins);
full_data_binned.ResponseType = responders_double_binned;


% glme = fitglme(full_data_binned, ...
%     ['Freezing ~ TimeBinStd + TrialStd + Group + ResponseType ' ...
%      '+ TimeBinStd:Group + TrialStd:Group + TimeBinStd:ResponseType ' ...
%      '+ ResponseType:TrialStd + (1|Mouse)'], ...
%     'Distribution', 'Binomial', ...
%     'Link', 'Logit');

glme = fitglme(full_data_binned, ...
    ['Freezing ~ TimeBinStd + TrialStd + Group + ResponseType ' ...
     '+ TimeBinStd:Group + TrialStd:Group + TimeBinStd:ResponseType ' ...
     '+ ResponseType:TrialStd + Group:ResponseType ' ...
     '+ TimeBinStd:Group:ResponseType + TrialStd:Group:ResponseType ' ...
     '+ (1|Mouse)'], ...
    'Distribution', 'Binomial', ...
    'Link', 'Logit');


% Display model summary
disp(glme);

%%
% Step 1: Collapse frame dimension (mean freezing per trial, per mouse)
mean_freeze_trials_psi = squeeze(mean(psi_loom_data, 2));  % [20 x 15]
mean_freeze_trials_veh = squeeze(mean(veh_loom_data, 2));  % [20 x 15]

% Step 2: Extract responders and non-responders by mouse ID
mean_psi_responders = mean(mean_freeze_trials_psi(:, extinction_responders_psi),2);
mean_veh_responders = mean(mean_freeze_trials_veh(:, extinction_responders_veh),2);
mean_psi_nonresponders = mean(mean_freeze_trials_psi(:, extinction_non_responders_psi),2);
mean_veh_nonresponders = mean(mean_freeze_trials_veh(:, extinction_non_responders_veh),2);


% Step 3: Plot responders
figure;
plot(1:20, mean_psi_responders, '-ro', 'LineWidth', 2); hold on;
plot(1:20, mean_veh_responders, '-ko', 'LineWidth', 2);
title('Extinction Responders'); xlabel('Trial'); ylabel('Mean Freezing');
legend('Psilocybin', 'Vehicle');
xlim([1 20])
ylim([0 1])
grid on
ax=gca;
ax.FontWeight = 'bold';

% Step 4: Plot non-responders
figure;
plot(1:20, mean_psi_nonresponders, '-ro', 'LineWidth', 2); hold on;
plot(1:20, mean_veh_nonresponders, '-ko', 'LineWidth', 2);
title('Extinction Non-Responders'); xlabel('Trial'); ylabel('Mean Freezing');
legend('Psilocybin', 'Vehicle');
xlim([1 20])
ylim([0 1])
grid on
ax=gca;
ax.FontWeight = 'bold';


%% Plotting all data together (no grouping of responders)
% Collapse frame dimension: mean freezing per trial per mouse
mean_freeze_trials_psi = squeeze(mean(psi_loom_data, 2));  % [20 x 15]
mean_freeze_trials_veh = squeeze(mean(veh_loom_data, 2));  % [20 x 15]

% Average across mice
mean_freeze_psi = mean(mean_freeze_trials_psi, 2);  % [20 x 1]
mean_freeze_veh = mean(mean_freeze_trials_veh, 2);  % [20 x 1]

% Plot both groups
figure;
plot(1:20, mean_freeze_psi, '-ro', 'LineWidth', 2); hold on;
plot(1:20, mean_freeze_veh, '-ko', 'LineWidth', 2);
title('Extinction Curve: All Mice'); xlabel('Trial'); ylabel('Mean Freezing');
legend('Psilocybin', 'Vehicle');
xlim([1 20]);
ylim([0 1]);
grid on;
%% Calculate slope of responders using linear model
% Create trial vector
trials = (1:20)';

% Linear fit for Psilocybin
lm_psi = fitlm(trials, mean_psi_responders);
slope_psi = lm_psi.Coefficients.Estimate(2);  % slope (rate of change)
r2_psi = lm_psi.Rsquared.Ordinary;           % R² value

% Linear fit for Vehicle
lm_veh = fitlm(trials, mean_veh_responders);
slope_veh = lm_veh.Coefficients.Estimate(2);
r2_veh = lm_veh.Rsquared.Ordinary;

% Display results
fprintf('Psilocybin slope: %.4f, R²: %.4f\n', slope_psi, r2_psi);
fprintf('Vehicle slope: %.4f, R²: %.4f\n', slope_veh, r2_veh);

%% Plotting overall freezing of non-learners v learners

% --- 1) Compute mean freezing per mouse (across all trials & frames) ---
% psi_loom_data(trial, frame, mouse)
mean_freeze_per_mouse_psi = squeeze(mean(mean(psi_loom_data, 1), 2));   % [15 x 1]
mean_freeze_per_mouse_veh = squeeze(mean(mean(veh_loom_data, 1), 2));   % [15 x 1]

% --- 2) Split into responders vs non-responders ---
psi_resp_vals    = mean_freeze_per_mouse_psi(extinction_responders_psi);
psi_nonresp_vals = mean_freeze_per_mouse_psi(extinction_non_responders_psi);
veh_resp_vals    = mean_freeze_per_mouse_veh(extinction_responders_veh);
veh_nonresp_vals = mean_freeze_per_mouse_veh(extinction_non_responders_veh);

% --- 3) Prepare data for a grouped bar chart ---
groupNames = {'Learners','Non‑Learners'};
psi_means    = [ mean(psi_resp_vals),    mean(psi_nonresp_vals)   ];
psi_sems     = [ std(psi_resp_vals)/sqrt(numel(psi_resp_vals)), ...
                 std(psi_nonresp_vals)/sqrt(numel(psi_nonresp_vals)) ];
veh_means    = [ mean(veh_resp_vals),    mean(veh_nonresp_vals)   ];
veh_sems     = [ std(veh_resp_vals)/sqrt(numel(veh_resp_vals)), ...
                 std(veh_nonresp_vals)/sqrt(numel(veh_nonresp_vals)) ];

dataMeans = [ psi_means; veh_means ];    % 2×2 matrix: rows = [Psilo; Veh], cols = [Resp; Non‑resp]
dataSEMs  = [ psi_sems;  veh_sems ];

% --- 4) Plot grouped bar with error bars ---
figure;
hb = bar(dataMeans, 'grouped');
hold on;

% number of groups and number of bars per group:
[ngroups, nbars] = size(dataMeans);
% Get the x-coordinate of each bar group
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = hb(i).XEndPoints;
end

% Plot error bars
for i = 1:nbars
    errorbar(x(i,:), dataMeans(:,i)', dataSEMs(:,i)', 'k', 'linestyle','none', 'LineWidth',1);
end

% Labels & styling
set(gca, 'XTickLabel', {'Psilocybin','Vehicle'});
legend(groupNames,'Location','Best');
ylabel('Mean freezing (all extinction trials)');
ylim([0 1]);
title('Overall freezing: learners vs non‑learners');
grid on;
ax=gca;
ax.FontWeight = 'bold';
hold off;

%% Plotting overall freezing of non-learners v learners

% --- 1) Compute mean freezing per mouse (across all trials & frames) ---
% psi_loom_data(trial, frame, mouse)
mean_freeze_per_mouse_psi_first5 = squeeze(mean(mean(psi_loom_data(1:5, :, :), 1), 2));   % [15 x 1]
mean_freeze_per_mouse_veh_first5 = squeeze(mean(mean(veh_loom_data(1:5, : , :), 1), 2));   % [15 x 1]

% --- 2) Split into responders vs non-responders ---
psi_resp_vals    = mean_freeze_per_mouse_psi_first5(extinction_responders_psi);
psi_nonresp_vals = mean_freeze_per_mouse_psi_first5(extinction_non_responders_psi);
veh_resp_vals    = mean_freeze_per_mouse_veh_first5(extinction_responders_veh);
veh_nonresp_vals = mean_freeze_per_mouse_veh_first5(extinction_non_responders_veh);

% --- 3) Prepare data for a grouped bar chart ---
groupNames = {'Learners','Non‑Learners'};
psi_means    = [ mean(psi_resp_vals),    mean(psi_nonresp_vals)   ];
psi_sems     = [ std(psi_resp_vals)/sqrt(numel(psi_resp_vals)), ...
                 std(psi_nonresp_vals)/sqrt(numel(psi_nonresp_vals)) ];
veh_means    = [ mean(veh_resp_vals),    mean(veh_nonresp_vals)   ];
veh_sems     = [ std(veh_resp_vals)/sqrt(numel(veh_resp_vals)), ...
                 std(veh_nonresp_vals)/sqrt(numel(veh_nonresp_vals)) ];

dataMeans = [ psi_means; veh_means ];    % 2×2 matrix: rows = [Psilo; Veh], cols = [Resp; Non‑resp]
dataSEMs  = [ psi_sems;  veh_sems ];

% --- 4) Plot grouped bar with error bars ---
figure;
hb = bar(dataMeans, 'grouped');
hold on;

% number of groups and number of bars per group:
[ngroups, nbars] = size(dataMeans);
% Get the x-coordinate of each bar group
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = hb(i).XEndPoints;
end

% Plot error bars
for i = 1:nbars
    errorbar(x(i,:), dataMeans(:,i)', dataSEMs(:,i)', 'k', 'linestyle','none', 'LineWidth',1);
end

% Labels & styling
set(gca, 'XTickLabel', {'Psilocybin','Vehicle'});
legend(groupNames,'Location','Best');
ylabel('Mean freezing (all extinction trials)');
ylim([0 1]);
title('Overall freezing: learners vs non‑learners (First 5 Trials Only)');
grid on;
ax=gca;
ax.FontWeight = 'bold';
hold off;

%% Plot mean across frames for responders v non-responders
% Step 1: Collapse trial dimension 
decay_period = 233:502;

mean_freeze_frames_psi = squeeze(mean(psi_loom_data(:, decay_period, :), 1));  
mean_freeze_frames_veh = squeeze(mean(veh_loom_data(:, decay_period, :), 1));  

% Step 2: Extract responders and non-responders by mouse ID
mean_psi_responders_frames = mean(mean_freeze_frames_psi(:, extinction_responders_psi),2);
mean_veh_responders_frames = mean(mean_freeze_frames_veh(:, extinction_responders_veh),2);
mean_psi_nonresponders_frames = mean(mean_freeze_frames_psi(:, extinction_non_responders_psi),2);
mean_veh_nonresponders_frames = mean(mean_freeze_frames_veh(:, extinction_non_responders_veh),2);

% Define time axis (frames 233 to 502)
time_axis = decay_period;  % Should be 1x270


% Define time axis
time_axis = decay_period;

% Create figure
figure;

% --- Plot responders ---
subplot(1,2,1);
plot(time_axis, mean_psi_responders_frames, 'r', 'LineWidth', 1.5); hold on;
plot(time_axis, mean_veh_responders_frames, 'k', 'LineWidth', 1.5);
xlabel('Frame');
ylabel('Mean Freezing (fraction of time)');
title('Responders');
legend({'Psilocybin', 'Vehicle'}, 'Location', 'best');
grid on;
ylim([0 1]);  % Adjust based on your data range
xlim([233 502])
ax=gca;
ax.FontWeight = 'bold';

% --- Plot non-responders ---
subplot(1,2,2);
plot(time_axis, mean_psi_nonresponders_frames, 'r', 'LineWidth', 1.5); hold on;
plot(time_axis, mean_veh_nonresponders_frames, 'k', 'LineWidth', 1.5);
xlabel('Frame');
ylabel('Mean Freezing (fraction of time)');
title('Non-responders');
legend({'Psilocybin', 'Vehicle'}, 'Location', 'best');
grid on;
ylim([0 1]);  % Match Y-axis for comparison
xlim([233 502])
ax=gca;
ax.FontWeight = 'bold';

%% Plot predicted v actual data for responder v non-responders across trials
trial_val = unique(full_data_binned.Trial);
group_val = unique(full_data_binned.Group);
response_val = unique(full_data_binned.ResponseType);  % [0 1]
Nval = numel(trial_val);

yhat = fitted(glme);

for r = 1:numel(response_val)

    % Preallocate for each response type
    ym_actual = zeros(Nval, numel(group_val));
    ym_pred = zeros(Nval, numel(group_val));
    sem_actual = zeros(Nval, numel(group_val));

    response_idx = full_data_binned.ResponseType == response_val(r);

    for g = 1:numel(group_val)
        for n = 1:Nval
            % Logical indices
            group_idx = full_data_binned.Group == group_val(g);
            trial_idx = full_data_binned.Trial == trial_val(n);
            idx = group_idx & trial_idx & response_idx;

            mice = unique(full_data_binned.Mouse(idx));
            n_mice = numel(mice);
            per_mouse_mean = zeros(n_mice, 1);

            for m = 1:n_mice
                mouse_idx = full_data_binned.Mouse == mice(m);
                mouse_trial_idx = idx & mouse_idx;
                per_mouse_mean(m) = mean(full_data_binned.Freezing(mouse_trial_idx));
            end

            ym_actual(n, g) = mean(per_mouse_mean);
            sem_actual(n, g) = std(per_mouse_mean) / sqrt(n_mice);
            ym_pred(n, g) = mean(yhat(idx));
        end
    end

    % Plot
    figure; hold on;

    errorbar(trial_val, ym_actual(:,1), sem_actual(:,1), 'ko-', ...
        'MarkerFaceColor', 'k', 'LineWidth', 1);
    plot(trial_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);

    errorbar(trial_val, ym_actual(:,2), sem_actual(:,2), 'ro-', ...
        'MarkerFaceColor', 'r', 'LineWidth', 1);
    plot(trial_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);

    xlabel('Trial Number', 'FontWeight', 'bold');
    ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
    legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
            'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
           'Location', 'Best');

    ylim([0 1]);
    xlim([1 20]);

    title(sprintf('Freezing Across Trials: %s', ...
        ternary(response_val(r)==1, 'Responders', 'Non-Responders')));

    ax = gca;
    ax.FontWeight = 'bold';
    ax.FontSize = 12;
    ax.LineWidth = 2;
    ax.GridLineWidth = 0.5;
    grid on;
    hold off;

end
%% Plot predicted v actual data for responder v non-responders across timebins
bin_val = unique(full_data_binned.TimeBin);
group_val = unique(full_data_binned.Group);
response_val = unique(full_data_binned.ResponseType);  % [0 1]
Nval = numel(bin_val);

yhat = fitted(glme);

for r = 1:numel(response_val)

    % Preallocate for each response type
    ym_actual = zeros(Nval, numel(group_val));
    ym_pred = zeros(Nval, numel(group_val));
    sem_actual = zeros(Nval, numel(group_val));

    response_idx = full_data_binned.ResponseType == response_val(r);

    for g = 1:numel(group_val)
        for n = 1:Nval
            group_idx = full_data_binned.Group == group_val(g);
            time_idx = full_data_binned.TimeBin == bin_val(n);
            idx = group_idx & time_idx & response_idx;

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

    % Plot
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

    title(sprintf('Freezing Across Time Bins: %s', ...
        ternary(response_val(r)==1, 'Responders', 'Non-Responders')));

    ylim([0 1]);
    xlim([min(bin_val), max(bin_val)]);
    ax = gca;
    ax.FontWeight = 'bold';
    ax.FontSize = 12;
    ax.LineWidth = 2;
    ax.GridLineWidth = 0.5;
    grid on;
    hold off;

end
%%
% Extract coefficient names and values
% Extract coefficient names and values
coef_names = glme.CoefficientNames;
estimates = glme.Coefficients.Estimate;

% Define helper to retrieve coefficient (returns 0 if not found)
get_coef = @(name) estimates(strcmp(coef_names, name));

% Main effect and interactions
b_base  = get_coef('TimeBinStd');                        % Vehicle + Non-responder
b_grp   = get_coef('Group_1:TimeBinStd');                % Psilocybin vs Vehicle
b_resp  = get_coef('TimeBinStd:ResponseType');           % Responder vs Non-responder
b_3way  = get_coef('Group_1:TimeBinStd:ResponseType');   % Interaction of all 3

% Compute simple slopes
slope_vehicle_nonresponder   = b_base;
slope_psilo_nonresponder     = b_base + b_grp;
slope_vehicle_responder      = b_base + b_resp;
slope_psilo_responder        = b_base + b_grp + b_resp + b_3way;

% Display the results
fprintf('\nSimple slopes for TimeBinStd (Freezing decay rate):\n');
fprintf('Vehicle + Non-responder:    %.4f\n', slope_vehicle_nonresponder);
fprintf('Psilocybin + Non-responder: %.4f\n', slope_psilo_nonresponder);
fprintf('Vehicle + Responder:        %.4f\n', slope_vehicle_responder);
fprintf('Psilocybin + Responder:     %.4f\n', slope_psilo_responder);

%%
% Coefficients and covariance matrix
coef_names = glme.CoefficientNames;
estimates = glme.Coefficients.Estimate;
cov_mat = glme.CoefficientCovariance;

% Helper to get index in vector/matrix
idx = @(name) find(strcmp(coef_names, name));

% Define contrasts (weights for each term in the slope)
% These vectors must match the length of coef_names
n_coef = numel(coef_names);
L_vehicle_nonresponder   = zeros(n_coef,1); L_vehicle_nonresponder(idx('TimeBinStd')) = 1;
L_psilo_nonresponder     = L_vehicle_nonresponder; L_psilo_nonresponder(idx('Group_1:TimeBinStd')) = 1;
L_vehicle_responder      = L_vehicle_nonresponder; L_vehicle_responder(idx('TimeBinStd:ResponseType')) = 1;
L_psilo_responder        = L_psilo_nonresponder; L_psilo_responder(idx('TimeBinStd:ResponseType')) = 1;
L_psilo_responder(idx('Group_1:TimeBinStd:ResponseType')) = 1;

% Pack into cell array for looping
L_all = {L_vehicle_nonresponder, L_psilo_nonresponder, L_vehicle_responder, L_psilo_responder};
labels = {'Vehicle + Non-responder', 'Psilocybin + Non-responder', ...
          'Vehicle + Responder', 'Psilocybin + Responder'};

% Degrees of freedom
df = glme.DFE;

% Compute slopes, SEs, t-stats, and p-values
fprintf('\nTest of simple slopes (vs zero):\n');
for i = 1:4
    L = L_all{i};
    slope = L' * estimates;
    SE = sqrt(L' * cov_mat * L);
    t = slope / SE;
    p = 2 * (1 - tcdf(abs(t), df));
    
    fprintf('%s:\n', labels{i});
    fprintf('  Slope = %.4f, SE = %.4f, t = %.2f, p = %.4f\n', slope, SE, t, p);
end

%%
% Contrast for Vehicle Responder - Vehicle Non-responder
L_diff = L_psilo_responder - L_psilo_nonresponder;

slope_diff = L_diff' * estimates;
SE_diff = sqrt(L_diff' * cov_mat * L_diff);
t_diff = slope_diff / SE_diff;
p_diff = 2 * (1 - tcdf(abs(t_diff), df));

fprintf('\nDifference in slope (PSilocybin Responder vs Psilocybin Non-responder):\n');
fprintf('  Diff = %.4f, SE = %.4f, t = %.2f, p = %.4f\n', slope_diff, SE_diff, t_diff, p_diff);

%%
function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end
