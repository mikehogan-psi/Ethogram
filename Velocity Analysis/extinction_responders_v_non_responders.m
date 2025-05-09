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


glme = fitglme(full_data_binned, ...
    ['Freezing ~ TimeBinStd + TrialStd + Group + ResponseType ' ...
     '+ TimeBinStd:Group + TrialStd:Group + TimeBinStd:ResponseType ' ...
     '+ ResponseType:TrialStd + ResponseType:Group:TimeBinStd + (1|Mouse)'], ...
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