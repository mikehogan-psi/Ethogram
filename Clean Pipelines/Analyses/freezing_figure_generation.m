% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!!Provide session to be analysed!!!
session = 'Renewal';

trial_type = 'looms'; % 'looms' or 'flashes'

% !!! Provide treatment mouse numbers for each treatment group !!!
received_psilocybin = [3; 5; 7; 8];
received_vehicle = [1; 2; 4; 6; 9];
mice_to_analyse = [1 2 3 4 5 6 7 8 9];
%%
% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

% Initialise 
freezing_folders = cell(length(mouse_files), 1);

% Get DLC folders for specified session
for mouse = mice_to_analyse
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    freezing_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Freezing');
end

freeze_data = cell(length(mouse_files), 1);

for mouse = mice_to_analyse
    current_freeze_data_path = freezing_folders{mouse};
    current_freeze_data = dir(fullfile(current_freeze_data_path, ['*_', trial_type, '_freezing.mat']));
    file_to_load = fullfile(current_freeze_data.folder, current_freeze_data.name);
    if strcmp(trial_type, 'looms')
        load(file_to_load, 'loom_freezing');
        freeze_data{mouse} = loom_freezing;
    elseif strcmp(trial_type, 'flashes')
        load(file_to_load, 'flash_freezing');
        freeze_data{mouse} = flash_freezing;
    end
end

freeze_data_psi = freeze_data(received_psilocybin);
freeze_data_veh = freeze_data(received_vehicle);

freeze_data_psi_cat = cat(3, freeze_data_psi{:});
freeze_data_veh_cat = cat(3, freeze_data_veh{:});

%%
% assume these exist and are trials x frames x mice
psi = freeze_data_psi_cat;
veh = freeze_data_veh_cat;

% basic checks
if ndims(psi)~=3 || ndims(veh)~=3
    error('Expect data to be 3-D arrays: [trials x frames x mice].');
end
[n_trials1, n_frames1, n_mice_psi] = size(psi);
[n_trials2, n_frames2, n_mice_veh] = size(veh);
if n_trials1~=n_trials2 || n_frames1~=n_frames2
    error('Psi and Veh arrays must have same trials and frames dims.');
end
n_trials = n_trials1;
n_frames = n_frames1;

% convert logical to numeric if needed
psi = double(psi);
veh = double(veh);

%% Trial-wise: average over frames, then compute mean +/- SEM across mice
% per-mouse trial averages: trials x mice
psi_trial_per_mouse = squeeze(mean(psi, 2));   % [trials x mice]
veh_trial_per_mouse = squeeze(mean(veh, 2));   % [trials x mice]

% mean across mice (trials x 1) and SEM across mice
psi_trial_mean = mean(psi_trial_per_mouse, 2);                 % trials x 1
psi_trial_sem  = std(psi_trial_per_mouse, 0, 2) ./ sqrt(n_mice_psi);

veh_trial_mean = mean(veh_trial_per_mouse, 2);
veh_trial_sem  = std(veh_trial_per_mouse, 0, 2) ./ sqrt(n_mice_veh);

% Plot trial-wise (as percent)
figure;
errorbar(1:n_trials, psi_trial_mean*100, psi_trial_sem*100, 'r-o', 'LineWidth', 1.5); hold on;
errorbar(1:n_trials, veh_trial_mean*100,  veh_trial_sem*100,  'b-o', 'LineWidth', 1.5);
xlabel('Trial'); ylabel('Mean freezing (%)'); title([session, ' ', trial_type, ' mean freezing across trials']);
legend('Psilocybin (± SEM)','Vehicle (± SEM)','Location','best');
grid on;
ylim([0 100]);
xlim([1 20])

%% Frame-wise: average over trials, then compute mean +/- SEM across mice
% per-mouse frame averages: frames x mice
psi_frame_per_mouse = squeeze(mean(psi, 1));   % [frames x mice]
veh_frame_per_mouse = squeeze(mean(veh, 1));   % [frames x mice]

% mean across mice (frames x 1) and SEM across mice
psi_frame_mean = mean(psi_frame_per_mouse, 2);                % frames x 1
psi_frame_sem  = std(psi_frame_per_mouse, 0, 2) ./ sqrt(n_mice_psi);

veh_frame_mean = mean(veh_frame_per_mouse, 2);
veh_frame_sem  = std(veh_frame_per_mouse, 0, 2) ./ sqrt(n_mice_veh);

% Plot frame-wise with shaded SEM
x = (1:n_frames)';

figure; hold on;
% psilocybin shading
psi_up = psi_frame_mean + psi_frame_sem;
psi_lo = psi_frame_mean - psi_frame_sem;
h_fill1 = fill([x; flipud(x)], [psi_up; flipud(psi_lo)]*100, [1 .7 .7], 'LineStyle','none'); 
set(h_fill1,'FaceAlpha',0.3);
plot(x, psi_frame_mean*100, 'r-', 'LineWidth', 1.5);
xlim([1 502])
ylim([0 100])

% vehicle shading
veh_up = veh_frame_mean + veh_frame_sem;
veh_lo = veh_frame_mean - veh_frame_sem;
h_fill2 = fill([x; flipud(x)], [veh_up; flipud(veh_lo)]*100, [0.7 0.7 1], 'LineStyle','none');
set(h_fill2,'FaceAlpha',0.3);
plot(x, veh_frame_mean*100, 'b-', 'LineWidth', 1.5);

xlabel('Frame'); ylabel('Mean freezing (%)'); title([session, ' ', trial_type, ' mean freezing across frames']);
legend('Psilocybin SEM','Psilocybin','Vehicle SEM','Vehicle'); grid on;

%% Heatmap: mean across mice (trials x frames)
psi_mean_over_mice = mean(psi, 3);  % trials x frames
veh_mean_over_mice = mean(veh, 3);  % trials x frames

figure;
subplot(1,2,1);
imagesc(psi_mean_over_mice); axis xy;
colorbar; xlabel('Frame'); ylabel('Trial'); title('Psilocybin: mean freezing');

subplot(1,2,2);
imagesc(veh_mean_over_mice); axis xy;
colorbar; xlabel('Frame'); ylabel('Trial'); title('Vehicle: mean freezing');
cb2 = colorbar; 
ylabel(cb2, 'Freezing Probability');
sgtitle([session, ' ', trial_type])
% Optional: smoothing across frames (uncomment to use)
% smooth_window = 5;
% psi_frame_mean = movmean(psi_frame_mean, smooth_window);
% veh_frame_mean = movmean(veh_frame_mean, smooth_window);

%%
% Calculate per-mouse freezing means
post_start = 233;
post_end   = size(psi, 2);

% psilocybin group
psi_post_per_mouse = squeeze(mean(mean(psi(:, post_start:post_end, :), 2), 1));

% vehicle group
veh_post_per_mouse = squeeze(mean(mean(veh(:, post_start:post_end, :), 2), 1));

[p,~,stats] = ranksum(psi_post_per_mouse, veh_post_per_mouse);
fprintf('Post-stim Mann–Whitney U test: p = %.4f\n', p);

% Convert to %
veh_vals = veh_post_per_mouse(:) * 100;
psi_vals = psi_post_per_mouse(:) * 100;

means = [mean(veh_vals), mean(psi_vals)];
sems  = [std(veh_vals)/sqrt(numel(veh_vals)), std(psi_vals)/sqrt(numel(psi_vals))];

figure; hold on;

% Bars
b = bar(1:2, means, 'FaceColor', 'flat');

% Vehicle (1) = blue, Psilocybin (2) = red
b.CData(1,:) = [0 0.4470 0.7410];
b.CData(2,:) = [1 0 0];

% Thick black outline on bars
b.EdgeColor = 'k';
b.LineWidth = 1.5;

% Error bars (black)
errorbar(1:2, means, sems, 'k.', 'LineWidth', 1.5);

% --- Individual datapoints (jittered) ---
jit = 0.10; % jitter amount (adjust 0.05–0.2)

xVeh = 1 + (rand(size(veh_vals)) - 0.5) * 2*jit;
xPsi = 2 + (rand(size(psi_vals)) - 0.5) * 2*jit;

scatter(xVeh, veh_vals, 35, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'k');
scatter(xPsi, psi_vals, 35, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'k');

% Cosmetics
ylim([0 100])
set(gca, 'XTick', 1:2, 'XTickLabel', {'Vehicle','Psilocybin'});
ylabel('Post-stim freezing (%)');
title('Mean post-stimulus freezing (all trials)');
box off;

ax = gca;
ax.LineWidth = 1.5;
ax.XColor    = 'k';
ax.YColor    = 'k';
ax.TickDir   = 'out';
ax.Box       = 'off';



%%
% Compute U statistic
n1 = length(psi_post_per_mouse);
n2 = length(veh_post_per_mouse);

U = stats.ranksum - n1*(n1+1)/2;   % convert ranksum to U

r_rb = 1 - (2*U)/(n1*n2);
fprintf("Rank-biserial correlation: %.3f\n", r_rb);


pooled_sd = sqrt( (var(psi_post_per_mouse)*(n1-1) + var(veh_post_per_mouse)*(n2-1)) / (n1+n2-2) );
d_cohen = (mean(psi_post_per_mouse) - mean(veh_post_per_mouse)) / pooled_sd;

fprintf("Cohen's d: %.3f\n", d_cohen);

%% Normality tests for post-stim freezing (built-in MATLAB only)

psi_vals = psi_post_per_mouse(:) * 100;   % convert to %
veh_vals = veh_post_per_mouse(:) * 100;

% Lilliefors normality tests
[psi_h, psi_p] = lillietest(psi_vals);
[veh_h, veh_p] = lillietest(veh_vals);

fprintf('\nNormality tests %s (Lilliefors):\n', session);
fprintf('  Psilocybin: p = %.4f (%s)\n', psi_p, ternary(psi_h==0, 'normal', 'non-normal'));
fprintf('  Vehicle:    p = %.4f (%s)\n', veh_p, ternary(veh_h==0, 'normal', 'non-normal'));


%% Unpaired t-test for post-stim freezing

% Convert to % for readability (optional)
psi_vals = psi_post_per_mouse(:) * 100;
veh_vals = veh_post_per_mouse(:) * 100;

% Normality tests with Lilliefors
[psi_h, psi_p] = lillietest(psi_vals);
[veh_h, veh_p] = lillietest(veh_vals);

fprintf('\nNormality tests Extinction (Lilliefors):\n');
fprintf('  Psilocybin: p = %.4f (%s)\n', psi_p, ternary(psi_h==0,'normal','non-normal'));
fprintf('  Vehicle:    p = %.4f (%s)\n', veh_p, ternary(veh_h==0,'normal','non-normal'));

% Unpaired t-test (Welch-corrected by default)
[~, p_ttest, ~, stats_ttest] = ttest2(psi_vals, veh_vals);

fprintf('\nUnpaired t-test (post-stim freezing):\n');
fprintf('  t(%d) = %.3f, p = %.4f\n', stats_ttest.df, stats_ttest.tstat, p_ttest);

% Cohen's d (independent groups)
mean_diff = mean(psi_vals) - mean(veh_vals);
pooled_sd = sqrt( ...
    ((numel(psi_vals)-1)*var(psi_vals) + (numel(veh_vals)-1)*var(veh_vals)) ...
    / (numel(psi_vals) + numel(veh_vals) - 2));

d_ttest = mean_diff / pooled_sd;

fprintf('  Cohen''s d (t-test) = %.3f\n', d_ttest);

%% Habituation: total freezing (full session)
num_mice = size(mouse_files, 1);

% Load habituation freezing for each mouse
hab_freeze_data_all = cell(num_mice, 1);

for mouse = 1:num_mice
    current_freeze_data_path = freezing_folders{mouse};
    hab_file = dir(fullfile(current_freeze_data_path, '*habituation_freezing.mat'));
    
    if isempty(hab_file)
        warning('No habituation freezing file found for mouse %d', mouse);
        continue;
    end
    
    hab_to_load = fullfile(hab_file(1).folder, hab_file(1).name);
    S = load(hab_to_load);   % contains hab_freeze_data (logical vector)
    hab_freeze_data_all{mouse} = logical(S.hab_freeze_data);
end

% Helper: mean freezing for full session
get_mean = @(vec) mean(vec);   % because vec is logical, mean = fraction freezing

%% Extract means per treatment group

psi_hab_means = nan(numel(received_psilocybin), 1);
for i = 1:numel(received_psilocybin)
    idx = received_psilocybin(i);
    psi_hab_means(i) = get_mean(hab_freeze_data_all{idx});
end

veh_hab_means = nan(numel(received_vehicle), 1);
for i = 1:numel(received_vehicle)
    idx = received_vehicle(i);
    veh_hab_means(i) = get_mean(hab_freeze_data_all{idx});
end

% Drop NaNs in case of missing mice
psi_hab_means = psi_hab_means(~isnan(psi_hab_means));
veh_hab_means = veh_hab_means(~isnan(veh_hab_means));

%% Compute group means & SEMs (convert to %)

mean_psi_hab = mean(psi_hab_means) * 100;
sem_psi_hab  = std(psi_hab_means) / sqrt(numel(psi_hab_means)) * 100;

mean_veh_hab = mean(veh_hab_means) * 100;
sem_veh_hab  = std(veh_hab_means) / sqrt(numel(veh_hab_means)) * 100;

%% Bar plot: total habituation freezing

figure;
bar_vals = [mean_veh_hab, mean_psi_hab];
bar_sems = [sem_veh_hab, sem_psi_hab];

b = bar(1:2, bar_vals, 'FaceColor', 'flat'); hold on;

% Vehicle = blue, Psilocybin = red
b.CData(1,:) = [0 0.4470 0.7410];  % vehicle
b.CData(2,:) = [1 0 0];            % psilocybin

% Thick black outline on bars
b.EdgeColor = 'k';
b.LineWidth = 1.5;

% Error bars
errorbar(1:2, bar_vals, bar_sems, 'k.', 'LineWidth', 1.5);

set(gca, 'XTick', 1:2, 'XTickLabel', {'Vehicle','Psilocybin'});
ylabel('Total freezing during habituation (%)');
title(['Habituation freezing (full session mean): ' session]);
ylim([0 100]);
box off;

ax = gca;
ax.LineWidth = 1.5;
ax.XColor    = 'k';
ax.YColor    = 'k';
ax.TickDir   = 'out';
ax.Box       = 'off';


%% Habituation statistics: normality + unpaired t-test (built-in MATLAB only)

% Use the (already computed) psi_hab_means and veh_hab_means
psi_vals = psi_hab_means * 100;   % convert to %
veh_vals = veh_hab_means * 100;

% Normality tests (Lilliefors / KS)
[psi_h, psi_p] = lillietest(psi_vals);
[veh_h, veh_p] = lillietest(veh_vals);

fprintf('\nNormality tests (hab data) (Lilliefors):\n');
fprintf('  Psilocybin: p = %.4f (%s)\n', psi_p, ternary(psi_h==0, 'normal', 'non-normal'));
fprintf('  Vehicle:    p = %.4f (%s)\n', veh_p, ternary(veh_h==0, 'normal', 'non-normal'));

% Unpaired t-test (MATLAB defaults to Welch if variances differ)
[~, p_ttest, ~, stats] = ttest2(psi_vals, veh_vals);

fprintf('\nUnpaired t-test (hab data):\n');
fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_ttest);

% Compute Cohen's d
mean_diff = mean(psi_vals) - mean(veh_vals);

pooled_sd = sqrt( ...
    ((numel(psi_vals)-1)*var(psi_vals) + (numel(veh_vals)-1)*var(veh_vals)) ...
    / (numel(psi_vals) + numel(veh_vals) - 2));

d = mean_diff / pooled_sd;

fprintf('Cohen''s d = %.3f\n', d);

%%
function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
