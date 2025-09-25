% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!!Provide session to be analysed!!!
session = 'Extinction';

trial_type = 'looms'; % 'looms' or 'flashes'

% !!! Provide treatment mouse numbers for each treatment group !!!
received_psilocybin = [3; 5; 7];
received_vehicle = [1; 2; 4; 6];
%%
% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

% Initialise 
freezing_folders = cell(length(mouse_files), 1);

% Get DLC folders for specified session
for mouse = 1:length(mouse_files)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    freezing_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Freezing');
end

freeze_data = cell(length(mouse_files), 1);

for mouse = 1:length(freezing_folders)
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
colorbar; xlabel('Frame'); ylabel('Trial'); title('Psilocybin: mean freezing (trial x frame)');

subplot(1,2,2);
imagesc(veh_mean_over_mice); axis xy;
colorbar; xlabel('Frame'); ylabel('Trial'); title('Vehicle: mean freezing (trial x frame)');
cb2 = colorbar; 
ylabel(cb2, 'Freezing Probability');
sgtitle([session, ' ', trial_type])
% Optional: smoothing across frames (uncomment to use)
% smooth_window = 5;
% psi_frame_mean = movmean(psi_frame_mean, smooth_window);
% veh_frame_mean = movmean(veh_frame_mean, smooth_window);
