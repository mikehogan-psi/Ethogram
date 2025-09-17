%% Directory setup

% Define folder path that contains neuronal spiking data 
% (i.e. firing rate matrices and response groups obtained from neural_data_analysis_pipeline.m)

  neuronal_data_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\processed_data\spiking_data\';

% Define folder path that contains kilosorted dats of this session
    kilosort_dir = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\kilosort4\'; % 


    triggers_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\processed_data\concatenated_triggers\';

% % Define filepath to save processed data to and name of saved data
%     savefile_name_glmm = 'extinction_neuronal_cofficients_poisson_';
%     savefile_name_data_tables = 'extinction_data_tables_';
%     save_folder = 'Z:\Abi\neuronal_data\mouse_2\Mike processed data\';

%% Load data

% get important variables
fs = 30000;  % neuronal data sampling rate  
tpre = 3;        % Time before each event to include in the window (i.e. pre-stimulus time)
tpost = 3;       % Time after each event to include.               (i.e. post-stimulus time)

% neuronal data - converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

spk = double(spk)/fs;   % converts sample number (by dividing by sampling rate) into seconds
clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
N_neurons = numel(clu_val); % get number of detected clusters i.e neurons

% load list of Loom-response-types 
load([neuronal_data_path 'mouse3_extinction_LOOMresp_neurons'], 'resp_duringLOOM_nr', 'resp_postLOOM_nr')

% load extracted loom and flash (i.e. stimuli onset) events
load([neuronal_data_path 'mouse_3_extinction_freezing_aligned_loom_neural_data'],'valid_onset_times', 't_loom')

% Create variable with all responsive neurons
resp_LOOM_nr = [resp_duringLOOM_nr; resp_postLOOM_nr];
resp_LOOM_nr = unique(resp_LOOM_nr);


% % Specify which type of responsive neurons you want to look at and save
% % name
current_group = resp_LOOM_nr;
% append_to_save = 'resp_LOOM_nr.mat';
% 
% % Generates filepaths for saving   
% savefile_name_glmm = [savefile_name_glmm,  append_to_save];
% save_path_glmm = fullfile(save_folder, savefile_name_glmm);
% 
% savefile_name_data_tables = [savefile_name_data_tables, append_to_save];
% save_path_data_tables = fullfile(save_folder, savefile_name_data_tables);
%% Preprocess data for fitting ensemble of learners      
N_trials  = numel(valid_onset_times); % total number of trials
num_bins = numel(t_loom);

% Set up empty matrix to store all neuron data
all_firing_rates = zeros(length(current_group), N_trials, num_bins);  % [neuron x trial x bin]          % seconds

% Get bin size
bin_size_s = (tpre + tpost)/num_bins;

% LOOP: For each neuron
for n = 1:length(current_group)
    idx = current_group(n);
    disp(['Processing neuron: ' num2str(idx-1)]);

    % Get spike times for this neuron
    tsp = spk(clu == clu_val(idx));

    % Get binned firing rates using raster_NM
    [~, ~, t_loom, fr, ~] = raster_NM(tsp, valid_onset_times, tpre, tpost, bin_size_s, false);

    % Store in master matrix
    all_firing_rates(n, :, :) = fr;
end




% AFTER NEURON LOOP — reshape firing data
disp(['Number of bins is ' num2str(num_bins)])  % display number of bins
disp(['Binsize is ' num2str(bin_size_s) ' seconds'])  % display binsize


freezing_vector = [zeros(1, round(num_bins/2)), ones(1, floor(num_bins/2))];
middle_idx = round(numel(freezing_vector)/2);
freezing_vector(middle_idx-2:middle_idx+2) = 2;

freezing_matrix = repmat(freezing_vector, N_trials, 1);


% Flatten for classification
y_state = reshape(freezing_matrix, [], 1);

% [neurons x trials x bins] → [trials x bins x neurons]
X_temp = permute(all_firing_rates, [2, 3, 1]);  

% Flatten trials and bins into rows, neurons into columns
X = reshape(X_temp, [], size(all_firing_rates, 1)); 

%%
% Fit an ensemble classifier (e.g., bagged trees)
Mdl = fitcensemble(X, y_state, 'Method', 'Bag');
%%
% 10-fold cross-validation
cvMdl = crossval(Mdl, 'KFold', 10);

% Predict and calculate classification error
classLoss = kfoldLoss(cvMdl);  % classification error
accuracy = 1 - classLoss;
disp(['Cross-validated accuracy: ' num2str(accuracy*100, '%.2f') '%'])

%%
% Get predictions
y_pred = kfoldPredict(cvMdl);

% Confusion matrix
figure;
confusionchart(y_state, y_pred);
title('Confusion Matrix');

%%
% Estimate predictor importance
imp = predictorImportance(Mdl);

% Plot importance
figure;
bar(imp);
xlabel('Neuron');
ylabel('Importance');
title('Neuron Importance in Freezing Prediction');

%%
[~, scores] = kfoldPredict(cvMdl);

% Only use the second column (probability for class 1)
[fpRate, tpRate, ~, AUC] = perfcurve(y_state, scores(:,2), 1);

figure;
plot(fpRate, tpRate, 'b-', 'LineWidth', 2);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC Curve (AUC = ' num2str(AUC, '%.2f') ')']);
grid on;

%%
% Refit on full dataset (not crossval)
Mdl_full = fitcensemble(X, y_state, 'Method', 'Bag');

% Get feature importance
importance = oobPermutedPredictorImportance(Mdl_full);

% Plot top contributing neurons
[~, sortIdx] = sort(importance, 'descend');
figure;
bar(importance(sortIdx));
xlabel('Neuron Index (sorted)');
ylabel('Importance Score');
title('Neuron Contributions to Freezing Prediction');

%%

% Define epochs relative to freeze onset (in seconds)
pre_epoch = [-3, -1];
transition_epoch = [-1, 1];
sustain_epoch = [1, 3];

% Convert epochs to bin indices (assuming t_loom in seconds)
pre_bins = find(t_loom >= pre_epoch(1) & t_loom <= pre_epoch(2));
trans_bins = find(t_loom >= transition_epoch(1) & t_loom <= transition_epoch(2));
sustain_bins = find(t_loom >= sustain_epoch(1) & t_loom <= sustain_epoch(2));

% Average firing rates in each epoch per neuron and trial
X_pre = squeeze(mean(all_firing_rates(:, :, pre_bins), 3))';      % [trial x neuron]
X_trans = squeeze(mean(all_firing_rates(:, :, trans_bins), 3))'; % [trial x neuron]
X_sustain = squeeze(mean(all_firing_rates(:, :, sustain_bins), 3))'; % [trial x neuron]

% Labels for decoding
y_pre_vs_trans = [zeros(N_trials,1); ones(N_trials,1)]; % 0=pre-freeze, 1=transition
X_pre_trans = [X_pre; X_trans];  % stack data for decoding

% Train classifier for transition detection
Mdl_trans = fitcensemble(X_pre_trans, y_pre_vs_trans, 'Method', 'Bag');
cvMdl_trans = crossval(Mdl_trans, 'KFold', 10);
acc_trans = 1 - kfoldLoss(cvMdl_trans);
fprintf('Transition decoding accuracy: %.2f%%\n', acc_trans*100);

% Similar for sustained freezing: decode pre vs sustained
y_pre_vs_sustain = [zeros(N_trials,1); ones(N_trials,1)];
X_pre_sustain = [X_pre; X_sustain];
Mdl_sustain = fitcensemble(X_pre_sustain, y_pre_vs_sustain, 'Method', 'Bag');
cvMdl_sustain = crossval(Mdl_sustain, 'KFold', 10);
acc_sustain = 1 - kfoldLoss(cvMdl_sustain);
fprintf('Sustained freezing decoding accuracy: %.2f%%\n', acc_sustain*100);

% You can also inspect feature importance from these models to find key neurons
imp_trans = predictorImportance(Mdl_trans);
imp_sustain = predictorImportance(Mdl_sustain);

%%

pre_idx = find(t_loom >= -3 & t_loom < -1);
trans_idx = find(t_loom >= -1 & t_loom <= 1);
post_idx = find(t_loom > 1 & t_loom <= 3);

% Initialize matrices: neurons x trials
mean_pre = squeeze(mean(all_firing_rates(:, :, pre_idx), 3));    % [neurons x trials]
mean_trans = squeeze(mean(all_firing_rates(:, :, trans_idx), 3));
mean_post = squeeze(mean(all_firing_rates(:, :, post_idx), 3));

p_values = zeros(length(current_group), 1);

for n = 1:length(current_group)
    % Data matrix for this neuron: trials x epochs (3)
    data = [mean_pre(n, :)', mean_trans(n, :)', mean_post(n, :)'];
    
    % Use Friedman test (nonparametric repeated measures)
    p = friedman(data, 1, 'off');
    
    p_values(n) = p;
end

% Find neurons significantly modulated by epoch
signif_neurons = find(p_values < 0.3);

fprintf('%d neurons significantly modulated by freezing epoch\n', length(signif_neurons));

for i = 1:length(signif_neurons)


    neuron_to_plot = signif_neurons(i); % first significant neuron
    
    mean_rates = [mean(mean_pre(neuron_to_plot, :)), ...
                  mean(mean_trans(neuron_to_plot, :)), ...
                  mean(mean_post(neuron_to_plot, :))];
    
    figure;
    bar(mean_rates);
    set(gca, 'XTickLabel', {'Pre-freeze', 'Freeze transition', 'Sustained freeze'});
    ylabel('Mean firing rate (Hz)');
    title(['Neuron ' num2str(neuron_to_plot) ' firing rate by freezing epoch']);

end

%%

figure;
for i = 1:length(signif_neurons)
    n = signif_neurons(i);
    
    subplot(ceil(length(signif_neurons)/2), 2, i)
    plot(t_loom, mean(squeeze(all_firing_rates(n, :, 1:end-1)), 1), 'k', 'LineWidth', 1.5);
    title(['Neuron ' num2str(n)]);
    xlabel('Time from freeze onset (s)');
    ylabel('Mean firing rate (Hz)');
    xline(0, '--r');  % freeze onset
    ylim([0 max(max(squeeze(all_firing_rates(n, :, :)))) * 1.2]);
end

%%
