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
N_trials  = 20; % total number of trials
fs = 30000;  % neuronal data sampling rate  
tpre = 0;        % Time before each event to include in the window (i.e. pre-stimulus time)
tpost = 20;       % Time after each event to include.               (i.e. post-stimulus time)

% neuronal data - converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

spk = double(spk)/fs;   % converts sample number (by dividing by sampling rate) into seconds
clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
N_neurons = numel(clu_val); % get number of detected clusters i.e neurons

% load list of Loom-response-types 
load([neuronal_data_path 'mouse3_extinction_LOOMresp_neurons'], 'resp_duringLOOM_nr', 'resp_postLOOM_nr')

% load extracted loom and flash (i.e. stimuli onset) events
load([triggers_path 'mouse3_extinction_extracted_events'],'evt_loom')

% Load freezing data
load('C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\freezing_looms_mouse3.mat')  
freeze_matrix = freezing_looms;

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

% define time variables to fit glmm for duing stim and post stim periods
period_length = 135;  
t_start = 0;                  

t_end = 9;            
num_bins = 120;    
binned_freezing_state = zeros(N_trials, num_bins);


% Set up empty matrix to store all neuron data
all_firing_rates = zeros(length(current_group), N_trials, num_bins);  % [neuron x trial x bin]

% Set bin size
bin_size = period_length / num_bins;   % frames
bin_size_s = bin_size / 15;            % seconds

% LOOP: For each neuron
for n = 1:length(current_group)
    idx = current_group(n);
    disp(['Processing neuron: ' num2str(idx-1)]);

    % Get spike times for this neuron
    tsp = spk(clu == clu_val(idx));

    % Get binned firing rates using raster_NM
    [~, ~, t_loom, fr, ~] = raster_NM(tsp, evt_loom, tpre, tpost, bin_size_s, false);

    % Keep only bins within specified analysis window
    idx_during = find(t_loom >= t_start & t_loom <= t_end);
    fr = fr(:, idx_during, :);  % [trial x bin]

    % Store in master matrix
    all_firing_rates(n, :, :) = fr;
end

% AFTER NEURON LOOP — reshape firing data
disp(['Binsize is ' num2str(bin_size_s) ' seconds'])  % display binsize



frames_start = round(t_start * 15) + 1;
frames_end   = round(t_end * 15);

cropped_freezing = freeze_matrix(:, frames_start:frames_end);  % [trial x frames]
frames_per_bin = round((frames_end - frames_start + 1) / num_bins);

for trial = 1:N_trials
    for bin = 1:num_bins
        start_idx = (bin-1)*frames_per_bin + 1;
        end_idx = min(bin*frames_per_bin, size(cropped_freezing, 2));
        binned_freezing_state(trial, bin) = ...
            mean(cropped_freezing(trial, start_idx:end_idx)) > 0.5;
    end
end


% Flatten for classification
y_state = reshape(binned_freezing_state, [], 1);

% [neurons x trials x bins] → [trials x bins x neurons]
X_temp = permute(all_firing_rates, [2, 3, 1]);   % size: [20, 60, 161]

% Flatten trials and bins into rows, neurons into columns
X = reshape(X_temp, [], size(all_firing_rates, 1));  % size: [1200, 161]

%%
% Mdl = fitcensemble(X, y_state, ...
%     'Method', 'Bag', ...
%     'NumLearningCycles', 100, ...
%     'Prior', 'uniform');  % Equal weight to both classes


%%

CVmdl = crossval(Mdl, 'KFold', 5);
loss = kfoldLoss(CVmdl);
fprintf('Cross-validated error: %.2f%%\n', loss*100);

%%
% Assume each 60 bins are in order from t = 0 to 9s
y_pred = kfoldPredict(CVmdl);


bin_labels = repmat(1:num_bins, N_trials, 1);  % Which time bin each row came from
bin_labels = bin_labels(:);

for b = 1:num_bins
    idx = bin_labels == b;
    acc_per_bin(b) = mean(y_pred(idx) == y_state(idx));
end


plot(1:num_bins, acc_per_bin, '-o');
xlabel('Time bin');
ylabel('Accuracy');
title('Decoding accuracy per time bin');

%%
freeze_per_bin = zeros(1, num_bins);
for b = 1:num_bins
    idx = bin_labels == b;
    freeze_per_bin(b) = mean(y_state(idx));
end

yyaxis left
plot(1:num_bins, acc_per_bin, '-o')
ylabel('Decoding Accuracy')

yyaxis right
plot(1:num_bins, freeze_per_bin, '-s')
ylabel('Proportion Freezing')

xlabel('Time bin')
title('Accuracy vs. Freezing Rate per Bin')
legend('Accuracy','% Freezing')

%%
confusionchart(y_state, y_pred);

%%
num_neurons = length(current_group);
accuracy_matrix = zeros(num_neurons, num_bins);  % [neurons x bins]

% [neurons x trials x bins] → [trials x bins x neurons]
X_temp = permute(all_firing_rates, [2, 3, 1]);  % [20, 60, neurons]

% Labels (freezing)
y = reshape(binned_freezing_state, [], 1);  % [1200 x 1]
bin_labels = repmat(1:num_bins, N_trials, 1);
bin_labels = bin_labels(:);  % [1200 x 1]

for n = 1:num_neurons
    fprintf('Neuron %d/%d\n', n, num_neurons);

    % Extract single neuron's data
    X_single = reshape(X_temp(:, :, n), [], 1);  % [1200 x 1]

    % Train classifier with class balancing
    try
        Mdl = fitcensemble(X_single, y, ...
            'Method', 'Bag', ...
            'NumLearningCycles', 100, ...
            'Prior', 'uniform');  % balance classes

        % Predict
        y_pred = predict(Mdl, X_single);

        % Accuracy per bin
        for b = 1:num_bins
            idx = bin_labels == b;
            accuracy_matrix(n, b) = mean(y_pred(idx) == y(idx));
        end
    catch
        warning('Could not train model for neuron %d (class imbalance or other issue)', n);
        accuracy_matrix(n, :) = NaN;
    end
end
%%
figure;
imagesc(accuracy_matrix);
xlabel('Time bin');
ylabel('Neuron');
title('Decoding accuracy per time bin (single-neuron models)');
colorbar;
caxis([0 1]);
colormap('winter');

%%
% Compare with null model:
% Calculate true accuracy per bin
acc_per_bin = zeros(1, num_bins);
for b = 1:num_bins
    idx = bin_labels == b;
    acc_per_bin(b) = mean(y_pred(idx) == y_state(idx));  % real decoding accuracy
end

% Compute shuffled (null) accuracy
num_shuffles = 100;
null_acc = zeros(num_shuffles, num_bins);

for s = 1:num_shuffles
    y_shuff = y_state(randperm(length(y_state)));  % shuffle labels
    for b = 1:num_bins
        idx = bin_labels == b;
        null_acc(s, b) = mean(y_shuff(idx) == y_state(idx));  % null accuracy
    end
end

% Mean and CI of shuffled distribution
null_mean = mean(null_acc, 1);
null_std  = std(null_acc, 0, 1);
ci_low  = null_mean - 1.96 * null_std / sqrt(num_shuffles);  % 95% CI lower
ci_high = null_mean + 1.96 * null_std / sqrt(num_shuffles);  % 95% CI upper

% Plot
figure;
hold on;
fill([1:num_bins, fliplr(1:num_bins)], ...
     [ci_low, fliplr(ci_high)], [0.8 0.8 0.8], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.5);  % 95% CI patch

plot(1:num_bins, null_mean, 'k--', 'LineWidth', 1.5);  % null mean
plot(1:num_bins, acc_per_bin, 'b-', 'LineWidth', 2);   % actual accuracy

xlabel('Time bin');
ylabel('Decoding accuracy');
legend({'Null 95% CI', 'Null mean', 'Actual decoding'}, 'Location', 'best');
title('Freezing decoding accuracy vs shuffled baseline');

% ------------------
% Report summary metrics
% ------------------
fprintf('\nAverage real accuracy: %.2f%%\n', mean(acc_per_bin)*100);
fprintf('Average null accuracy: %.2f%% ± %.2f%% (std)\n', ...
    mean(null_mean)*100, mean(null_std)*100);
%%