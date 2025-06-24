%% Neuronal spike decoding analysis

%% Directory Setup

% Define folder path that contains neronal spiking data 
% (i.e. firing rate matrices and response groups obtained from neural_data_analysis_pipeline.m)

  neuronal_data_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\spiking_data\';

% Define folder path that contains kilosorted dats of this session
    kilosort_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\Extinction\kilosort4\';  

% Define folder path that contains loom and flash event trigger timnestamps
    triggers_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\concatinated_triggers\';

%% Loading neuronal spike and trial data

% converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

    spk = double(spk)/30000;% converts sample number (by dividing through sampling rate) into seconds
    clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
    Ncell = numel(clu_val); % get number of detected clusters 

% load extracted loom and flash (i.e. stimuli onset) events
load([triggers_path 'mouse2_extinction_extracted_events'],'evt_loom', 'evt_flash')

    Ntrial = length(evt_loom);

%% Extracting number of spikes per neuron in loom and flash trials

% time window options
    tpre = 10;         % Time before each event to include in the window (i.e. pre-stimulus time)
    tpost = 23.3;       % Time after each event to include.               (i.e. post-stimulus time)
    bin_size = 0.1;     % Bin size (time resolution of histogram)


% initialize variables to hold spikecounts [trials x time-bins x cell] 
    binned_skp_counts_loom  = [];  % Nr of spikes in each time-bin of loom  trials 
    binned_skp_counts_flash = [];  % Nr of spikes in each time-bin of loom  trials 
   

for n = 1:Ncell
    tsp = spk(clu == clu_val(n));

    [~, ~, ~, ~, ~, ~, ~, binned_skp_counts_loom(:,:,n)]  = raster_NM(tsp, evt_loom, tpre, tpost, bin_size, false);
    [~, ~, ~, ~, ~, ~, ~, binned_skp_counts_flash(:,:,n)] = raster_NM(tsp, evt_flash, tpre, tpost, bin_size, false);
    
end

% Trim edges (Removes the first and last bins of the histogram to  avoid edge effects)
    binned_skp_counts_loom = binned_skp_counts_loom(:,2:end-1,:);
    binned_skp_counts_flash = binned_skp_counts_flash(:,2:end-1,:);
    
%% Create Classifier input 

% concatinate loom and flash matrices to combines binned population spike matrices 
    X_trials = cat(1, binned_skp_counts_loom, binned_skp_counts_flash);  % [40 × n_bins × n_cells]

% Flatten time X_trials into a 2D feature matrix (trials × features):    
   [n_trials, n_bins, n_cells] = size(X_trials);
   X = reshape(X_trials, [n_trials, n_bins * n_cells]);  % [n_trials × (n_bins × n_cells)]

% create labels to train classifier
    y = [ones(Ntrial,1); zeros(Ntrial,1)];  % 1 = loom, 0 = flash

%% train and test classifier (whith standardizing between trials)

cv = cvpartition(y, 'KFold', 5);
accuracy = zeros(cv.NumTestSets, 1);

for i = 1:cv.NumTestSets
    trainIdx = training(cv, i);
    testIdx = test(cv, i);

    % Compute mean and std only from training data
    mu = mean(X(trainIdx,:), 1);
    sigma = std(X(trainIdx,:), [], 1);
    
    % Avoid division by zero in case of zero variance
    sigma(sigma == 0) = 1;

    % Standardize
    X_train = (X(trainIdx,:) - mu) ./ sigma;
    X_test  = (X(testIdx,:)  - mu) ./ sigma;

    % Train and test SVM
    model = fitcsvm(X_train, y(trainIdx), 'KernelFunction', 'linear');
    y_pred = predict(model, X_test);
    
    accuracy(i) = mean(y_pred == y(testIdx));
end


mean_acc = mean(accuracy);

fprintf('Mean decoding accuracy: %.2f%%\n', mean_acc*100);

%% train and test classifier (whithout standardizing between trials)

cv = cvpartition(y, 'KFold', 5);
accuracy = zeros(cv.NumTestSets, 1);

for i = 1:cv.NumTestSets
    trainIdx = training(cv, i);
    testIdx = test(cv, i);
    
    model = fitcsvm(X(trainIdx,:), y(trainIdx), 'KernelFunction', 'linear');
    y_pred = predict(model, X(testIdx,:));
    
    accuracy(i) = mean(y_pred == y(testIdx));
end

mean_acc = mean(acc);

fprintf('Mean decoding accuracy: %.2f%%\n', mean_acc*100);


%% Interpret Feature Weights

w = model.Beta;  % Weight vector from linear SVM
W = reshape(w, [size(X_trials,2), size(X_trials,3)]);  % [time_bin × neuron]
figure 
imagesc(W);
xlabel('Neuron'); ylabel('Time Bin');
title('Decoding weights');
colorbar;

%%
% Define threshold for significance (e.g., > 0 to push towards loom)
thresh = 0.00025;

% Find neurons with any time bin weight > threshold
signif_neurons = any(W > thresh, 1); % 1 x neurons logical vector

% Extract weights for those neurons only
W_signif = W(:, signif_neurons);

% Sort neurons by max weight or sum of weights to order by importance
[~, order] = sort(max(W_signif,[],1), 'descend');
W_signif_sorted = W_signif(:, order);

figure;
imagesc(W_signif_sorted);
xlabel('Neuron (sorted by max weight)');
ylabel('Time Bin');
title('Decoding weights: neurons pushing towards loom');
colorbar;




















