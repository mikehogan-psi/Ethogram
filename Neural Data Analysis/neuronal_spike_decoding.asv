%% Neuronal spike decoding analysis - FLASH vs LOOM
% this script .... (give summary)

%% Directory Setup

% Define folder path that contains neronal spiking data 
% (i.e. firing rate matrices and response groups obtained from neural_data_analysis_pipeline.m)

  neuronal_data_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\spiking_data\';

% Define folder path that contains kilosorted dats of this session
    kilosort_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\Extinction\kilosort4\';  

% Define folder path that contains loom and flash event trigger timnestamps
    triggers_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\concatinated_triggers\';

% Define folder path in which decoding results shall be saved
    savepath = 'Z:\Abi\neuronal_data\mouse_2\spike_decoding_analysis\';


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
   X = zscore(X, 0, 1);  % z-score across trials (column-wise)


% create labels to train classifier
    y = [ones(Ntrial,1); zeros(Ntrial,1)];  % 1 = loom, 0 = flash



%% train and test classifier

cv = cvpartition(y, 'KFold', 5);
accuracy = zeros(cv.NumTestSets, 1);

for i = 1:cv.NumTestSets
    trainIdx = training(cv, i);
    testIdx = test(cv, i);
    
    model = fitcsvm(X(trainIdx,:), y(trainIdx), 'KernelFunction', 'linear');
    
    y_pred = predict(model, X(testIdx,:));
    
    accuracy(i) = mean(y_pred == y(testIdx));
end

mean_acc = mean(accuracy);

fprintf('Mean decoding accuracy: %.2f%%\n', mean_acc*100);

% Visualiz Feature Weights

w = model.Beta;  % Weight vector from linear SVM
W = reshape(w, [size(X_trials,2), size(X_trials,3)]);  % [time_bin × neuron]

% Sort neurons by max weight or sum of weights to order by importance
[~, order] = sort(max(W,[],1), 'descend');
W_sorted = W(:, order);

figure;
imagesc(W_sorted);
xlabel('Neuron (sorted by max weight)');
ylabel('Time Bin');
title('Decoding weights - Radial Basis Function Kernel svm');
colorbar;

%% Permutation Test to Identify Significant decoding weights

n_permutations = 1000;
w_perm = zeros(n_permutations, size(X,2));

for p = 1:n_permutations
    y_perm = y(randperm(length(y)));  % shuffle labels
    model_perm = fitcsvm(X, y_perm, 'KernelFunction', 'linear');
    w_perm(p,:) = model_perm.Beta;
end

% Reshape actual weights and permuted weights
W_perm = reshape(w_perm, n_permutations, size(W,1), size(W,2));  % [n_perm × time_bins × neurons]
max_abs_perm = squeeze(max(abs(W_perm),[],2));  % max over time bins, size: [n_perm x neurons]
max_abs_true = max(abs(W), [], 1);              % [1 x neurons]

% Compute p-values per neuron
p_vals = mean(max_abs_perm >= max_abs_true, 1);

% Identify neurons significantly contributing (e.g., p < 0.05)
alpha = 0.05;
signif_neurons = p_vals < alpha;
signif_neurons_idx = find(signif_neurons);
fprintf('%d neurons significantly contribute to decoding (p < %.2f)\n', sum(signif_neurons), alpha);

% Optional: Visualize only significant neurons
W_signif = W(:, signif_neurons);
[~, order] = sort(max(abs(W_signif),[],1), 'descend');
W_signif_sorted = W_signif(:, order);

figure;
imagesc(W_signif_sorted);
xlabel('Neuron (sorted, significant only)');
ylabel('Time Bin');
title('Decoding Weights (Significant Neurons Only)');
colorbar;

% save p-values and significant neurons idexes
   save([savepath 'population_decoding'], 'W_signif', 'signif_neurons_idx' );


%% Neuron-by-Neuron Decoding with Permutation to test decoding accuracy

n_perm = 1000;      % Number of permutations for significance testing
p_thresh = 0.05;    % Significance threshold
n_neurons = size(X_trials, 3);
n_bins = size(X_trials, 2);

real_acc = zeros(n_neurons, 1);
perm_acc = zeros(n_neurons, n_perm);
p_vals = zeros(n_neurons, 1);

fprintf('Running decoding for %d neurons...\n', n_neurons);

for n = 1:n_neurons
    % Extract time-binned activity for neuron n: trials x timebins
    Xn = squeeze(X_trials(:, :, n));  
    
    % Flatten into trials x features (i.e. one feature vector per trial)
    Xn_flat = zscore(Xn, 0, 1);           % Z-score across trials
    Xn_flat = Xn_flat(:,:) ;              % trials × timebins

    % Cross-validated real decoding accuracy
    cv = cvpartition(y, 'KFold', 5);
    acc_cv = zeros(cv.NumTestSets, 1);
    for i = 1:cv.NumTestSets
        trainIdx = training(cv, i);
        testIdx = test(cv, i);
        
        mdl = fitcsvm(Xn_flat(trainIdx,:), y(trainIdx), 'KernelFunction', 'linear');
        y_pred = predict(mdl, Xn_flat(testIdx,:));
        acc_cv(i) = mean(y_pred == y(testIdx));
    end
    real_acc(n) = mean(acc_cv);
    
    % Permutation testing
    for p = 1:n_perm
        y_perm = y(randperm(length(y)));
        acc_perm_cv = zeros(cv.NumTestSets, 1);
        
        for i = 1:cv.NumTestSets
            trainIdx = training(cv, i);
            testIdx = test(cv, i);
            
            mdl = fitcsvm(Xn_flat(trainIdx,:), y_perm(trainIdx), 'KernelFunction', 'linear');
            y_pred = predict(mdl, Xn_flat(testIdx,:));
            acc_perm_cv(i) = mean(y_pred == y_perm(testIdx));
        end
        perm_acc(n, p) = mean(acc_perm_cv);
    end
    
    % Compute p-value
    p_vals(n) = mean(real_acc(n) <= perm_acc(n, :));

    fprintf('\n neuron %d permutation test: p = %.4g', n, p_vals(n));

end

% Identify significant neurons
signif_neurons = find(p_vals < p_thresh);

fprintf('\n%d out of %d neurons decode above chance (p < %.2f)\n', ...
    numel(signif_neurons), n_neurons, p_thresh);

% Optional: plot decoding accuracy vs. null distribution
figure;
histogram(real_acc, 'BinWidth', 0.01, 'FaceColor', 'b');
hold on;
histogram(perm_acc(:), 'BinWidth', 0.01, 'FaceColor', 'r', 'FaceAlpha', 0.4);
xlabel('Decoding Accuracy');
ylabel('Count');
legend('Real Neuron Accuracy', 'Permuted Labels');
title('Neuron-by-Neuron Decoding Accuracy Distribution');

% save p-values and significant neurons idexes
   save([savepath 'single_neuron_decoding'], 'real_acc', 'perm_acc', 'p_vals', 'signif_neurons');


%% Time-Resolved Decoding for Each Significant Neuron

% Assume:
% X_trials: [n_trials x n_timebins x n_neurons]
% y: trial labels (1 = loom, 0 = flash)
% signif_neurons: indices of significant neurons from previous analysis
% bin_labels: vector of time points or bin indices for plotting

n_timebins = size(X_trials, 2);
n_trials = size(X_trials, 1);
n_sig = length(signif_neurons);

bin_acc = zeros(n_sig, n_timebins);  % decoding accuracy per bin per neuron

for i = 1:n_sig
    n_idx = signif_neurons(i);
    for t = 1:n_timebins
        % Extract feature: single time bin for this neuron across trials
        X_t = squeeze(X_trials(:, t, n_idx));  % [n_trials x 1]
        X_t = zscore(X_t);  % normalize across trials
        
        % Train/test decoding with cross-validation
        cv = cvpartition(y, 'KFold', 5);
        acc = zeros(cv.NumTestSets, 1);
        for k = 1:cv.NumTestSets
            trainIdx = training(cv, k);
            testIdx = test(cv, k);
            
            model = fitcsvm(X_t(trainIdx), y(trainIdx), 'KernelFunction', 'linear');
            y_pred = predict(model, X_t(testIdx));
            acc(k) = mean(y_pred == y(testIdx));
        end
        bin_acc(i, t) = mean(acc);
    end
end

% Plot time-resolved decoding accuracy for all significant neurons (heatmap)
figure;
imagesc(bin_acc);
xlabel('Time Bin');
ylabel('Neuron (significant only)');
title('Time-Resolved Decoding Accuracy per Neuron');
colorbar;

% save time-binned accuracies of significant neurons 
   save([savepath 'single_neuron_decoding'], 'bin_acc', '-append');




%% saving individual graphs accuracy over time-bins for one neuron

for n = 1:n_sig 

sig_neuron_idx = n;
neuron_acc = bin_acc(sig_neuron_idx,:);
smoothed_acc = movmean(neuron_acc, 10);  % smoothing

figure;
hold on
plot(neuron_acc, 'LineWidth', 0.5);
plot(smoothed_acc, 'LineWidth', 2);
yline(0.5, '--k');  % Chance level
xlabel('Time Bin');
ylabel('Decoding Accuracy');
title(sprintf('Time-Resolved Decoding (Neuron %d)', signif_neurons(sig_neuron_idx)));
legend('raw', 'smoothed')
grid on
hold off

fig = gcf; 
saveas(fig, fullfile('Z:\Abi\neuronal_data\mouse_2\spike_decoding_analysis\figures\single_neuron_time_res_decoding', ...
            sprintf('Neuron_%d.png', signif_neurons(sig_neuron_idx))));
 

close

end








