%% Neuronal spike decoding analysis - BEHAVIOUR vs NO-BEHAVIOUR
% this script .... (give summary)

%% Directory Setup

% Define folder path that contains neronal spiking data 
% (i.e. firing rate matrices and response groups obtained from neural_data_analysis_pipeline.m)

  neuronal_data_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\spiking_data\';

% Define folder path that contains kilosorted dats of this session
    kilosort_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\Extinction\kilosort4\';  

% Define folder path that contains loom and flash event trigger timnestamps
    triggers_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\concatinated_triggers\';   

% Define folder path that contains behavioural labels
    behaviour_path = 'Z:\Abi\behavioral_analysis\Implanted_mice\mouse_2_behaviours\';

% Define folder path in which decoding results shall be saved
    savepath = 'Z:\Abi\neuronal_data\mouse_2\spike_decoding_analysis\';

%% define variables

 % define target behaviour that shall be analysed
      behaviour = 'freezing';

 % define session 
      sesh = 'extinction';
  %   sesh = 'renewal';


%% Loading neuronal spike, trial and behavioural data

% converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

    spk = double(spk)/30000;% converts sample number (by dividing through sampling rate) into seconds
    clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
    Ncell = numel(clu_val); % get number of detected clusters 

% load extracted triggers for the session part you want to analyse
    % load part 1
    load([triggers_path 'mouse2_' sesh '_p1.mat'],'evt');
        evt_p1 = evt(2:end); % remove first trigger (artefact)

    % load part 2
    load([triggers_path 'mouse2_' sesh '_p2.mat'],'evt');
        evt_p2 = evt(2:end); % remove first trigger (artefact)

    % concatinate part 1 and part 2
    evt = [evt_p1;evt_p2];   

    % convert into seconds
    evt = double(evt)/30000; % converts sample number (by dividing through sampling rate) into seconds
                             % [20080x1] - contains timestamp (sec) of every frame

% load behavioural labels 
    % load part 1
    load([ behaviour_path behaviour '_labels_mouse2_' sesh '_p1.mat'], 'predicted_labels');
            labels_p1 = predicted_labels;
    % load part 2        
    load([ behaviour_path behaviour '_labels_mouse2_' sesh '_p2.mat'], 'predicted_labels');
            labels_p2 = predicted_labels;
    % concatinate part 1 and part 2
        labels = [labels_p1;labels_p2]; % [20080x1] - contains behaviour label (0 = no behaviour, 1 = behaviour) for every frame

% clear unecessary variables
clearvars evt_p1 evt_p2  labels_p1 labels_p2 predicted_labels 

%% Extracting number of spikes per neuron in each frame

% Define window duration per frame 
frame_duration = median(diff(evt));  % ~0.0667 sec for 15 Hz sampling rate

% Pre-allocate
X = zeros(length(evt), Ncell);  % Rows = frames, Cols = neurons
y = labels;  % Already aligned to evt

% For each frame, count spikes per neuron
for i = 1:length(evt)
    t_start = evt(i);
    t_end = t_start + frame_duration;
    
    % Get spikes in this window
    in_window = spk >= t_start & spk < t_end;
    spikes_now = clu(in_window);
    
    % Count spikes per neuron
    counts = histcounts(spikes_now, [clu_val; max(clu_val)+1]);
    
    % load into firing rate matrix
    X(i, :) = counts; % [20080 x Nneurons]

end

%% identify which neurons can significantly decode behaviour
%  Neuron-by-Neuron Decoding with Permutation to test decoding accuracy

% Parameters
n_perm = 1000;             % Number of permutations for significance testing
alpha = 0.05;              % Significance threshold
n_neurons = size(X, 2);    % Total number of neurons
real_acc = zeros(n_neurons, 1);
perm_acc = zeros(n_neurons, n_perm);

fprintf('Running single-neuron decoding and permutation test...\n');

for n = 4:n_neurons

    fprintf('Running single-neuron decoding and permutation test for neuron %d\n', n);

    % Extract spike counts for neuron n
    Xn = X(:, n);  % [frames x 1]
    
    % Normalize (optional, but may help SVM performance)
    Xn = zscore(Xn);

    % 5-fold cross-validation accuracy on real data
    cv = cvpartition(y, 'KFold', 5); % s
    acc_cv = zeros(cv.NumTestSets, 1);

    for fold = 1:cv.NumTestSets
        trainIdx = training(cv, fold);
        testIdx = test(cv, fold);

        mdl = fitcsvm(Xn(trainIdx), y(trainIdx), 'KernelFunction', 'linear');
        y_pred = predict(mdl, Xn(testIdx));
        acc_cv(fold) = mean(y_pred == y(testIdx));
    end

    real_acc(n) = mean(acc_cv);

    % Permutation test
    for p = 1:n_perm

        fprintf('Running permutation %d\n', p);


        y_perm = y(randperm(length(y)));
        acc_p = zeros(cv.NumTestSets, 1);

        for fold = 1:cv.NumTestSets
            trainIdx = training(cv, fold);
            testIdx = test(cv, fold);

            mdl = fitcsvm(Xn(trainIdx), y_perm(trainIdx), 'KernelFunction', 'linear');
            y_pred = predict(mdl, Xn(testIdx));
            acc_p(fold) = mean(y_pred == y_perm(testIdx));
        end

        perm_acc(n, p) = mean(acc_p);
    end
end

% Compute p-values
p_vals = mean(bsxfun(@ge, perm_acc, real_acc), 2);  % For each neuron

% Find significant neurons
signif_neurons = find(p_vals < alpha);

fprintf('\n%d out of %d neurons decode behaviour significantly (p < %.2f)\n', ...
        numel(signif_neurons), n_neurons, alpha);

% save p-values and significant neurons idexes
   save([savepath behaviour '_signle_neuron_freezing_decoding'], 'real_acc', 'perm_acc', 'p_vals', 'signif_neurons');




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train and test classifier: behavior vs. no-behavior

% Balance the data
inds_behaviour = find(y == 1);
inds_nobehaviour = find(y == 0);
n = min(numel(inds_behaviour), numel(inds_nobehaviour));

% Subsample to balance
inds_use = [inds_behaviour(randperm(numel(inds_behaviour), n)); ...
            inds_nobehaviour(randperm(numel(inds_nobehaviour), n))];

% Shuffle
inds_use = inds_use(randperm(numel(inds_use)));

% Extract balanced data
X_bal = X(inds_use, :);
y_bal = y(inds_use);

% Create 5-fold cross-validation on balanced labels
cv = cvpartition(y_bal, 'KFold', 5);
accuracy = zeros(cv.NumTestSets, 1);

for i = 1:cv.NumTestSets
    trainIdx = training(cv, i);
    testIdx = test(cv, i);
    
    % Train linear SVM
    model = fitcsvm(X_bal(trainIdx,:), y_bal(trainIdx), 'KernelFunction', 'linear');
    
    % Predict on test set
    y_pred = predict(model, X_bal(testIdx,:));
    
    % Compute accuracy
    accuracy(i) = mean(y_pred == y_bal(testIdx));
end

% Mean accuracy across folds
mean_acc = mean(accuracy);
fprintf('Mean decoding accuracy: %.2f%%\n', mean_acc*100);

%% Visualize SVM Feature Weights

% Get weights (only valid for linear SVM)
w = model.Beta;  % [neurons x 1]

% % Sort neurons by absolute weight
% [~, order] = sort(abs(w), 'descend');
% w_sorted = w(order);

% Plot
figure;
bar(w);
xlabel('Kilosort cluster');
ylabel('SVM Weight');
title(['Neuron Importance for Decoding ' behaviour]);



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
   save([savepath 'population_decoding_' behaviour], 'W_signif', 'signif_neurons_idx' );


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

% extract decoding accuracy for current neuron
sig_neuron_idx = n;
neuron_acc = bin_acc(sig_neuron_idx,:);
smoothed_acc = movmean(neuron_acc, 10);  % smoothing

% convert time-bins back into seconds
time = (1:length(neuron_acc)) * bin_size - tpre;  % Align to stimulus time

% plot figure
figure;
hold on
    plot(time, neuron_acc, 'LineWidth', 0.5);
    plot(time, smoothed_acc, 'LineWidth', 2);
    xline(0, '-k');     % stimulus start
    yline(0.5, '--k');  % Chance level
    xlabel('Time (s)');
    ylabel('Decoding Accuracy');
    title(sprintf('Time-Resolved Decoding (Neuron %d)', signif_neurons(sig_neuron_idx)));
    legend('raw', 'smoothed')
    grid on
hold off

%save figure
% fig = gcf; 
% saveas(fig, fullfile('Z:\Abi\neuronal_data\mouse_2\spike_decoding_analysis\figures\single_neuron_time_res_decoding', ...
%             sprintf('Neuron_%d.png', signif_neurons(sig_neuron_idx))));
% close

end


%% Computing statistical significance of time-bins

% Settings
neuron = 2; 
n_perms = 1000;
chance_level = 0.5;
alpha = 0.05;  % Significance threshold

% Extract time-resolved features for the neuron
X_neuron = squeeze(X_trials(:,:,neuron));  % [trials × time-bins]
[n_trials, n_bins] = size(X_neuron);

% Real accuracies
real_acc = zeros(1, n_bins);
null_acc = zeros(n_perms, n_bins);

for t = 1:n_bins
    real_acc(t) = mean(crossval(@(Xtrain,Ytrain,Xtest,Ytest) ...
        mean(predict(fitcsvm(Xtrain,Ytrain,'KernelFunction','linear'), Xtest) == Ytest), ...
        X_neuron(:,t), y, 'KFold', 5));
    
    for p = 1:n_perms
        y_perm = y(randperm(length(y)));
        null_acc(p,t) = mean(crossval(@(Xtrain,Ytrain,Xtest,Ytest) ...
            mean(predict(fitcsvm(Xtrain,Ytrain,'KernelFunction','linear'), Xtest) == Ytest), ...
            X_neuron(:,t), y_perm, 'KFold', 5));
    end
end

% Compute p-values
p_vals = mean(null_acc >= real_acc, 1);  % one-tailed test
significant_bins = p_vals < alpha;

% Plot
time = (1:n_bins)*bin_size - tpre;  % adjust to actual time axis
figure;
plot(time, real_acc, 'b', 'LineWidth', 2); hold on;
plot(time(significant_bins), real_acc(significant_bins), 'ro', 'MarkerFaceColor','r');
yline(chance_level, '--k', 'Chance');
xlabel('Time (s)');
ylabel('Decoding Accuracy');
title(sprintf('Time-Resolved Decoding w/ Significance (Neuron %d)', sig_neuron_idx));
legend('Accuracy', 'Significant Bin');
grid on;

%% Computing statistical significance of time-bins (faster)

% Settings
neuron_idx = 2; 
n_perms = 1000;
chance_level = 0.5;
alpha = 0.05;  % Significance threshold


X_neuron = squeeze(X_trials(:,:,neuron_idx));
[n_trials, n_bins] = size(X_neuron);
real_acc = zeros(1, n_bins);
null_acc = zeros(n_perms, n_bins);

cv = cvpartition(y, 'KFold', 5);
folds = arrayfun(@(i) struct('trainIdx', training(cv,i), 'testIdx', test(cv,i)), 1:cv.NumTestSets);

y_perms = zeros(length(y), n_perms);
for p = 1:n_perms
    y_perms(:,p) = y(randperm(length(y)));
end

parfor t = 1:n_bins
    X_t = X_neuron(:, t);

    % Real accuracy
    acc_real = zeros(1, cv.NumTestSets);
    for i = 1:cv.NumTestSets
        mdl = fitcsvm(X_t(folds(i).trainIdx), y(folds(i).trainIdx), 'KernelFunction', 'linear');
        yhat = predict(mdl, X_t(folds(i).testIdx));
        acc_real(i) = mean(yhat == y(folds(i).testIdx));
    end
    real_acc(t) = mean(acc_real);

    % Null distribution
    for p = 1:n_perms
        acc_perm = zeros(1, cv.NumTestSets);
        for i = 1:cv.NumTestSets
            y_perm = y_perms(:,p);
            mdl = fitcsvm(X_t(folds(i).trainIdx), y_perm(folds(i).trainIdx), 'KernelFunction', 'linear');
            yhat = predict(mdl, X_t(folds(i).testIdx));
            acc_perm(i) = mean(yhat == y_perm(folds(i).testIdx));
        end
        null_acc(p,t) = mean(acc_perm);
    end
end


% Compute p-values
p_vals = mean(null_acc >= real_acc, 1);  % one-tailed test
significant_bins = p_vals < alpha;

% Plot
time = (1:n_bins)*bin_size - tpre;  % adjust to actual time axis
figure;
plot(time, real_acc, 'b', 'LineWidth', 2); hold on;
plot(time(significant_bins), real_acc(significant_bins), 'ro', 'MarkerFaceColor','r');
yline(chance_level, '--k', 'Chance');
xlabel('Time (s)');
ylabel('Decoding Accuracy');
title(sprintf('Time-Resolved Decoding w/ Significance (Neuron %d)', sig_neuron_idx));
legend('Accuracy', 'Significant Bin');
grid on;

