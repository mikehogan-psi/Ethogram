% CLASSIFY_CELLS_MIKE
% Clusters neurons based on their GLM coefficients and compares
% cluster composition between treatment groups.
%
% Steps:
% 1) Load GLM results for each mouse
% 2) Keep only cells with pseudoR2 above threshold
% 3) Build coefficient matrix B (cells x coeffs), optionally drop intercept
% 4) Normalise / standardise coefficients
% 5) PCA for dimensionality reduction
% 6) Cluster cells (k-means) in PCA space
% 7) Visualise clusters with t-SNE
% 8) Compare cluster composition for treatment vs control

%----------------- USER INPUTS -----------------%
fpath      = 'D:\PhD 3rd Year\poisson_GLM_data_21_01_26';
mouse_num  = 2:9;

% TODO: replace with real treatment labels per mouse (0 = ctrl, 1 = treated)
mouse_treatment = [0,1,0,1,0,1,1,0];

TH_pR2     = 0.05;  % pseudoR2 threshold
dropIntercept = true; % set true if first coefficient is intercept
do_shuffle = false;
maxK       = 10;    % maximum K to evaluate for silhouette
rng(24);            % reproducibility for t-SNE and k-means
%------------------------------------------------%

B = [];               % coefficient matrix (cells x coeffs)
treatment_vec = [];   % treatment label per cell
mouse_vec     = [];   % mouse ID per cell (for later per-mouse stats)
cellID_vec   = [];   % local cell index within each mouse file


for im = 1:numel(mouse_num)
    mID = mouse_num(im);
    datafile = fullfile(fpath, ...
        ['mouse' num2str(mID) '_extinction_poisson_GLM_data_and_results.mat']);
    load(datafile);  % should load pseudoR2_test, weights_per_cell, etc.

    % Keep only cells with sufficient pseudoR2
    indok = find(pseudoR2_test > TH_pR2);   % indices into original cell list

    % Concatenate coefficients for these cells
      % [num_coeff x Ncells]
    temp = weights_per_cell(indok);                  % keep good cells only
    temp = horzcat(temp{:});  

    if dropIntercept
        temp = temp(2:end, :);
    end
    
    % Append as additional cells (transpose: cells x coeffs)
    B = cat(1, B, temp.');
   

    % Labels for these cells
    nOK = numel(indok);
    treatment_vec = cat(2, treatment_vec, mouse_treatment(im) * ones(1, nOK));
    mouse_vec     = cat(2, mouse_vec,     mID               * ones(1, nOK));
    cellID_vec = cat(2, cellID_vec, (indok - 1).');   % 0-based CellID


end



if do_shuffle
    nShuf = 200;
    bestSil_null = nan(nShuf,1);
    B0 = B;

    for s = 1:nShuf
        Bsh = B0;
        for i = 1:size(B0,1)
            Bsh(i,:) = B0(i, randperm(size(B0,2)));
        end
    
        % normalise exactly as real
        Bsh = Bsh ./ vecnorm(Bsh,2,2);
    
        [~, scoreS, lambdaS] = pca(Bsh, 'Centered', true);
        explVarS = cumsum(lambdaS)/sum(lambdaS);
        numPCS = find(explVarS>0.9,1,'first');
        Xs = scoreS(:,1:min(numPCS,10));
    
        Ks = 2:maxK;
        ms = nan(size(Ks));
        for ii = 1:numel(Ks)
            idxS = kmeans(Xs, Ks(ii), 'Replicates', 20, 'MaxIter', 1000, 'Display','off');
            ms(ii) = mean(silhouette(Xs, idxS));
        end
    
        bestSil_null(s) = max(ms);
    end
    
    fprintf('Row-shuffle null best silhouette: mean=%.3f, sd=%.3f\n', mean(bestSil_null), std(bestSil_null));

end



%----------------- NORMALISATION -----------------%
% Option 1: L2-normalise each cell's coefficient vector
for iCell = 1:size(B,1)
    nrm = norm(B(iCell,:));
    if nrm > 0
        B(iCell,:) = B(iCell,:) ./ nrm;
    end
end

% Option 2 (alternative): z-score each coefficient across cells
% B = zscore(B, 0, 1);

%----------------- PCA -----------------%
% Centre the data across cells
[coeff, score, lambda] = pca(B, 'Centered', true);

% Number of PCs explaining 90% variance
explVar = cumsum(lambda) / sum(lambda);
num_PC  = find(explVar > 0.9, 1, 'first');

% Use up to 10 PCs for clustering (or all if fewer)
num_PC_cluster = min(num_PC, 10);
X_cluster = score(:, 1:num_PC_cluster);

%----------------- CHOOSE K BY SILHOUETTE (OPTIONAL) -----------------%
Ks = 2:maxK;
mean_sil = nan(size(Ks));

for ii = 1:numel(Ks)
    % k-means for this K
    idx_tmp = kmeans(X_cluster, Ks(ii), ...
        'Replicates', 50, 'MaxIter', 1000, 'Display', 'off');
    
    % silhouette returns one value per point; take the mean
    s = silhouette(X_cluster, idx_tmp);
    mean_sil(ii) = mean(s);
end

if do_shuffle
    real_bestSil = max(mean_sil);
    p = mean(bestSil_null >= real_bestSil);
    fprintf('Real best silhouette = %.3f; p(null >= real) = %.5f\n', real_bestSil, p);
end

% Plot mean silhouette vs K
figure; 
plot(Ks, mean_sil, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Number of clusters (K)');
ylabel('Mean silhouette value');
title('Silhouette analysis for k-means clustering');
grid on;

% Pick K with highest mean silhouette
[~, bestIdx] = max(mean_sil);
K = Ks(bestIdx);
fprintf('Chosen K = %d (max mean silhouette = %.3f).\n', K, mean_sil(bestIdx));

% Final k-means with chosen K
idx = kmeans(X_cluster, K, 'Replicates', 100, 'MaxIter', 1000);

% Map clusters back to original cells
% One row per cell kept in B
cluster_assignments = table( ...
    mouse_vec.', ...
    cellID_vec.', ...
    treatment_vec.', ...
    idx, ...
    'VariableNames', {'MouseID','CellID','Treatment','ClusterID'});

% % Optional: save for later indexing
% save(fullfile(fpath, 'cluster_assignments_extinction.mat'), 'cluster_assignments');


%----------------- t-SNE FOR VISUALISATION -----------------%
B_tsne = tsne(X_cluster); % t-SNE on PCA space for plotting only

% Figure 1 - t-SNE clustering
figure; hold on;
colors = lines(K);
for k = 1:K
    scatter(B_tsne(idx==k,1), B_tsne(idx==k,2), 10, colors(k,:), 'filled');
end
xlabel('t-SNE1'); ylabel('t-SNE2');
title('t-SNE embedding coloured by k-means cluster');

% Figure 2 - mean coefficient pattern for each cluster
figure;
for k = 1:K
    subplot(ceil(K/2), 2, k); hold on;
    plot(mean(B(idx==k,:),1), 'LineWidth', 2, 'Marker', 'x');
    ylim([-1 1]); % adjust as needed
    xlim([1 10])
    xlabel('Coefficient #');
    ylabel('Coeff value (normalised)');
    title(sprintf('Cluster %d (n=%d)', k, sum(idx==k)));
end

%----------------- PER-MOUSE CLUSTER PROPORTIONS -----------------%
mouse_ids = unique(mouse_vec);
nMice     = numel(mouse_ids);

cluster_prop = nan(nMice, K);      % each row: one mouse; each column: cluster
treat_per_mouse = nan(nMice, 1);   % 0 = ctrl, 1 = treated

for iM = 1:nMice
    mID = mouse_ids(iM);

    % cells belonging to this mouse
    mask_m = (mouse_vec == mID);

    % cluster labels for this mouse's cells
    idx_m  = idx(mask_m);

    % total number of cells (after pR2 threshold etc.) for this mouse
    nCells_m = numel(idx_m);

    % proportions per cluster
    for k = 1:K
        cluster_prop(iM, k) = sum(idx_m == k) / nCells_m;
    end

    % treatment label for this mouse (mode across its cells)
    treat_per_mouse(iM) = mode(treatment_vec(mask_m));
end

is_ctrl    = (treat_per_mouse == 0);
is_treated = (treat_per_mouse == 1);

prop_ctrl    = cluster_prop(is_ctrl, :);    % [nCtrl x K]
prop_treated = cluster_prop(is_treated, :); % [nTreat x K]

nCtrl   = sum(is_ctrl);
nTreat  = sum(is_treated);

mean_ctrl = mean(prop_ctrl, 1);
sem_ctrl  = std(prop_ctrl, [], 1) / sqrt(nCtrl);

mean_treated = mean(prop_treated, 1);
sem_treated  = std(prop_treated, [], 1) / sqrt(nTreat);

figure; hold on;

x = 1:K;
bar_width = 0.4;

% bars: ctrl and treated
b1 = bar(x - bar_width/2, mean_ctrl,  bar_width);
b2 = bar(x + bar_width/2, mean_treated, bar_width);

% error bars
errorbar(x - bar_width/2, mean_ctrl,  sem_ctrl, 'k', 'linestyle', 'none');
errorbar(x + bar_width/2, mean_treated, sem_treated, 'k', 'linestyle', 'none');

xlabel('Cluster ID');
ylabel('Proportion of cells per mouse');
legend({'Ctrl', 'Treated'}, 'Location', 'best');
title('Per-mouse cluster proportions');
ylim([0 1]); % proportions

p_vals = nan(1, K);

for k = 1:K
    x_ctrl    = prop_ctrl(:, k);
    x_treated = prop_treated(:, k);

    % use t-test if you’re happy with normality;
    % use ranksum if n is tiny / distributions dodgy.
    % [~, p_vals(k)] = ttest2(x_ctrl, x_treated);
    p_vals(k) = ranksum(x_ctrl, x_treated);
end

disp('p-values per cluster (ctrl vs treated, per-mouse proportions):');
disp(p_vals);


