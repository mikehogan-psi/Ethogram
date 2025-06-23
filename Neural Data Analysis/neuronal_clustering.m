%% Load coefficient matrix 

load("Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\logistic_regression_data\extinction_neuronal_cofficients.mat")
feature_vectors = glmm_output;
clearvars("glmm_output");


%%
% Detect neurons with any feature > 5 SD from mean
is_outlier = any(abs(feature_vectors) > 5, 2);  % returns logical vector of outliers

% Remove them
features_clean = feature_vectors(~is_outlier, :);
features_clean_z = zscore(features_clean);
%% Cluster data and plot clusters
K = 8; % Change number of clusters as needed

% Generate clusters
[cluster_idx, centroids] = kmeans(features_clean_z, K, 'Replicates', 50);


% Plot clusters
figure;
for k = 1:K
    colours = lines(k);
    plot(features_clean_z(cluster_idx==k,1), features_clean_z(cluster_idx==k, 2), ...
         '.', 'MarkerSize', 12, 'Color', colours(k,:));
    hold on
    plot(centroids(k,1), centroids(k,2), 'o', 'Color', colours(k,:));
end
hold off

%% Calculate means for each cluster

% Creates c x f matrix with mean of cofficients for each cluster, where c
% is the clusters and f is the features
cluster_means = grpstats(features_clean_z, cluster_idx); 
feature_labels = ["Intercept (Stim)", "TimeBin (Stim)", "Trial (Stim)",...
    "Intercept (Post)", "TimeBin (Post)", "Trial (Post)"];

figure;
imagesc(cluster_means);
xticks(1:8)
xticklabels(feature_labels)

ylabel('Cluster')
yticks(1:5)
colorbar

cells_per_cluster = zeros(1, k);
for c = 1:K
    current_cells_per_cluster = sum(cluster_idx == c);
    cells_per_cluster(1, c) = current_cells_per_cluster;
end

figure
bar(cells_per_cluster);
ylabel('Number of Cells')
xlabel('Cluster')
%%

figure;
for f = 1:length(feature_labels)
    subplot(2, 3, f)
    boxplot(features_clean_z(:, f), cluster_idx)
    title(feature_labels(f))
    xlabel('Cluster')
    ylabel('Coefficient (z-scored)')
end