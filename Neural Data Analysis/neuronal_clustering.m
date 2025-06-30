%% Load coefficient matrix 

load("Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\logistic_regression_data\extinction_neuronal_cofficients_poisson_resp_LOOM_nr.mat")
feature_vectors = glmm_output;
clearvars("glmm_output");


%% Remove outliers
%Save neuron numbers first
neuron_numbers = feature_vectors(:,end);

%Remove neuron numbers for analysis
feature_vectors = feature_vectors(:, 1:end-1);

% Detect neurons with any feature > n SD from mean
is_outlier = any(abs(feature_vectors) > 5, 2);  % Returns logical vector of outliers

% Remove them
features_clean = feature_vectors(~is_outlier, :);
features_clean_z = zscore(features_clean);


%% Cluster data and plot clusters
max_K = 10; % Change number of clusters as needed

avg_silhouette = zeros(max_K, 1);

% Find silhouette score for each possible cluster value in range 2:max_K 
for k = 2:max_K % silhouette doesn't make sense for K=1
    [cluster_idx] = kmeans(features_clean_z, k, 'Replicates', 50);
    s = silhouette(features_clean_z, cluster_idx);
    avg_silhouette(k) = mean(s);
end

% Plot silhouette scores
figure;
plot(2:max_K, avg_silhouette(2:end), '-o');
xlabel('Number of Clusters (K)');
ylabel('Mean Silhouette Score');

% Set K to value with highest silouette score
[~, K] = max(avg_silhouette);

[cluster_idx, centroids] = kmeans(features_clean_z, K, 'Replicates', 50);

% Plot clusters
figure;
for k = 1:K
    colours = lines(k);
    plot(features_clean_z(cluster_idx==k,1), features_clean_z(cluster_idx==k, 2), ...
         '.', 'MarkerSize', 12, 'Color', colours(k,:));
    hold on
    plot(centroids(k,1), centroids(k,2), 'x', 'MarkerSize', 10, 'Color', colours(k,:));
end
hold off
%% Save neuron number with cluster assignment
neuron_cluster_assignments = [neuron_numbers, cluster_idx];
%% Calculate means for each cluster

% Creates c x f matrix with mean of cofficients for each cluster, where c
% is the clusters and f is the features
cluster_means = grpstats(features_clean_z, cluster_idx); 
feature_labels = ["TimeBin (Stim)", "Trial (Stim)", "TimeBin (Post)", "Trial (Post)", "Freezing (Stim)", "Freezing (Post)"];

figure;
imagesc(cluster_means);
xticks(1:length(feature_labels))
xticklabels(feature_labels)

ylabel('Cluster')
yticks(1:K)
colorbar

cells_per_cluster = zeros(1, K);
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

%%
[coeff, score] = pca(features_clean_z);
figure;
scatter(score(:,1), score(:,2), 20, cluster_idx, 'filled');
%%
load("Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\logistic_regression_data\all_data_responsive_neurons.mat")

for row = 1:size(all_data_tables, 1)
    for col = 1:size(all_data_tables,2)
        current_table = all_data_tables{row, col};
        neuron_cluster = repelem(cluster_idx(row), 400)';
        current_table.ClusterID = neuron_cluster;
        all_data_tables{row, col} = current_table;
    end
end
%%
all_post_tables = vertcat(all_data_tables{:, 2});

% Combine post-stimulus tables
T_post = vertcat(all_data_tables{:, 2});

% Get unique cluster IDs
clusters = unique(T_post.ClusterID);

% Set up figure with appropriate number of subplots
num_clusters = length(clusters);
n_cols = ceil(sqrt(num_clusters));
n_rows = ceil(num_clusters / n_cols);

figure;

for i = 1:num_clusters
    c = clusters(i);
    T_c = T_post(T_post.ClusterID == c, :);

    % Compute mean firing rate per trial
    mean_FR = grpstats(T_c.Firing, T_c.Trial, 'mean');

    % Plot in subplot
    subplot(n_rows, n_cols, i);
    plot(1:length(mean_FR), mean_FR, 'LineWidth', 2);
    title("Cluster " + string(c));
    xlabel('Trial');
    xlim([1 20])
    ylabel('Mean FR (Hz)');
    ylim([0, max(mean_FR)*1.2]); % Optional: consistent y-axis scaling
end

sgtitle('Post-stimulus: Mean Firing Rate across Trials by Cluster');  % overall title

%%
% Combine post-stimulus tables
T_post = vertcat(all_data_tables{:, 2});

% Get unique cluster IDs
clusters = unique(T_post.ClusterID);

% Set up figure with appropriate number of subplots
num_clusters = length(clusters);
n_cols = ceil(sqrt(num_clusters));
n_rows = ceil(num_clusters / n_cols);

figure;

for i = 1:num_clusters
    c = clusters(i);
    T_c = T_post(T_post.ClusterID == c, :);

    % Compute mean firing rate per TimeBin
    mean_FR = grpstats(T_c.Firing, T_c.TimeBin, 'mean');

    % Plot in subplot
    subplot(n_rows, n_cols, i);
    plot(1:length(mean_FR), mean_FR, 'LineWidth', 2);
    title("Cluster " + string(c));
    xlabel('Time Bin (1s)');
    ylabel('Mean FR');
    xlim([1 20])
    ylim([0, max(mean_FR)*1.2]); % Optional: consistent y-axis
end

sgtitle('Post-stimulus: Mean Firing Rate across Time Bins by Cluster');
