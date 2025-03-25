% Calculate mean freezing time per trial for each mouse
post_stim_freezing = validated_freeze_matrix(:,post_stim_period,:);
mean_freezing_post_stim = squeeze(mean(validated_freeze_matrix, 2)); % Mean over frames
% Size: 40 x 30
%%
% Calculate the pairwise distances between mice based on total freezing
distances = pdist(mean_freezing_post_stim'); % Transpose to mice x trials

% Perform hierarchical clustering
Z = linkage(distances, 'ward'); % 'ward' minimizes within-cluster variance

% Plot the dendrogram
figure;
dendrogram(Z);
title('Hierarchical Clustering of Mice');
xlabel('Mouse ID');
ylabel('Euclidean Distance');
%%
num_clusters = 4; % Set the desired number of clusters
clusters = cluster(Z, 'maxclust', num_clusters);

% Visualise the clustering
disp('Cluster Assignments:');
disp(clusters);

disp('Vehicle Mice:')
for i = 1:2:length(clusters)
    disp(['Mouse ' num2str(i) ': ' 'cluster ' num2str(clusters(i))]);
end

disp('Psilocybin Mice:')
for i = 2:2:length(clusters)
    disp(['Mouse ' num2str(i) ': ' 'cluster ' num2str(clusters(i))]);
end
%% isolating data from clusters to determine what they mean

%creating logical masks for indexing
cluster_1_idx = clusters' == 1;
cluster_2_idx = clusters' == 2;
cluster_3_idx = clusters' == 3;
cluster_4_idx = clusters' == 4;

cluster_1_data = mean_freezing_post_stim(:, cluster_1_idx);
cluster_2_data = mean_freezing_post_stim(:, cluster_2_idx);
cluster_3_data = mean_freezing_post_stim(:, cluster_3_idx);
cluster_4_data = mean_freezing_post_stim(:, cluster_4_idx);

figure;
hold on;

% Combine data for boxplot
combined_data = [cluster_1_data(:); cluster_2_data(:); cluster_3_data(:); cluster_4_data(:)];
group_labels = [repmat(1, size(cluster_1_data(:))); ...
                repmat(2, size(cluster_2_data(:))); ...
                repmat(3, size(cluster_3_data(:))); ...
                repmat(4, size(cluster_4_data(:)))];

% Create the boxplot
boxplot(combined_data, group_labels, 'Labels', {'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'});
ylabel('Mean Proportion of Time Spent Frozen');
title('Comparison of Freezing Time Across Clusters (Post-Loom)');

hold off;
