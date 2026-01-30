% Compares mean coefficients between extinction and retention per cluster
% after retention data has been projected onto clusters formed based only
% on extinction data (shows how well they fit existing clusters)
%
% Loaded files: cluster assignments for extinction and retention +
% coefficient matrices for extinction and retention

load("D:\PhD 3rd Year\renewal_glm_data\cluster_assignments_retention_to_extinction.mat")
load("D:\PhD 3rd Year\poisson_GLM_data_21_01_26\cluster_assignments_extinction.mat")
load("D:\PhD 3rd Year\poisson_GLM_data_21_01_26\extinction_coefficients.mat")
load("D:\PhD 3rd Year\renewal_glm_data\renewal_coefficients.mat")

num_weights = size(B_ret, 2);
num_clusters = max(cluster_assignments_retention.ClusterID);

mean_weights_retention = zeros(num_weights, num_clusters);


for cluster = 1:num_clusters
    % Filter coefficients matrix for current cluster data only
    current_clu_mask = cluster_assignments_retention.ClusterID == cluster;
    clu_coefs = B_ret(current_clu_mask, :);
    % For each weight, calculate the mean 
    for weight = 1:num_weights
        current_col = clu_coefs(:, weight);
        mean_current_col = mean(current_col,'omitnan');
        mean_weights_retention(weight, cluster) = mean_current_col;
    end
end

%%
load("D:\PhD 3rd Year\poisson_GLM_data_21_01_26\cluster_assignments_extinction.mat")

num_weights = size(B_ret, 2);
num_clusters = max(cluster_assignments.ClusterID);

mean_weights_extinction = zeros(num_weights, num_clusters);


for cluster = 1:num_clusters
    % Filter coefficients matrix for current cluster data only
    current_clu_mask = cluster_assignments.ClusterID == cluster;
    clu_coefs = B(current_clu_mask, :);
    % For each weight, calculate the mean 
    for weight = 1:num_weights
        current_col = clu_coefs(:, weight);
        mean_current_col = mean(current_col,'omitnan');
        mean_weights_extinction(weight, cluster) = mean_current_col;
    end
end

%%
% ---- inputs you should have ----
% mean_ext: [num_weights x num_clusters]
% mean_ret: [num_weights x num_clusters]
% weight_names: string/cellstr [num_weights x 1] (e.g., your GLM covariate names)

weight_names = coef_names;
mean_ext = mean_weights_extinction;
mean_ret = mean_weights_retention;

figure('Color','w');

nRows = ceil(num_clusters/2);
nCols = 2;

for c = 1:num_clusters
    subplot(nRows, nCols, c); hold on;

    x = 1:num_weights;

    plot(x, mean_ext(:,c), 'LineWidth', 2);        % extinction (template)
    plot(x, mean_ret(:,c), '--', 'LineWidth', 1.8); % retention-assigned

    yline(0, ':');

    xticks(x);
    xticklabels(weight_names);
    xtickangle(45);

    title(sprintf('Cluster %d', c));
    ylabel('Mean coefficient');
    xlim([0.5 num_weights+0.5]);

    ax = gca;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.LineWidth = 1.2;

    if c == 1
        legend({'Ext (template)','Ret (assigned)'}, 'Box','off', 'Location','best');
    end
end

sgtitle('Mean GLM coefficients per cluster: Extinction vs Retention-assigned');
