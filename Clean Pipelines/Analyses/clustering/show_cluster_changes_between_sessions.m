% Load assignment tables
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ren\cluster_assignments_renewal_to_extinction.mat")
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ext\cluster_assignments_extinction.mat")

mouse_ids = unique(cluster_assignments_renewal.MouseID);
cluster_ids = unique(cluster_assignments_renewal.ClusterID);

num_mice = numel(mouse_ids);
num_clusters = numel(cluster_ids);

ext_clu_number_per_mouse = zeros(num_mice + 1, num_clusters);

for mouse_idx = mouse_ids(:)'
    % Process data for each mouse
    current_data = cluster_assignments_extinction(cluster_assignments_extinction.MouseID == mouse_idx, :);
    % Perform analysis or calculations on current_data here
    for cluster_idx = cluster_ids(:)'
        cluster_data = current_data(current_data.ClusterID == cluster_idx, :);
        % Perform further analysis on cluster_data
        num_current_cluster = height(cluster_data);
        ext_clu_number_per_mouse(mouse_idx, cluster_idx) = num_current_cluster;
    end
end


ren_clu_number_per_mouse = zeros(num_mice + 1, num_clusters);

for mouse_idx = mouse_ids(:)'
    % Process data for each mouse
    current_data = cluster_assignments_renewal(cluster_assignments_renewal.MouseID == mouse_idx, :);
    % Perform analysis or calculations on current_data here
    for cluster_idx = cluster_ids(:)'
        cluster_data = current_data(current_data.ClusterID == cluster_idx, :);
        % Perform further analysis on cluster_data
        num_current_cluster = height(cluster_data);
        ren_clu_number_per_mouse(mouse_idx, cluster_idx) = num_current_cluster;
    end
end

%%
delta = ren_clu_number_per_mouse - ext_clu_number_per_mouse;  % [nMice x nClusters]

figure;
imagesc(delta);
colorbar;
xlabel('Cluster');
ylabel('Mouse');
xticks(1:numel(cluster_ids)); xticklabels(string(cluster_ids));
yticks(1:numel(mouse_ids));   yticklabels(string(mouse_ids));
title('Change in # cells per cluster (Renewal - Extinction)');

%%
% Drop row 1 (MouseID==1 had no data)
ext_counts = ext_clu_number_per_mouse(2:end, :);
ren_counts = ren_clu_number_per_mouse(2:end, :);
mouse_ids_plot = (2:size(ext_clu_number_per_mouse,1))';  % 2..9

% Convert to proportions per mouse (rows sum to 1)
ext_tot = sum(ext_counts, 2);
ren_tot = sum(ren_counts, 2);

ext_plot = ext_counts ./ ext_tot;
ren_plot = ren_counts ./ ren_tot;

% If any mouse has 0 total cells, avoid Inf/NaN explosions
ext_plot(ext_tot == 0, :) = NaN;
ren_plot(ren_tot == 0, :) = NaN;

% Consistent colours across subplots
C = lines(numel(mouse_ids_plot));

nClu  = numel(cluster_ids);
nCols = ceil(sqrt(nClu));
nRows = ceil(nClu / nCols);

figure;
for c = 1:nClu
    subplot(nRows, nCols, c); hold on;

    for m = 1:numel(mouse_ids_plot)
        y = [ext_plot(m,c) ren_plot(m,c)];
        if all(isnan(y)), continue; end

        plot([1 2], y, '-o', ...
            'Color', C(m,:), ...
            'MarkerFaceColor', C(m,:));
    end

    xlim([0.8 2.2]);
    ylim([0 1]);
    xticks([1 2]);
    xticklabels({'Ext','Ren'});
    ylabel('Proportion of cells');
    title(sprintf('Cluster %s', string(cluster_ids(c))));
end
sgtitle('Per-mouse change by cluster (proportions)');
