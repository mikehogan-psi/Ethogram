% PSTH_ACROSS_BINS
% Generates mean PSTHs across bins for all neurons and for individual
% clusters, as assigned by classify_cells_temp. Set do_normalise to
% true to perform a baseline subtraction on data (mean pre-stimulus firing)

master_directory = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)";

session= 'Extinction';

mice_to_analyse = [2 3 4 5 6 7 8 9];

num_mice = max(mice_to_analyse);

trial_type = 0; % 1 = loom, 0 = flash

received_psi = [3 5 7 8];
received_veh = [2 4 6 9];
do_normalise = false;


if strcmp(session, 'Extinction')
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ext\cluster_assignments_extinction.mat")
    cluster_assignments = cluster_assignments_extinction;
elseif strcmp(session, 'Renewal')
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ren\cluster_assignments_renewal_to_extinction.mat")
    cluster_assignments = cluster_assignments_renewal;
end

%% Load processed neural and behavioural data

mouse_files = dir(fullfile(master_directory, 'Mouse*'));

mouse_paths = cell(numel(mouse_files), 1);

for mouse = mice_to_analyse
    current_mouse_folder = mouse_files(mouse).folder;
    current_mouse_name = mouse_files(mouse).name;
    current_mouse_path = fullfile(current_mouse_folder, current_mouse_name,...
        session, 'Combined Data');
    mouse_paths{mouse} = current_mouse_path;
end

combined_data_matrices = cell(numel(mouse_files), 1);

for mouse = mice_to_analyse
    current_file_to_load = dir(fullfile(mouse_paths{mouse}, '*GLM_data*'));
    current_data = load(fullfile(current_file_to_load.folder, current_file_to_load.name));
    combined_data_matrices{mouse} = current_data;
end


%%
psth_data = cell(num_mice, 1);


for mouse = mice_to_analyse
    current_data = combined_data_matrices{mouse};
    current_data = current_data.X_table_full;
    spiking_of_interest_idx = current_data.TrialIdentifier == trial_type;
    spiking_of_interest.SpikeCount = current_data.SpikeCount(spiking_of_interest_idx);
    spiking_of_interest.TimeBin = repmat(1:66, 1, length(spiking_of_interest.SpikeCount)/66)';
    spiking_of_interest.CellID = current_data.CellID(spiking_of_interest_idx);
    psth_data{mouse} = spiking_of_interest;
end

num_bins = 66;
pre_stim_bins = 1:19;

if do_normalise

    for mouse = mice_to_analyse
        current_data = psth_data{mouse};
        num_cells = length(unique(current_data.CellID));
    
        for neuron = 1:num_cells
            current_neuron = neuron-1;
            current_neuron_idx = current_data.CellID == current_neuron;
            current_spikes = current_data.SpikeCount(current_neuron_idx);
            current_bins = current_data.TimeBin(current_neuron_idx);
    
            pre_idx = current_bins >= pre_stim_bins(1) & current_bins <= pre_stim_bins(end);
            baseline = mean(current_spikes(pre_idx));
            norm_spikes = current_spikes - baseline;

            current_data.SpikeCount(current_neuron_idx) = norm_spikes;
        end
    
        psth_data{mouse} = current_data;
    
    end

end
 

avg_spiking_per_bin = nan(num_mice, num_bins);

for mouse = mice_to_analyse
    current_data = psth_data{mouse};
    for bin = 1:num_bins
        current_bin_data_idx = current_data.TimeBin == bin;
        current_bin_data = mean(current_data.SpikeCount(current_bin_data_idx));
        avg_spiking_per_bin(mouse, bin) = current_bin_data;
    end
end  
        

%%
time_bin_width = 0.5;              % s per bin
t = (0:num_bins-1) * time_bin_width;

% Group means 
psi_spiking = avg_spiking_per_bin(received_psi,:);
veh_spiking = avg_spiking_per_bin(received_veh,:);

mean_psi_spiking = mean(psi_spiking, 1);
mean_veh_spiking = mean(veh_spiking, 1);

% SEM across mice
sem_psi = std(psi_spiking, 0, 1) ./ sqrt(size(psi_spiking,1));
sem_veh = std(veh_spiking, 0, 1) ./ sqrt(size(veh_spiking,1));

%% Simple line plot with errorbars
veh_col = [0 0.4470 0.7410];   % blue
psi_col = [1 0 0];             % red

veh_upper = mean_veh_spiking + sem_veh;
veh_lower = mean_veh_spiking - sem_veh;

psi_upper = mean_psi_spiking + sem_psi;
psi_lower = mean_psi_spiking - sem_psi;

figure('Color','w'); hold on;

% Vehicle shaded SEM
fill([t fliplr(t)], [veh_upper fliplr(veh_lower)], veh_col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Vehicle mean line
plot(t, mean_veh_spiking, 'Color', veh_col, 'LineWidth', 2);

% Psilocybin shaded SEM
fill([t fliplr(t)], [psi_upper fliplr(psi_lower)], psi_col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Psilocybin mean line
plot(t, mean_psi_spiking, 'Color', psi_col, 'LineWidth', 2);

xline(10, 'k--', 'LineWidth', 1.5);
xline(13.3, 'k--', 'LineWidth', 1.5)

xlabel('Time (s)');
if do_normalise
    ylabel('Mean spike count per bin (normalised to pre-stim)');
    ylim([-1 2])
else
    ylabel('Mean spike count per bin');
end
title(['PSTH: Vehicle vs Psilocybin Loom ' session]);
legend({'Veh SEM','Vehicle','Psi SEM','Psilocybin'}, 'Location','best');


ax = gca;
ax.LineWidth = 1.5;
ax.TickDir   = 'out';
ax.Box       = 'off';
xlim([0 33])



%%

post_idx = t >= 15 & t <= 30;   % logical index for post-stim bins

veh_post = mean(veh_spiking(:, post_idx), 2);   % mean over bins, per mouse
psi_post = mean(psi_spiking(:, post_idx), 2);

[~,p,~,stats] = ttest2(veh_post, psi_post);
fprintf('Post-stim t-test: p = %.4f, t(%d) = %.3f\n', ...
        p, stats.df, stats.tstat);

%% Plot PSTHs for individual clusters

num_clusters = max(cluster_assignments.ClusterID);
sorted_clusters = cell(num_mice, num_clusters);
num_bins = 66;
pre_stim_bins = 1:19;

    
for mouse = mice_to_analyse

    current_data = combined_data_matrices{mouse}.X_table_full;
    for cluster = 1:num_clusters
        spiking_temp = [];
        timebin_holder = [];
        neuron_ID_holder = [];
        current_cluster_idx = cluster_assignments.ClusterID == cluster...
            & cluster_assignments.MouseID == mouse;
        current_cluster_neurons = cluster_assignments.CellID(current_cluster_idx);
        for neuron_idx = 1:length(current_cluster_neurons)
            neuron = current_cluster_neurons(neuron_idx);
            current_neuron_idx = current_data.CellID == neuron & current_data.TrialIdentifier == trial_type;

            spiking_temp = [spiking_temp; current_data.SpikeCount(current_neuron_idx)];
            timebin_holder = [timebin_holder; repmat(1:66, 1, sum(current_neuron_idx)/num_bins)'];
            neuron_ID_holder = [neuron_ID_holder; repelem(neuron, sum(current_neuron_idx), 1)];
        
        end        
        sorted_clusters{mouse, cluster} = table(spiking_temp, timebin_holder, neuron_ID_holder,...
           'VariableNames', {'SpikeCount', 'TimeBin', 'CellID'});
    end
    
end



if do_normalise
    for mouse = mice_to_analyse
        for cluster = 1:num_clusters
            current_data = sorted_clusters{mouse, cluster};            
            cellIDs = unique(current_data.CellID);
            num_cells = length(cellIDs);
            for neuron_idx = 1:num_cells
                neuron = cellIDs(neuron_idx);
                all_neuron_idx = current_data.CellID == neuron;
                current_spikes = current_data.SpikeCount(all_neuron_idx);
                current_bins = current_data.TimeBin(all_neuron_idx);
    
                pre_idx = current_bins >= pre_stim_bins(1) & current_bins <= pre_stim_bins(end);
                baseline = mean(current_spikes(pre_idx));
                norm_spikes = current_spikes - baseline;

                current_data.SpikeCount(all_neuron_idx) = norm_spikes;
            end
        
        sorted_clusters{mouse, cluster} = current_data;

        end  

    end
end


avg_spiking_per_bin_clusters = cell(num_clusters, 1);

for cluster = 1:num_clusters
    current_cluster_data = sorted_clusters(:, cluster);

    avg_spiking_per_bin = nan(num_mice, num_bins);
    
    for mouse = mice_to_analyse
        current_data = current_cluster_data{mouse};
        for bin = 1:num_bins
            current_bin_data_idx = current_data.TimeBin == bin;
            current_bin_data = mean(current_data.SpikeCount(current_bin_data_idx));
            avg_spiking_per_bin(mouse, bin) = current_bin_data;
        end
    end
    
    avg_spiking_per_bin_clusters{cluster} = avg_spiking_per_bin;
end


psi_clusters = cell(num_clusters, 1);
veh_clusters = cell(num_clusters, 1);

for cluster = 1:num_clusters
    current_cluster = avg_spiking_per_bin_clusters{cluster};
    psi_clusters_temp = current_cluster(received_psi, :);
    veh_clusters_temp = current_cluster(received_veh, :);
    psi_clusters{cluster} = psi_clusters_temp;
    veh_clusters{cluster} = veh_clusters_temp;
end

%%
time_bin_width = 0.5;              % s per bin
t = (0:num_bins-1) * time_bin_width;

veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

figure('Color','w');

for cluster = 1:num_clusters
    subplot(2,2,cluster); hold on;  % 4x2 grid for 8 clusters

    psi_mat = psi_clusters{cluster};   % [nPsi x num_bins]
    veh_mat = veh_clusters{cluster};   % [nVeh x num_bins]

    % skip if cluster empty for both groups
    if isempty(psi_mat) && isempty(veh_mat)
        title(sprintf('Cluster %d (no data)', cluster));
        axis off;
        continue;
    end

    % handle possible NaNs (e.g. empty mice)
    if ~isempty(veh_mat)
        valid_veh = ~all(isnan(veh_mat), 2);
        veh_mat   = veh_mat(valid_veh, :);
        mean_veh  = mean(veh_mat, 1, 'omitnan');
        sem_veh   = std(veh_mat, 0, 1, 'omitnan') ./ sqrt(size(veh_mat,1));
    else
        mean_veh = nan(1, num_bins);
        sem_veh  = nan(1, num_bins);
    end

    if ~isempty(psi_mat)
        valid_psi = ~all(isnan(psi_mat), 2);
        psi_mat   = psi_mat(valid_psi, :);
        mean_psi  = mean(psi_mat, 1, 'omitnan');
        sem_psi   = std(psi_mat, 0, 1, 'omitnan') ./ sqrt(size(psi_mat,1));
    else
        mean_psi = nan(1, num_bins);
        sem_psi  = nan(1, num_bins);
    end

        % --- Per-mouse lines (overlay) ---
    % Vehicle mice (thin lines)
    if ~isempty(veh_mat)
        for i = 1:size(veh_mat,1)
            plot(t, veh_mat(i,:), 'Color', [veh_col 0.3], 'LineWidth', 1);
        end
    end

    % Psilocybin mice (thin lines)
    if ~isempty(psi_mat)
        for i = 1:size(psi_mat,1)
            plot(t, psi_mat(i,:), 'Color', [psi_col 0.3], 'LineWidth', 1);
        end
    end


    % shaded SEM – vehicle
    veh_upper = mean_veh + sem_veh;
    veh_lower = mean_veh - sem_veh;
    fill([t fliplr(t)], [veh_upper fliplr(veh_lower)], veh_col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(t, mean_veh, 'Color', veh_col, 'LineWidth', 1.5);

    % shaded SEM – psilocybin
    psi_upper = mean_psi + sem_psi;
    psi_lower = mean_psi - sem_psi;
    fill([t fliplr(t)], [psi_upper fliplr(psi_lower)], psi_col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(t, mean_psi, 'Color', psi_col, 'LineWidth', 1.5);



    % cosmetics
    title(sprintf('Cluster %d', cluster));
    xlabel('Time (s)');
    % if do_normalise
    %     ylabel('Normalised Mean Spike Count/Bin');
    % else
        ylabel('Mean Spike Count/Bin');
    % end
    xlim([t(1) t(end)]);
    ax = gca;
    ax.LineWidth = 1;
    ax.TickDir   = 'out';
    ax.Box       = 'off';
    xline(10, 'k--', 'LineWidth', 1.5);
    xline(13.3, 'k--', 'LineWidth', 1.5)
    % ylim([0 15])

    if cluster == 1
        legend({'Veh SEM','Vehicle','Psi SEM','Psilocybin'}, ...
               'Location','best', 'Box','off');
    end
end

%% Statistical tests: compare baselines between treatment groups
baseline_bins = 1:19;

p_vals = nan(num_clusters,1);
baseline_rates_all = nan(9, num_clusters);  % 9 mice, row 1 empty

for cluster = 1:num_clusters
    baseline_rates = mean(avg_spiking_per_bin_clusters{cluster}(:, baseline_bins), 2, 'omitnan'); % [9x1]
    baseline_rates_all(:,cluster) = baseline_rates;

    x_ctrl = baseline_rates(received_veh);
    x_psi  = baseline_rates(received_psi);

    p_vals(cluster) = ranksum(x_ctrl, x_psi);
end

disp(p_vals)
