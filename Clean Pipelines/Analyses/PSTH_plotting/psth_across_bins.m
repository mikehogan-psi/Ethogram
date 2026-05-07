% PSTH_ACROSS_BINS
% -------------------------------------------------------------------------
% This script generates peri-stimulus time histograms (PSTHs) from binned
% spike-count data. It can plot:
%
%   1) Mean PSTHs across all neurons, averaged per mouse, then compared
%      between treatment groups.
%
%   2) Mean PSTHs separately for each GLM-derived cluster.
%
%   3) Trial-block PSTHs, showing how cluster responses evolve across
%      blocks of trials.
%
% Optional baseline subtraction can be applied using the mean pre-stimulus
% firing rate for each neuron.
%
% Assumptions:
%   - SpikeCount has already been binned into 0.5 s bins.
%   - Each trial has 66 bins.
%   - TrialIdentifier == 1 corresponds to loom trials.
%   - TrialIdentifier == 0 corresponds to flash trials.
%   - CellIDs are zero-indexed in the GLM table.
%   - mouse_files indexing matches the mouse numbers in mice_to_analyse.
% -------------------------------------------------------------------------


%% PSTH settings

mice_to_analyse = [2 3 4 5 6 7 8 9];

trial_type = 1; % 1 = loom, 0 = flash

received_psi = [3 5 7 8];
received_veh = [2 4 6 9];

do_normalise = true;

num_bins = 66;
pre_stim_bins = 1:19;

time_bin_width = 0.5;
t = (0:num_bins-1) * time_bin_width;

num_mice = numel(hist_and_spiking_matrices);

%% Extract PSTH data for the selected trial type

% This cell array stores one table-like structure per mouse containing:
%   - SpikeCount
%   - TimeBin
%   - CellID
%% Extract PSTH data from hist_and_spiking_matrices

psth_data = cell(num_mice, 1);

for mouse = mice_to_analyse

    current_data = hist_and_spiking_matrices{mouse};

    if isempty(current_data) || height(current_data) == 0
        warning(['Mouse ' num2str(mouse) ' has no hist/spiking data'])
        continue
    end

    % Keep only the trial type of interest
    current_data = current_data(current_data.TrialIdentifier == trial_type, :);

    % Reconstruct TimeBin
    % This assumes rows are ordered by cell/trial/bin and each trial has 66 bins.
    if mod(height(current_data), num_bins) ~= 0
        warning(['Mouse ' num2str(mouse) ': row count not divisible by num_bins'])
        continue
    end

    current_data.TimeBin = repmat((1:num_bins)', height(current_data)/num_bins, 1);

    psth_data{mouse} = current_data;

end

%% Optional baseline subtraction

if do_normalise

    for mouse = mice_to_analyse
        current_data = psth_data{mouse};

        if isempty(current_data) || height(current_data) == 0
            continue
        end

        % This assumes CellIDs run from 0 to num_cells-1.
        % If some CellIDs are missing or non-contiguous, using unique(CellID)
        % would be safer.
        cellIDs = unique(current_data.CellID);

        for neuron_idx = 1:numel(cellIDs)
            current_neuron = cellIDs(neuron_idx);
            current_neuron_idx = current_data.CellID == current_neuron;

            current_spikes = current_data.SpikeCount(current_neuron_idx);
            current_bins = current_data.TimeBin(current_neuron_idx);

            % Identify pre-stimulus bins for this neuron
            pre_idx = current_bins >= pre_stim_bins(1) & current_bins <= pre_stim_bins(end);

            % Mean baseline firing across all pre-stim bins/trials
            baseline = mean(current_spikes(pre_idx));

            % Subtract baseline from all bins for this neuron
            norm_spikes = current_spikes - baseline;
            current_data.SpikeCount(current_neuron_idx) = norm_spikes;
        end

        psth_data{mouse} = current_data;
    end
end

%% Average spike count per time bin, per mouse

% Rows are mice, columns are time bins.
% Row 1 may remain NaN if Mouse 1 is not included.
avg_spiking_per_bin = nan(num_mice, num_bins);

for mouse = mice_to_analyse
    current_data = psth_data{mouse};

    if isempty(current_data) || ~istable(current_data) || height(current_data) == 0
        continue
    end

    for bin = 1:num_bins
        current_bin_data_idx = current_data.TimeBin == bin;
        current_bin_data = mean(current_data.SpikeCount(current_bin_data_idx));
        avg_spiking_per_bin(mouse, bin) = current_bin_data;
    end
end

%% Compute group means and SEMs

time_bin_width = 0.5;              % seconds per bin
t = (0:num_bins-1) * time_bin_width;

% Extract mouse-level PSTHs for each treatment group
psi_spiking = avg_spiking_per_bin(received_psi, :);
veh_spiking = avg_spiking_per_bin(received_veh, :);

% Remove mice with no data
valid_psi = ~all(isnan(psi_spiking), 2);
valid_veh = ~all(isnan(veh_spiking), 2);

psi_spiking_valid = psi_spiking(valid_psi, :);
veh_spiking_valid = veh_spiking(valid_veh, :);

% Means across valid mice only
mean_psi_spiking = mean(psi_spiking_valid, 1, 'omitnan');
mean_veh_spiking = mean(veh_spiking_valid, 1, 'omitnan');

% SEM across valid mice only
sem_psi = std(psi_spiking_valid, 0, 1, 'omitnan') ./ sqrt(size(psi_spiking_valid, 1));
sem_veh = std(veh_spiking_valid, 0, 1, 'omitnan') ./ sqrt(size(veh_spiking_valid, 1));
%% Plot overall vehicle vs psilocybin PSTH

% Global plotting defaults
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultTextFontWeight', 'bold');
set(groot, 'defaultAxesLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 12);

% Plot colours
veh_col = [0 0.4470 0.7410];   % blue
psi_col = [1 0 0];             % red

% SEM bounds for shaded error regions
veh_upper = mean_veh_spiking + sem_veh;
veh_lower = mean_veh_spiking - sem_veh;

psi_upper = mean_psi_spiking + sem_psi;
psi_lower = mean_psi_spiking - sem_psi;

figure('Color', 'w'); hold on;

% Vehicle shaded SEM
fill([t fliplr(t)], [veh_upper fliplr(veh_lower)], veh_col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Vehicle mean PSTH
plot(t, mean_veh_spiking, 'Color', veh_col, 'LineWidth', 2);

% Psilocybin shaded SEM
fill([t fliplr(t)], [psi_upper fliplr(psi_lower)], psi_col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Psilocybin mean PSTH
plot(t, mean_psi_spiking, 'Color', psi_col, 'LineWidth', 2);

% Stimulus onset/offset markers.
% These values assume stimulus starts around 10 s and ends around 13.3 s.
xline(10, 'k--', 'LineWidth', 1.5);
xline(13.3, 'k--', 'LineWidth', 1.5);

xlabel('Time (s)');

if do_normalise
    ylabel('Mean spike count per bin (normalised to pre-stim)');
    ylim([-1 2]);
else
    ylabel('Mean spike count per bin');
end

title(['PSTH: Vehicle vs Psilocybin Loom ' session]);
legend({'Veh SEM', 'Vehicle', 'Psi SEM', 'Psilocybin'}, 'Location', 'best');

ax = gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.Box = 'off';
xlim([0 33]);

%% Simple post-stimulus group comparison

% Define a post-stimulus window for statistical comparison.
% Here this is 15 to 30 s, meaning the test is not on the immediate stimulus
% window, but on later post-stimulus activity.
post_idx = t >= 15 & t <= 30;

% Mean over post-stim bins, producing one value per mouse
veh_post = mean(veh_spiking_valid(:, post_idx), 2, 'omitnan');
psi_post = mean(psi_spiking_valid(:, post_idx), 2, 'omitnan');

[~, p, ~, stats] = ttest2(veh_post, psi_post);

fprintf('Post-stim t-test: p = %.4f, t(%d) = %.3f\n', ...
        p, stats.df, stats.tstat);

%% Sort spikes by ClusterID already stored in hist_and_spiking_matrices

num_clusters = max(cluster_assignments.ClusterID);

sorted_clusters = cell(num_mice, num_clusters);

for mouse = mice_to_analyse

    current_data = hist_and_spiking_matrices{mouse};

    if isempty(current_data) || height(current_data) == 0
        continue
    end

    % Keep only selected trial type
    current_data = current_data(current_data.TrialIdentifier == trial_type, :);

    % Reconstruct TimeBin
    if mod(height(current_data), num_bins) ~= 0
        warning(['Mouse ' num2str(mouse) ': row count not divisible by num_bins'])
        continue
    end

    current_data.TimeBin = repmat((1:num_bins)', height(current_data)/num_bins, 1);

    for cluster = 1:num_clusters

        cluster_idx = current_data.ClusterID == cluster;

        sorted_clusters{mouse, cluster} = current_data(cluster_idx, :);

    end
end

%% Optional baseline subtraction within each cluster

if do_normalise
    for mouse = mice_to_analyse
        for cluster = 1:num_clusters
            current_data = sorted_clusters{mouse, cluster};

            if isempty(current_data) || ~istable(current_data) || height(current_data) == 0
                continue
            end

            cellIDs = unique(current_data.CellID);
            num_cells = length(cellIDs);

            for neuron_idx = 1:num_cells
                neuron = cellIDs(neuron_idx);
                all_neuron_idx = current_data.CellID == neuron;

                current_spikes = current_data.SpikeCount(all_neuron_idx);
                current_bins = current_data.TimeBin(all_neuron_idx);

                % Compute neuron-specific baseline using pre-stimulus bins
                pre_idx = current_bins >= pre_stim_bins(1) & current_bins <= pre_stim_bins(end);
                baseline = mean(current_spikes(pre_idx));

                % Subtract baseline from that neuron's full PSTH
                norm_spikes = current_spikes - baseline;
                current_data.SpikeCount(all_neuron_idx) = norm_spikes;
            end

            sorted_clusters{mouse, cluster} = current_data;
        end
    end
end

%% Add trial numbers and trial-block labels

% This section reconstructs trial identity for each neuron based on row
% order. It assumes that rows for each CellID are ordered sequentially by
% trial and bin.
%
% TrialBlock groups trials into blocks of 5:
%   Block 1 = trials 1-5
%   Block 2 = trials 6-10
%   Block 3 = trials 11-15
%   Block 4 = trials 16-20

num_bins = 66;
trials_per_block = 5;

for mouse = mice_to_analyse
    for cluster = 1:num_clusters
        T = sorted_clusters{mouse, cluster};

        % Skip empty mouse/cluster combinations
        if isempty(T) || height(T) == 0
            sorted_clusters{mouse, cluster} = T;
            continue
        end

        % Preallocate trial labels
        T.TrialNum = nan(height(T), 1);
        T.TrialBlock = nan(height(T), 1);

        cellIDs = unique(T.CellID);

        for c = 1:numel(cellIDs)
            cid = cellIDs(c);

            % Indices for this neuron. find() preserves row order.
            idx = find(T.CellID == cid);

            n = numel(idx);
            nTrials = n / num_bins;

            % Check that row count is compatible with complete trials
            if abs(nTrials - round(nTrials)) > 1e-9
                warning('Mouse %d Cluster %d CellID %d: rows (%d) not divisible by %d bins.', ...
                        mouse, cluster, cid, n, num_bins);
                continue
            end

            nTrials = round(nTrials);

            % Assign trial numbers and block numbers for this neuron
            T.TrialNum(idx) = repelem((1:nTrials)', num_bins, 1);
            T.TrialBlock(idx) = ceil(T.TrialNum(idx) / trials_per_block);
        end

        sorted_clusters{mouse, cluster} = T;
    end
end

%% Average spike count per bin for each cluster and mouse

% avg_spiking_per_bin_clusters{cluster} is a matrix:
%   rows = mice
%   columns = time bins
avg_spiking_per_bin_clusters = cell(num_clusters, 1);

for cluster = 1:num_clusters
    current_cluster_data = sorted_clusters(:, cluster);

    avg_spiking_per_bin = nan(num_mice, num_bins);

    for mouse = mice_to_analyse
        current_data = current_cluster_data{mouse};

        if isempty(current_data) || ~istable(current_data) || height(current_data) == 0
            continue
        end

        for bin = 1:num_bins
            current_bin_data_idx = current_data.TimeBin == bin;
            current_bin_data = mean(current_data.SpikeCount(current_bin_data_idx));
            avg_spiking_per_bin(mouse, bin) = current_bin_data;
        end
    end

    avg_spiking_per_bin_clusters{cluster} = avg_spiking_per_bin;
end

%% Split cluster PSTHs by treatment group

psi_clusters = cell(num_clusters, 1);
veh_clusters = cell(num_clusters, 1);

for cluster = 1:num_clusters
    current_cluster = avg_spiking_per_bin_clusters{cluster};

    psi_clusters_temp = current_cluster(received_psi, :);
    veh_clusters_temp = current_cluster(received_veh, :);

    psi_clusters{cluster} = psi_clusters_temp;
    veh_clusters{cluster} = veh_clusters_temp;
end

%% Plot PSTHs separately for each cluster

time_bin_width = 0.5;
t = (0:num_bins-1) * time_bin_width;

veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

figure('Color', 'w');

for cluster = 1:num_clusters
    subplot(2, 1, cluster); hold on;

    psi_mat = psi_clusters{cluster};   % [nPsi mice x num_bins]
    veh_mat = veh_clusters{cluster};   % [nVeh mice x num_bins]

    % Skip if cluster has no data in either group
    if isempty(psi_mat) && isempty(veh_mat)
        title(sprintf('Cluster %d (no data)', cluster));
        axis off;
        continue;
    end

    % Remove mice with no valid data and compute vehicle mean/SEM
    if ~isempty(veh_mat)
        valid_veh = ~all(isnan(veh_mat), 2);
        veh_mat = veh_mat(valid_veh, :);
        mean_veh = mean(veh_mat, 1, 'omitnan');
        sem_veh = std(veh_mat, 0, 1, 'omitnan') ./ sqrt(size(veh_mat, 1));
    else
        mean_veh = nan(1, num_bins);
        sem_veh = nan(1, num_bins);
    end

    % Remove mice with no valid data and compute psilocybin mean/SEM
    if ~isempty(psi_mat)
        valid_psi = ~all(isnan(psi_mat), 2);
        psi_mat = psi_mat(valid_psi, :);
        mean_psi = mean(psi_mat, 1, 'omitnan');
        sem_psi = std(psi_mat, 0, 1, 'omitnan') ./ sqrt(size(psi_mat, 1));
    else
        mean_psi = nan(1, num_bins);
        sem_psi = nan(1, num_bins);
    end

    % Plot individual vehicle mouse PSTHs as faint lines
    if ~isempty(veh_mat)
        for i = 1:size(veh_mat, 1)
            plot(t, veh_mat(i, :), 'Color', [veh_col 0.3], 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
    end

    % Plot individual psilocybin mouse PSTHs as faint lines
    if ~isempty(psi_mat)
        for i = 1:size(psi_mat, 1)
            plot(t, psi_mat(i, :), 'Color', [psi_col 0.3], 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
    end

    % Vehicle shaded SEM and mean line
    veh_upper = mean_veh + sem_veh;
    veh_lower = mean_veh - sem_veh;

    fill([t fliplr(t)], [veh_upper fliplr(veh_lower)], veh_col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(t, mean_veh, 'Color', veh_col, 'LineWidth', 1.5);

    % Psilocybin shaded SEM and mean line
    psi_upper = mean_psi + sem_psi;
    psi_lower = mean_psi - sem_psi;

    fill([t fliplr(t)], [psi_upper fliplr(psi_lower)], psi_col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(t, mean_psi, 'Color', psi_col, 'LineWidth', 1.5);

    title(sprintf('Cluster %d', cluster));
    xlabel('Time (s)');
    ylabel('Mean Spike Count/Bin');
    xlim([t(1) t(end)]);

    ax = gca;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    ax.Box = 'off';

    % Stimulus onset/offset markers
    xline(10, 'k--', 'LineWidth', 1.5);
    xline(13.3, 'k--', 'LineWidth', 1.5);

    % Figure-level title depending on session and normalisation
    if do_normalise && strcmp(session, 'Extinction')
        sgtitle('Extinction Mean PSTHs per Cluster (Baseline Subtracted)');
    elseif ~do_normalise && strcmp(session, 'Extinction')
        sgtitle('Extinction Mean PSTHs per Cluster');
    elseif do_normalise && strcmp(session, 'Renewal')
        sgtitle('Renewal Mean PSTHs per Cluster (Baseline Subtracted)');
    elseif ~do_normalise && strcmp(session, 'Renewal')
        sgtitle('Renewal Mean PSTHs per Cluster');
    end

    % Optional fixed y-axis
    % ylim([0 15])

    if cluster == 1
        legend({'Veh SEM', 'Vehicle', 'Psi SEM', 'Psilocybin'}, ...
               'Location', 'best', 'Box', 'off');
    end
end

%% Statistical tests: compare baseline firing between treatment groups

% This tests whether baseline firing differs between vehicle and psilocybin
% groups for each cluster. This is useful because differences in the PSTH
% after stimulus presentation may partly reflect baseline firing differences.

baseline_bins = 1:19;

p_vals = nan(num_clusters, 1);
baseline_rates_all = nan(num_mice, num_clusters);  % row 1 empty if Mouse 1 is unused

for cluster = 1:num_clusters

    % Mean baseline firing per mouse for this cluster
    baseline_rates = mean(avg_spiking_per_bin_clusters{cluster}(:, baseline_bins), 2, 'omitnan');
    baseline_rates_all(:, cluster) = baseline_rates;

    % Extract treatment groups
    x_ctrl = baseline_rates(received_veh);
    x_psi = baseline_rates(received_psi);

    % Non-parametric comparison between treatment groups
    p_vals(cluster) = ranksum(x_ctrl, x_psi);
end

disp(p_vals);

%% Compute block-wise PSTHs for each cluster

% This section examines how PSTHs evolve over trial blocks.
% With 20 trials and 5 trials per block, there are 4 blocks.

num_blocks = 4;
avg_spiking_per_bin_byBlock_clusters = cell(num_clusters, 1);

for cluster = 1:num_clusters

    % A has dimensions:
    %   mouse x time bin x trial block
    A = nan(num_mice, num_bins, num_blocks);

    for mouse = mice_to_analyse
        T = sorted_clusters{mouse, cluster};

        % Skip empty data or data without TrialBlock labels
        if isempty(T) || height(T) == 0 || ~ismember('TrialBlock', T.Properties.VariableNames)
            continue
        end

        for b = 1:num_blocks
            Tb = T(T.TrialBlock == b, :);
            if isempty(Tb)
                continue
            end

            for bin = 1:num_bins
                A(mouse, bin, b) = mean(Tb.SpikeCount(Tb.TimeBin == bin), 'omitnan');
            end
        end
    end

    avg_spiking_per_bin_byBlock_clusters{cluster} = A;
end

%% Plot block-wise PSTHs for selected clusters

% Choose which clusters to display
clusters_to_plot = [1 2];

num_blocks = 4;
veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

figure('Color', 'w');

% Create two large side-by-side panels, one for each cluster.
p1 = uipanel('Parent', gcf, 'Units', 'normalized', ...
    'Position', [0.02 0.08 0.47 0.87], 'BorderType', 'none');

p2 = uipanel('Parent', gcf, 'Units', 'normalized', ...
    'Position', [0.51 0.08 0.47 0.87], 'BorderType', 'none');

panels = {p1, p2};

for s = 1:2
    cluster = clusters_to_plot(s);

    % A is mouse x bin x block for this cluster
    A = avg_spiking_per_bin_byBlock_clusters{cluster};

    % 2x2 layout inside each cluster panel: one tile per trial block
    tl = tiledlayout(panels{s}, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Cluster %d', cluster), 'FontWeight', 'bold');

    for b = 1:num_blocks
        ax = nexttile(tl); hold(ax, 'on');

        % Extract block-specific PSTHs for each group
        veh_mat = squeeze(A(received_veh, :, b));
        psi_mat = squeeze(A(received_psi, :, b));

        % Mean and SEM across mice
        mean_veh = mean(veh_mat, 1, 'omitnan');
        mean_psi = mean(psi_mat, 1, 'omitnan');

        sem_veh = std(veh_mat, 0, 1, 'omitnan') ./ sqrt(size(veh_mat, 1));
        sem_psi = std(psi_mat, 0, 1, 'omitnan') ./ sqrt(size(psi_mat, 1));

        % Vehicle shaded SEM and mean line
        fill(ax, [t fliplr(t)], [mean_veh + sem_veh fliplr(mean_veh - sem_veh)], ...
            veh_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(ax, t, mean_veh, 'Color', veh_col, 'LineWidth', 2);

        % Psilocybin shaded SEM and mean line
        fill(ax, [t fliplr(t)], [mean_psi + sem_psi fliplr(mean_psi - sem_psi)], ...
            psi_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(ax, t, mean_psi, 'Color', psi_col, 'LineWidth', 2);

        % Stimulus onset/offset markers
        xline(ax, 10, 'k--', 'LineWidth', 1.2);
        xline(ax, 13.3, 'k--', 'LineWidth', 1.2);

        % Trial range represented by this block
        block_start = (b - 1) * 5 + 1;
        block_end = b * 5;
        title(ax, sprintf('Trials %d–%d', block_start, block_end));

        xlim(ax, [0 33]);
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Mean spike count/bin');

        ax.TickDir = 'out';
        ax.Box = 'off';
        ax.LineWidth = 1.2;

        % Optional fixed y-axis for raw, non-normalised data
        if ~do_normalise
            ylim([1 8]);
        end

        if b == 1
            legend(ax, {'Veh SEM', 'Vehicle', 'Psi SEM', 'Psilocybin'}, ...
                'Location', 'best', 'Box', 'off');
        end
    end
end
%%
%% Plot PSTHs by brain region: SC pooled versus PAG

sc_regions = ["SCdg","SCdw","SCig","SCiw","SCop","SCsg","SCzo"];
region_names = ["SC", "PAG"];
num_regions = numel(region_names);

region_psth_data = cell(num_mice, num_regions);
n_cells_per_mouse_region = nan(num_mice, num_regions);

for mouse = mice_to_analyse

    current_data = hist_and_spiking_matrices{mouse};

    if isempty(current_data) || height(current_data) == 0
        continue
    end

    % Keep only selected trial type
    current_data = current_data(current_data.TrialIdentifier == trial_type, :);

    if isempty(current_data) || height(current_data) == 0
        continue
    end

    % Reconstruct TimeBin
    if mod(height(current_data), num_bins) ~= 0
        warning(['Mouse ' num2str(mouse) ': row count not divisible by num_bins'])
        continue
    end

    current_data.TimeBin = repmat((1:num_bins)', height(current_data)/num_bins, 1);

    % Define region membership
    sc_idx = ismember(current_data.Brain_Region, sc_regions);
    pag_idx = current_data.Brain_Region == "PAG";

    region_idx = {sc_idx, pag_idx};

    for region = 1:num_regions

        region_data = current_data(region_idx{region}, :);

        if isempty(region_data) || height(region_data) == 0
            region_psth_data{mouse, region} = region_data;
            n_cells_per_mouse_region(mouse, region) = 0;
            continue
        end

        n_cells_per_mouse_region(mouse, region) = numel(unique(region_data.CellID));

        % Optional baseline subtraction per cell
        if do_normalise

            cellIDs = unique(region_data.CellID);

            for c = 1:numel(cellIDs)

                cid = cellIDs(c);
                cell_idx = region_data.CellID == cid;

                pre_idx = cell_idx & ...
                    region_data.TimeBin >= pre_stim_bins(1) & ...
                    region_data.TimeBin <= pre_stim_bins(end);

                baseline = mean(region_data.SpikeCount(pre_idx), 'omitnan');

                region_data.SpikeCount(cell_idx) = ...
                    region_data.SpikeCount(cell_idx) - baseline;

            end
        end

        region_psth_data{mouse, region} = region_data;

    end
end

%% Average spike count per bin for each region and mouse

avg_spiking_per_bin_regions = cell(num_regions, 1);

for region = 1:num_regions

    A = nan(num_mice, num_bins);

    for mouse = mice_to_analyse

        region_data = region_psth_data{mouse, region};

        if isempty(region_data) || height(region_data) == 0
            continue
        end

        for bin = 1:num_bins
            A(mouse, bin) = mean(region_data.SpikeCount(region_data.TimeBin == bin), 'omitnan');
        end
    end

    avg_spiking_per_bin_regions{region} = A;

end

%%
%% Plot region PSTHs: Vehicle vs Psilocybin

veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

figure('Color', 'w');

for region = 1:num_regions

    subplot(1, num_regions, region); hold on;

    A = avg_spiking_per_bin_regions{region};

    veh_mat = A(received_veh, :);
    psi_mat = A(received_psi, :);

    % Remove mice with no data for this region
    veh_mat = veh_mat(~all(isnan(veh_mat), 2), :);
    psi_mat = psi_mat(~all(isnan(psi_mat), 2), :);

    mean_veh = mean(veh_mat, 1, 'omitnan');
    mean_psi = mean(psi_mat, 1, 'omitnan');

    sem_veh = std(veh_mat, 0, 1, 'omitnan') ./ sqrt(size(veh_mat, 1));
    sem_psi = std(psi_mat, 0, 1, 'omitnan') ./ sqrt(size(psi_mat, 1));

    % Individual mouse lines
    for i = 1:size(veh_mat, 1)
        plot(t, veh_mat(i, :), 'Color', [veh_col 0.25], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end

    for i = 1:size(psi_mat, 1)
        plot(t, psi_mat(i, :), 'Color', [psi_col 0.25], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % Vehicle mean ± SEM
    fill([t fliplr(t)], ...
        [mean_veh + sem_veh fliplr(mean_veh - sem_veh)], ...
        veh_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    plot(t, mean_veh, 'Color', veh_col, 'LineWidth', 2);

    % Psilocybin mean ± SEM
    fill([t fliplr(t)], ...
        [mean_psi + sem_psi fliplr(mean_psi - sem_psi)], ...
        psi_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    plot(t, mean_psi, 'Color', psi_col, 'LineWidth', 2);

    % Count n cells and n mice
    veh_cell_n = sum(n_cells_per_mouse_region(received_veh, region), 'omitnan');
    psi_cell_n = sum(n_cells_per_mouse_region(received_psi, region), 'omitnan');

    veh_mouse_n = sum(n_cells_per_mouse_region(received_veh, region) > 0);
    psi_mouse_n = sum(n_cells_per_mouse_region(received_psi, region) > 0);

    title({
        char(region_names(region))
        sprintf('Veh: n=%d cells/%d mice | Psi: n=%d cells/%d mice', ...
        veh_cell_n, veh_mouse_n, psi_cell_n, psi_mouse_n)
        });

    xline(10, 'k--', 'LineWidth', 1.5);
    xline(13.3, 'k--', 'LineWidth', 1.5);

    xlabel('Time (s)');

    if do_normalise
        ylabel('Mean spike count/bin baseline-subtracted');
    else
        ylabel('Mean spike count/bin');
    end

    xlim([0 33]);

    ax = gca;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.LineWidth = 1.2;

    if region == 1
        legend({'Veh SEM', 'Vehicle', 'Psi SEM', 'Psilocybin'}, ...
            'Location', 'best', 'Box', 'off');
    end

end

if do_normalise
    sgtitle([session ' PSTHs by Brain Region Baseline Subtracted']);
else
    sgtitle([session ' PSTHs by Brain Region']);
end

%%
%% Plot PSTHs by SC subregion plus PAG

% Define each anatomical region separately
region_names = ["SCdg", "SCdw", "SCig", "SCiw", "SCop", "SCsg", "SCzo", "PAG"];
num_regions = numel(region_names);

region_psth_data = cell(num_mice, num_regions);
n_cells_per_mouse_region = nan(num_mice, num_regions);

for mouse = mice_to_analyse

    current_data = hist_and_spiking_matrices{mouse};

    if isempty(current_data) || ~istable(current_data) || height(current_data) == 0
        continue
    end

    % Keep only selected trial type
    current_data = current_data(current_data.TrialIdentifier == trial_type, :);

    if isempty(current_data) || height(current_data) == 0
        continue
    end

    % Reconstruct TimeBin
    if mod(height(current_data), num_bins) ~= 0
        warning(['Mouse ' num2str(mouse) ': row count not divisible by num_bins'])
        continue
    end

    current_data.TimeBin = repmat((1:num_bins)', height(current_data)/num_bins, 1);

    for region = 1:num_regions

        region_label = region_names(region);

        % Select rows for this specific region
        region_idx = current_data.Brain_Region == region_label;
        region_data = current_data(region_idx, :);

        if isempty(region_data) || height(region_data) == 0
            region_psth_data{mouse, region} = table();
            n_cells_per_mouse_region(mouse, region) = 0;
            continue
        end

        % Number of unique cells from this mouse in this region
        n_cells_per_mouse_region(mouse, region) = numel(unique(region_data.CellID));

        % Optional baseline subtraction per cell
        if do_normalise

            cellIDs = unique(region_data.CellID);

            for c = 1:numel(cellIDs)

                cid = cellIDs(c);
                cell_idx = region_data.CellID == cid;

                pre_idx = cell_idx & ...
                    region_data.TimeBin >= pre_stim_bins(1) & ...
                    region_data.TimeBin <= pre_stim_bins(end);

                baseline = mean(region_data.SpikeCount(pre_idx), 'omitnan');

                region_data.SpikeCount(cell_idx) = ...
                    region_data.SpikeCount(cell_idx) - baseline;

            end
        end

        region_psth_data{mouse, region} = region_data;

    end
end

%% Average spike count per bin for each SC/PAG region and mouse

avg_spiking_per_bin_regions = cell(num_regions, 1);

for region = 1:num_regions

    A = nan(num_mice, num_bins);

    for mouse = mice_to_analyse

        region_data = region_psth_data{mouse, region};

        if isempty(region_data) || ~istable(region_data) || height(region_data) == 0
            continue
        end

        for bin = 1:num_bins
            bin_idx = region_data.TimeBin == bin;
            A(mouse, bin) = mean(region_data.SpikeCount(bin_idx), 'omitnan');
        end
    end

    avg_spiking_per_bin_regions{region} = A;

end

%% Plot SC subregion/PAG PSTHs: Vehicle vs Psilocybin

veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

figure('Color', 'w');

num_rows = 2;
num_cols = 4;

for region = 1:num_regions

    subplot(num_rows, num_cols, region); hold on;

    A = avg_spiking_per_bin_regions{region};

    veh_mat = A(received_veh, :);
    psi_mat = A(received_psi, :);

    % Remove mice with no data for this region
    veh_mat = veh_mat(~all(isnan(veh_mat), 2), :);
    psi_mat = psi_mat(~all(isnan(psi_mat), 2), :);

    % Skip completely empty regions
    if isempty(veh_mat) && isempty(psi_mat)
        title(char(region_names(region)));
        axis off
        continue
    end

    mean_veh = mean(veh_mat, 1, 'omitnan');
    mean_psi = mean(psi_mat, 1, 'omitnan');

    sem_veh = std(veh_mat, 0, 1, 'omitnan') ./ sqrt(size(veh_mat, 1));
    sem_psi = std(psi_mat, 0, 1, 'omitnan') ./ sqrt(size(psi_mat, 1));

    % Individual mouse lines
    for i = 1:size(veh_mat, 1)
        plot(t, veh_mat(i, :), 'Color', [veh_col 0.25], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end

    for i = 1:size(psi_mat, 1)
        plot(t, psi_mat(i, :), 'Color', [psi_col 0.25], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % Vehicle mean ± SEM
    if ~isempty(veh_mat)
        fill([t fliplr(t)], ...
            [mean_veh + sem_veh fliplr(mean_veh - sem_veh)], ...
            veh_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        plot(t, mean_veh, 'Color', veh_col, 'LineWidth', 2);
    end

    % Psilocybin mean ± SEM
    if ~isempty(psi_mat)
        fill([t fliplr(t)], ...
            [mean_psi + sem_psi fliplr(mean_psi - sem_psi)], ...
            psi_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        plot(t, mean_psi, 'Color', psi_col, 'LineWidth', 2);
    end

    % Cell and mouse counts
    veh_cell_n = sum(n_cells_per_mouse_region(received_veh, region), 'omitnan');
    psi_cell_n = sum(n_cells_per_mouse_region(received_psi, region), 'omitnan');

    veh_mouse_n = sum(n_cells_per_mouse_region(received_veh, region) > 0);
    psi_mouse_n = sum(n_cells_per_mouse_region(received_psi, region) > 0);

    title({
        char(region_names(region))
        sprintf('Veh: %d cells/%d mice | Psi: %d cells/%d mice', ...
        veh_cell_n, veh_mouse_n, psi_cell_n, psi_mouse_n)
        });

    xline(10, 'k--', 'LineWidth', 1.2);
    xline(13.3, 'k--', 'LineWidth', 1.2);

    xlabel('Time (s)');

    if do_normalise
        ylabel('Baseline-subtracted spike count/bin');
    else
        ylabel('Mean spike count/bin');
    end

    xlim([0 33]);

    ax = gca;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.LineWidth = 1.2;

    if region == 1
        legend({'Veh SEM', 'Vehicle', 'Psi SEM', 'Psilocybin'}, ...
            'Location', 'best', 'Box', 'off');
    end

end

if do_normalise
    sgtitle([session ' PSTHs by SC Subregion/PAG Baseline Subtracted']);
else
    sgtitle([session ' PSTHs by SC Subregion/PAG']);
end

%%
%% Plot PSTHs by broad SC functional groups: optic vs behaviour SC

% Functional/anatomical grouping
optic_SC_regions     = ["SCzo", "SCsg", "SCop"];
behaviour_SC_regions = ["SCig", "SCdg"];

group_names = ["Optic SC", "Behaviour SC"];
num_groups = numel(group_names);

group_region_lists = {
    optic_SC_regions
    behaviour_SC_regions
};

group_psth_data = cell(num_mice, num_groups);
n_cells_per_mouse_group = nan(num_mice, num_groups);

for mouse = mice_to_analyse

    current_data = hist_and_spiking_matrices{mouse};

    if isempty(current_data) || ~istable(current_data) || height(current_data) == 0
        continue
    end

    % Keep only selected trial type
    current_data = current_data(current_data.TrialIdentifier == trial_type, :);

    if isempty(current_data) || height(current_data) == 0
        continue
    end

    % Reconstruct TimeBin
    if mod(height(current_data), num_bins) ~= 0
        warning(['Mouse ' num2str(mouse) ': row count not divisible by num_bins'])
        continue
    end

    current_data.TimeBin = repmat((1:num_bins)', height(current_data)/num_bins, 1);

    for group = 1:num_groups

        current_region_list = group_region_lists{group};

        group_idx = ismember(current_data.Brain_Region, current_region_list);
        group_data = current_data(group_idx, :);

        if isempty(group_data) || height(group_data) == 0
            group_psth_data{mouse, group} = table();
            n_cells_per_mouse_group(mouse, group) = 0;
            continue
        end

        % Count unique cells in this mouse/group
        n_cells_per_mouse_group(mouse, group) = numel(unique(group_data.CellID));

        % Optional baseline subtraction per cell
        if do_normalise

            cellIDs = unique(group_data.CellID);

            for c = 1:numel(cellIDs)

                cid = cellIDs(c);
                cell_idx = group_data.CellID == cid;

                pre_idx = cell_idx & ...
                    group_data.TimeBin >= pre_stim_bins(1) & ...
                    group_data.TimeBin <= pre_stim_bins(end);

                baseline = mean(group_data.SpikeCount(pre_idx), 'omitnan');

                group_data.SpikeCount(cell_idx) = ...
                    group_data.SpikeCount(cell_idx) - baseline;

            end
        end

        group_psth_data{mouse, group} = group_data;

    end
end

%% Average spike count per bin for each SC group and mouse

avg_spiking_per_bin_groups = cell(num_groups, 1);

for group = 1:num_groups

    A = nan(num_mice, num_bins);

    for mouse = mice_to_analyse

        group_data = group_psth_data{mouse, group};

        if isempty(group_data) || ~istable(group_data) || height(group_data) == 0
            continue
        end

        for bin = 1:num_bins
            bin_idx = group_data.TimeBin == bin;
            A(mouse, bin) = mean(group_data.SpikeCount(bin_idx), 'omitnan');
        end
    end

    avg_spiking_per_bin_groups{group} = A;

end

%% Plot optic SC vs behaviour SC PSTHs

veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

figure('Color', 'w');

for group = 1:num_groups

    subplot(1, num_groups, group); hold on;

    A = avg_spiking_per_bin_groups{group};

    veh_mat = A(received_veh, :);
    psi_mat = A(received_psi, :);

    % Remove mice with no data for this group
    veh_mat = veh_mat(~all(isnan(veh_mat), 2), :);
    psi_mat = psi_mat(~all(isnan(psi_mat), 2), :);

    % Skip empty group
    if isempty(veh_mat) && isempty(psi_mat)
        title(char(group_names(group)));
        axis off
        continue
    end

    % Mean and SEM across mice
    mean_veh = mean(veh_mat, 1, 'omitnan');
    mean_psi = mean(psi_mat, 1, 'omitnan');

    sem_veh = std(veh_mat, 0, 1, 'omitnan') ./ sqrt(size(veh_mat, 1));
    sem_psi = std(psi_mat, 0, 1, 'omitnan') ./ sqrt(size(psi_mat, 1));

    % Individual mouse lines
    for i = 1:size(veh_mat, 1)
        plot(t, veh_mat(i, :), 'Color', [veh_col 0.25], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end

    for i = 1:size(psi_mat, 1)
        plot(t, psi_mat(i, :), 'Color', [psi_col 0.25], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % Vehicle mean ± SEM
    if ~isempty(veh_mat)
        fill([t fliplr(t)], ...
            [mean_veh + sem_veh fliplr(mean_veh - sem_veh)], ...
            veh_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        plot(t, mean_veh, 'Color', veh_col, 'LineWidth', 2);
    end

    % Psilocybin mean ± SEM
    if ~isempty(psi_mat)
        fill([t fliplr(t)], ...
            [mean_psi + sem_psi fliplr(mean_psi - sem_psi)], ...
            psi_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        plot(t, mean_psi, 'Color', psi_col, 'LineWidth', 2);
    end

    % Count cells and mice
    veh_cell_n = sum(n_cells_per_mouse_group(received_veh, group), 'omitnan');
    psi_cell_n = sum(n_cells_per_mouse_group(received_psi, group), 'omitnan');

    veh_mouse_n = sum(n_cells_per_mouse_group(received_veh, group) > 0);
    psi_mouse_n = sum(n_cells_per_mouse_group(received_psi, group) > 0);

    title({
        char(group_names(group))
        sprintf('Veh: %d cells/%d mice | Psi: %d cells/%d mice', ...
        veh_cell_n, veh_mouse_n, psi_cell_n, psi_mouse_n)
        });

    xline(10, 'k--', 'LineWidth', 1.5);
    xline(13.3, 'k--', 'LineWidth', 1.5);

    xlabel('Time (s)');

    if do_normalise
        ylabel('Baseline-subtracted spike count/bin');
    else
        ylabel('Mean spike count/bin');
    end

    xlim([0 33]);

    ax = gca;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.LineWidth = 1.2;

    if group == 1
        legend({'Veh SEM', 'Vehicle', 'Psi SEM', 'Psilocybin'}, ...
            'Location', 'best', 'Box', 'off');
    end

end

if do_normalise
    sgtitle([session ' PSTHs: Optic SC vs Behaviour SC Baseline Subtracted']);
else
    sgtitle([session ' PSTHs: Optic SC vs Behaviour SC']);
end