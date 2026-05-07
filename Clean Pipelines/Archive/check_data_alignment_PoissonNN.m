% pick a cell to visualize
num_bins = 66;
ycol = y(:);
ids = unique(cell_ids(:))';   % e.g., [0 1 2 ... N]

for neuron_id = ids      % <- loops over 0,1,2,... exactly as stored
    idx = find(cell_ids == neuron_id);
    if isempty(idx), continue; end

    n_trials = numel(idx)/num_bins;
    if mod(numel(idx), num_bins) ~= 0, warning('skip %d', neuron_id); continue; end

    Y = reshape(ycol(idx), num_bins, n_trials)';   % [n_trials x 66]
    imagesc(Y); set(gca,'YDir','normal'); colormap hot; colorbar
    title(sprintf('Neuron %d', neuron_id));

        % annotate stimulus window and loom/flash boundary
    fps = 15; bin_size = 0.5;
    stim_start_bin = ceil((152/fps) / bin_size);
    stim_end_bin   = ceil((202/fps) / bin_size);
    
    hold on
    xline([stim_start_bin stim_end_bin] + 0.5, 'w--', 'LineWidth', 1);   % stim window
    yline(20.5, 'w:');  % looms (1..20) vs flashes (21..40), per your stacking
    hold off
    % pause(0.1)
    ginput(); close(gcf);
end
%%
% Ensure X is features x samples
if size(X,1) ~= 10 && size(X,2) == 10
    X = X.';
end

num_bins = 66;  % 0.5 s bins
fps = 15; bin_size = 0.5;
stim_start_bin = ceil((152/fps) / bin_size);
stim_end_bin   = ceil((202/fps) / bin_size);
feat_names = {'Grooming','Rearing','Darting','Freezing','Velocity'};

ids = unique(cell_ids(:))';

for neuron_id = ids
    idx = find(cell_ids == neuron_id);
    if isempty(idx), continue; end
    if mod(numel(idx), num_bins) ~= 0
        warning('Neuron %d: |idx|=%d not multiple of %d bins, skipping', neuron_id, numel(idx), num_bins);
        continue
    end
    n_trials = numel(idx) / num_bins;

    figure('Name', sprintf('Neuron %d – features 1–5', neuron_id));
    t = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

    for f = 1:5
        % Select row f (feature) and the columns for this neuron, then columnize
        col = X(f, idx).';
        Yf  = reshape(col, num_bins, n_trials).';

        nexttile;
        imagesc(Yf);
        set(gca,'YDir','normal');
        colormap hot; colorbar;
        title(feat_names{f});
        xlabel('Time bin (0.5 s)'); ylabel('Trial');

        hold on
        xline([stim_start_bin stim_end_bin] + 0.5, 'w--', 'LineWidth', 1);
        if n_trials == 40, yline(20.5, 'w:'); end
        hold off
    end

    title(t, sprintf('Neuron %d – X features 1–5', neuron_id));
    ginput(); close(gcf);
end
