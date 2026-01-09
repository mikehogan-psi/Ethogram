% Calculates y position of kilosort clusters along the length of the probe
% starting from the region closest to the tip
% 
% OUTPUT: A table containing 3 rows: 
% 1. Relative cell depth- depth of each cluster relative to channels
% closest to the tip 
% 2. Cell depth- depth of each cluster starting from the first recording
% site
% 3. Cell ID: Cell ID of associated cluster
%
% To calculates cluster depth, this script finds the channel with
% largest peak-to-peak value for each template, finds the depth of that
% template and then assigns each cluster the depth associated with its
% modal template (i.e. template assigned to the majority of spikes assigned
% to a cluster)


%%  Setup (input required)

master_directory = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)";

mice_to_analyse = [2 3 4 5 6 7 8 9];

session = 'Renewal';
%% Get filepaths for neural data

mouse_files = dir(fullfile(master_directory, 'Mouse*'));

num_mice = numel(mice_to_analyse);

% Extract kilosort filepaths
neural_data_folders = cell(num_mice, 1);

for mouse = mice_to_analyse
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    neural_data_folders{mouse} = fullfile(mouse_path, mouse_name, ...
        session, 'Neural Data', 'Concatenated Data', 'kilosort4\');           
end



%% Get save paths for cell positions

% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

data_folders = cell(length(mice_to_analyse), 1);

% Get behaviour folders for specified session
for mouse = mice_to_analyse
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    data_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Neural Data', 'Concatenated Data');
end
%%

for mouse = mice_to_analyse

    kilosort_path = neural_data_folders{mouse};

    chan_map      = readNPY(fullfile(kilosort_path,'channel_map.npy'));        % length = nCh
    chan_pos      = readNPY(fullfile(kilosort_path,'channel_positions.npy'));  % [nCh x 2], columns: x,y (µm)
    spike_clusters     = readNPY(fullfile(kilosort_path,'spike_clusters.npy'));     % [nSpikes x 1]
    spike_templates    = readNPY(fullfile(kilosort_path,'spike_templates.npy'));    % [nSpikes x 1]
    templates    = readNPY(fullfile(kilosort_path,'templates.npy'));          % [nTemplates x nTime x nChan]
   
 
    % Calculate peak to peak amplitude across time for each template (putative
    % cell) for each channel
    ptp = squeeze(max(templates, [], 2) - min(templates, [], 2));
    
    % Get indices of channels with peak template amplitudes
    [max_ptp, best_channels_idx] = max(ptp, [], 2);
    
    % Find position of peak amplitude templates on probe
    best_positions_templates_x = chan_pos(best_channels_idx, 1);
    best_positions_templates_y = chan_pos(best_channels_idx, 2);
    
    % Channel ID for each template
    best_channel_for_template = chan_map(best_channels_idx);   

    % Get unique cell IDs
    cluster_values = double(unique(spike_clusters));
    n_clusters = numel(cluster_values);
    depth_clusters = zeros(n_clusters, 1);
    lat_clusters = zeros(n_clusters, 1);
    peak_channel_ids = zeros(n_clusters, 1);
    
    % Find most common template associated with cluster and assign cluster
    % x, y and channel associated with that template
    for i = 1:numel(cluster_values)
        cell_id = cluster_values(i);
        % Find all template ids allocated to this cluster
        template_idxs = spike_templates(spike_clusters == cell_id);
        % Find most common template id associated with cluster
        common_template = mode(template_idxs) + 1; % +1 because template ids start at 0
        % Assign cluster depth of most common template
        current_cell_depth = best_positions_templates_y(common_template);
        depth_clusters(i) = current_cell_depth;
        current_cell_lat = best_positions_templates_x(common_template);
        lat_clusters(i) = current_cell_lat;
        % Assign cluster channel associated with template
        peak_channel_ids(i) = best_channel_for_template(common_template);
    end
    
    
    % Subtract minimum value to get positions of cells relative to position
    % of first recording site
    y0 = min(chan_pos(:, 2));
    x0 = min(chan_pos(:, 1));
    
    depth_clusters = depth_clusters-y0;
    lat_clusters = lat_clusters-x0;

    % Generate table
    depth_table = table(lat_clusters, depth_clusters,  cluster_values,...
        peak_channel_ids, 'VariableNames', {'x', 'y', 'CellID', 'Channel'});
    
 
    mouse_name = mouse_files(mouse).name;
    mouse_name = lower(strrep(mouse_name, ' ', ''));
    save_name = [mouse_name, '_', lower(session), '_', 'cell_depths.mat'];
    save_path = fullfile(data_folders{mouse}, save_name);
    save(save_path, "depth_table");
    disp([save_name ' saved to ' data_folders{mouse} '!'])

end
