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

mice_to_analyse = [2 3 4 5 6 7];

session = 'Extinction';
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
    [max_ptp, best_channels_templates] = max(ptp, [], 2);
    
    % Find position of peak amplitude templates on probe
    best_positions_templates = chan_pos(best_channels_templates, 2);
    
    % Get unique cell IDs
    cluster_values = double(unique(spike_clusters));
    
    depth_clusters = zeros(size(cluster_values, 1), 1);
    
    % Find most common template associated with cluster and assign cluster
    % depth associated with that template
    for i = 1:numel(cluster_values)
        cell_id = cluster_values(i);
        % Find all template ids allocated to thi cluster
        template_idxs = spike_templates(spike_clusters == cell_id);
        % Find most common template id associated with cluster
        common_template = mode(template_idxs) + 1; % +1 because template ids start at 0
        % Assign cluster depth of most common template
        current_cell_depth = best_positions_templates(common_template);
        depth_clusters(i) = current_cell_depth;
    end
    
    min_depth = min(depth_clusters);
    
    non_relative_depth = depth_clusters - min_depth;
    
    % Generate table
    depth_table = table(depth_clusters, non_relative_depth, cluster_values,...
        'VariableNames', {'Relative Depth (µm)', 'Depth (µm)', 'Cell ID'});
    
    mouse_name = mouse_files(mouse).name;
    mouse_name = lower(strrep(mouse_name, ' ', ''));
    save_name = [mouse_name, '_', lower(session), '_', 'cell_depths.mat'];
    save_path = fullfile(data_folders{mouse}, save_name);
    save(save_path, "depth_table");
    disp([save_name ' saved to ' data_folders{mouse} '!'])

end