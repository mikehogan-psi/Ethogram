% EXTRACT_FIRING_HABITUATION
%
% Extracts and plots habituation firing-rate time courses from Kilosort output.
% For each mouse, uses the first habituation trigger as t0, bins spikes into
% 0.5 s bins over the whole hab session (300s), averages across cells to get a per-mouse
% mean firing-rate trace, then aggregates mice into VEH vs PSI groups and
% plots group mean ± SEM (plus an optional moving-average smoothed trace).
%% Directory Setup

master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';
 
session = 'Renewal';

mice_to_analyse = [2 3 4 5 6 7 8 9];
max_mouse_id    = max(mice_to_analyse); 

%% Load neural data

% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

num_mice = numel(mice_to_analyse);

n_trials  = 150; % total number of trials
fs = 30000;  % Ephys sampling rate (30kHz)

% Extract trigger file filepaths
evt_folders = cell(max_mouse_id, 1);

for mouse = mice_to_analyse
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        evt_folders{mouse} = fullfile(mouse_path, mouse_name, ...
            session, 'Neural Data', 'Triggers'); 
end

% Load trigger files
evt_data = cell(num_mice, 1);
for mouse = mice_to_analyse
    current_evt_folder = evt_folders{mouse};  
    current_evt_file_1 = dir(fullfile(current_evt_folder, '*_checkerboard_triggers.mat'));
    tmp1 = load(fullfile(current_evt_file_1.folder, current_evt_file_1.name));
    evt_data{mouse} = tmp1.evt(2:end);
end


% Extract trigger timings, convert to seconds, 
% and then add extra trigger at the end as edge for input to histc

for mouse = mice_to_analyse
    current_mouse_evt = evt_data{mouse};
    current_mouse_evt = double(current_mouse_evt);
    current_mouse_evt = current_mouse_evt/fs;
    bin_duration = median(diff(current_mouse_evt));  % (nTrials x 1) per-trial bin width
    % Append one extra edge per trial equal to last bin + avg bin duration
    % for that trial
    evt_data{mouse} = [current_mouse_evt; current_mouse_evt(end) + bin_duration];  

end

% Extract kilosort filepaths
neural_data_folders = cell(num_mice, 1);

for mouse = mice_to_analyse
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    neural_data_folders{mouse} = fullfile(mouse_path, mouse_name, ...
        session, 'Neural Data', 'Concatenated Data', 'kilosort4\');           
    
end

% Load neural data from kilosort paths
spike_matrix = cell(num_mice, 1);
cluster_value_matrix = cell(num_mice, 1);

for mouse = mice_to_analyse
    % Loads n_spikes x 1 vector containing cell IDs associated with each
    % spike recorded
    clusters = readNPY(fullfile(neural_data_folders{mouse}, 'spike_clusters.npy'));
    % Loads n_spikes x 1 vector containing spike times (in raw 30kHz sample
    % format) of each spike recorded
    spike_times = readNPY(fullfile(neural_data_folders{mouse}, 'spike_times.npy'));
    % Convert spike times to seconds
    spike_times = double(spike_times) / fs;  
    % Find unique cell IDs
    cluster_values = unique(clusters);  
    Ncell = numel(cluster_values);

    % Store both the spike times and the cluster labels
    spike_matrix{mouse} = spike_times;
    cluster_id_matrix{mouse} = clusters;
    cluster_value_matrix{mouse} = double(cluster_values);
end


%% Checkerboard spiking in 0.5 s bins (per cell)

trial_length = 300;   % 2 s

bin_size = 2;                              % seconds
num_bins = trial_length / bin_size;          % 600 bins

% Store one matrix per mouse: [nCells x num_bins]
checkerboard_binned_spiking = cell(num_mice, 1);

for mouse = mice_to_analyse

    % --- habituation start time ---
    % evt_data{mouse} should be in SECONDS already
    % (if not, divide by fs before this step)
    t0 = evt_data{mouse}(1);                 % START of habituation in seconds

    % define bin edges from t0 to t0 + 300 s
    edges = t0 : bin_size : (t0 + trial_length);   % 0.5 s bins → 601 edges

    % --- spikes for this mouse ---
    spike_times   = spike_matrix{mouse};        % seconds
    clusters      = cluster_id_matrix{mouse};
    cluster_vals  = cluster_value_matrix{mouse};
    nCells        = numel(cluster_vals);

    counts = zeros(nCells, num_bins);  % [cells x bins]

    for c = 1:nCells
        st = spike_times(clusters == cluster_vals(c));   % spikes for this cell
        % bin into 0.5 s bins over the habituation window
        counts(c, :) = histcounts(st, edges);            % 1 x num_bins
    end

    checkerboard_binned_spiking{mouse} = counts;   % [nCells x 600] for this mouse
end

mean_checkerboard_spiking = cell(num_mice, 1);

for mouse = mice_to_analyse
    temp = checkerboard_binned_spiking{mouse};
    temp = mean(temp, 2);
    mean_checkerboard_spiking{mouse} = temp;
end

mean_checkerboard_spiking_norm = cell(num_mice, 1);

checker_bin_size = 2;     % seconds
hab_bin_size     = 0.5;   % seconds

for mouse = mice_to_analyse
    check = checkerboard_binned_spiking{mouse};
    hab = mean_fr_habituation{mouse};
    check_Hz = check / checker_bin_size;  % spikes/s for each 2 s bin
    hab_Hz   = hab   / hab_bin_size;      % spikes/s baseline per cell
    resp_norm = check_Hz - hab_Hz;
    mean_checkerboard_spiking_norm{mouse} = resp_norm;
end
%%

for mouse = mice_to_analyse
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        evt_folders{mouse} = fullfile(mouse_path, mouse_name, ...
            session, 'Neural Data', 'Concatenated Data'); 
end


depth_tables = cell(num_mice, 1);

for mouse = mice_to_analyse
    temp = evt_folders{mouse};
    current_depth_struc = dir(fullfile(temp, 'mouse*_renewal_cell_depths.mat'));
    current_depths_paths = fullfile(current_depth_struc.folder, current_depth_struc.name);
    current_depth = load(current_depths_paths);
    temp = current_depth.depth_table;
    temp = sortrows(temp, {'y','x'}, 'descend');
    depth_tables{mouse} = temp;
end

sorted_spiking = cell(num_mice, 1);

for mouse = mice_to_analyse
    current_spiking = mean_checkerboard_spiking_norm{mouse};
    current_depth = depth_tables{mouse};
    depth_order = current_depth.CellID + 1;
    temp = current_spiking(depth_order);
    sorted_spiking{mouse} = temp;
end
