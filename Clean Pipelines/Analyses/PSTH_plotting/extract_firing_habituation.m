% EXTRACT_FIRING_HABITUATION
%
% Extract habituation data and plots it across whole session, showing raw
% and smoothed plots
% Also extract mean spike count per bin for each cell for use elsewhere
%% Directory Setup

master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';
 
session = 'Extinction';

mice_to_analyse = [2 3 4 5 6 7 8 9];
max_mouse_id    = max(mice_to_analyse); 

%% Load neural data

% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

num_mice = numel(mice_to_analyse);

n_trials  = 40; % total number of trials
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
    current_evt_file_1 = dir(fullfile(current_evt_folder, '*_habituation_triggers.mat'));
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


%% Habituation spiking in 0.5 s bins (per cell)

fps              = 15;
frames_per_trial = 4500;
trial_length     = frames_per_trial / fps;   % 300 s

bin_size = 0.5;                              % seconds
num_bins = trial_length / bin_size;          % 600 bins

% Store one matrix per mouse: [nCells x num_bins]
hab_binned_spiking = cell(num_mice, 1);

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

    hab_binned_spiking{mouse} = counts;   % [nCells x 600] for this mouse
end

%% Find mean spiking rate per cell across bins during habituation

mean_fr_habituation = cell(size(hab_binned_spiking));

for mouse = mice_to_analyse
    current_hab_data = hab_binned_spiking{mouse};
    current_mean_spiking = zeros(size(current_hab_data, 1), 1);
    for neuron = 1:size(current_hab_data, 1)
        current_mean_spiking(neuron, :) = mean(current_hab_data(neuron, :))/bin_size;
    end
    mean_fr_habituation{mouse} = current_mean_spiking;
end
        
%%

received_psi = [3 5 7 8];
received_veh = [2 4 6];

mean_per_bin_hab_fr = cell(num_mice, 1);

for mouse = mice_to_analyse
    mean_per_bin_hab_fr{mouse} = mean(hab_binned_spiking{mouse}, 1);
end


psi_mean_per_bin_hab_fr = mean_per_bin_hab_fr(received_psi);
veh_mean_per_bin_hab_fr = mean_per_bin_hab_fr(received_veh);

%%
% Turn cell arrays into [nMice x nBins] matrices
psi_mat = vertcat(psi_mean_per_bin_hab_fr{:});   % rows = mice, cols = bins
veh_mat = vertcat(veh_mean_per_bin_hab_fr{:});

% Mean and SEM across mice
mean_psi = mean(psi_mat, 1);
mean_veh = mean(veh_mat, 1);

sem_psi = std(psi_mat, 0, 1) ./ sqrt(size(psi_mat,1));
sem_veh = std(veh_mat, 0, 1) ./ sqrt(size(veh_mat,1));

nBins = size(psi_mat, 2);
bin_width = 0.5;                 % change if not 0.5 s per bin
t = (0:nBins-1) * bin_width;     % time axis

veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

veh_upper = mean_veh + sem_veh;
veh_lower = mean_veh - sem_veh;

psi_upper = mean_psi + sem_psi;
psi_lower = mean_psi - sem_psi;

figure('Color','w'); hold on;

% Vehicle shaded SEM
fill([t fliplr(t)], [veh_upper fliplr(veh_lower)], veh_col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t, mean_veh, 'Color', veh_col, 'LineWidth', 2);

% Psilocybin shaded SEM
fill([t fliplr(t)], [psi_upper fliplr(psi_lower)], psi_col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t, mean_psi, 'Color', psi_col, 'LineWidth', 2);

xlabel('Time in habituation (s)');
ylabel('Mean firing rate (Hz)');
title(['Mean Habituation firing rate: ' session]);
legend({'Veh SEM','Vehicle','Psi SEM','Psilocybin'}, 'Location','best');
ylim([1.5 5])

ax = gca;
ax.LineWidth = 1.5;
ax.TickDir   = 'out';
ax.Box       = 'off';

%%
window = 30;  % number of bins in the moving window, tweak this

mean_veh_smooth = smoothdata(mean_veh, 'movmean', window);
mean_psi_smooth = smoothdata(mean_psi, 'movmean', window);

figure('Color','w'); hold on;
plot(t, mean_veh_smooth, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
plot(t, mean_psi_smooth, 'LineWidth', 2, 'Color', [1 0 0]);
xlabel('Time (s)');
ylabel('Mean firing rate (Spike Count)');
title(['Mean Habituation FR (Smoothed: ' num2str(window*0.5) 's Window): ' session]);
legend({'Vehicle','Psilocybin'}, 'Location','best');
ylim([2 4.5])