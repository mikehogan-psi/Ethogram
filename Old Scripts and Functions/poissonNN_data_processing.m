%% Directory Setup

master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';
 
session = 'Extinction';

mice_to_analyse = [2 3 4 5 6 7];

%% Load behavioural data

behaviours = {'Grooming', 'Rearing', 'Darting'};

trial_types = {'looms', 'flashes'};

% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

num_mice = numel(mice_to_analyse);
num_behaviours = numel(behaviours);

behaviour_folders = cell(num_mice, num_behaviours);

for mouse = mice_to_analyse
    for behaviour = 1:num_behaviours
        current_behaviour = behaviours{behaviour};
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        behaviour_folders{mouse, behaviour} = fullfile(mouse_path, mouse_name, ...
            session, 'Behavioural Data', 'Extracted Behaviours', current_behaviour);           
    end
end


% Load complex behaviours
behaviour_data_looms = cell(num_mice, num_behaviours);
behaviour_data_flashes = cell(num_mice, num_behaviours);

for mouse = mice_to_analyse
    for behaviour = 1:num_behaviours
        for trial_type = 1:numel(trial_types)
            current_trial_type = trial_types{trial_type};
            current_behaviour_path = behaviour_folders{mouse, behaviour};
            if strcmp(current_trial_type, 'looms')
                current_behaviour_file = dir(fullfile(current_behaviour_path, ...
            ['*_', behaviours{behaviour}, '_probabilities_looms.mat']));            
            tmp =  load(fullfile(current_behaviour_file.folder, ...
                current_behaviour_file.name), 'looms_scores');
            behaviour_data_looms{mouse, behaviour} = tmp.looms_scores;
            elseif strcmp(current_trial_type, 'flashes')
                current_behaviour_file = dir(fullfile(current_behaviour_path, ...
            ['*_', behaviours{behaviour}, '_probabilities_flashes.mat']));
            tmp =  load(fullfile(current_behaviour_file.folder, ...
                current_behaviour_file.name), 'flashes_scores');
            behaviour_data_flashes{mouse, behaviour} = tmp.flashes_scores;
            end
        end
    end
end

% Load freezing data
freeze_folders = cell(num_mice, 1);
freezing_data_looms = cell(num_mice, 1);
freezing_data_flashes = cell(num_mice, 1);

for mouse = mice_to_analyse
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    freeze_folders{mouse} = fullfile(mouse_path, mouse_name, ...
                 session, 'Behavioural Data', 'Extracted Behaviours', 'Freezing');
end

for mouse = mice_to_analyse 
    freeze_folder = freeze_folders{mouse};  
    % Looms
    current_freezing_paths = dir(fullfile(freeze_folder, '*_looms_freezing.mat'));
    tmp = load(fullfile(current_freezing_paths.folder, current_freezing_paths.name), 'loom_freezing');
    freezing_data_looms{mouse} = tmp.loom_freezing;
    
    % Flashes
    current_freezing_paths = dir(fullfile(freeze_folder, '*_flashes_freezing.mat'));
    tmp = load(fullfile(current_freezing_paths.folder, current_freezing_paths.name), 'flash_freezing');
    freezing_data_flashes{mouse} = tmp.flash_freezing;
end


% Load velocity data
velocity_folders = cell(num_mice, 1);
velocity_data_looms = cell(num_mice, 1);
velocity_data_flashes = cell(num_mice, 1);

for mouse = mice_to_analyse
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    velocity_folders{mouse} = fullfile(mouse_path, mouse_name, ...
                 session, 'Behavioural Data', 'Extracted Behaviours', 'Velocity');
end

for mouse = mice_to_analyse 
    velocity_folder = velocity_folders{mouse};  
    % Looms
    current_velocity_paths = dir(fullfile(velocity_folder, '*_looms_velocity.mat'));
    tmp = load(fullfile(current_velocity_paths.folder, current_velocity_paths.name), 'loom_velocity');
    velocity_data_looms{mouse} = tmp.loom_velocity;
    
    % Flashes
    current_velocity_paths = dir(fullfile(velocity_folder, '*_flashes_velocity.mat'));
    tmp = load(fullfile(current_velocity_paths.folder, current_velocity_paths.name), 'flash_velocity');
    velocity_data_flashes{mouse} = tmp.flash_velocity;
end

% Concatenate velocity, position and freezing with other behaviours
% Final order of behaviours is grooming, rearing, darting, freezing
behaviour_data_looms = cat(2, behaviour_data_looms, freezing_data_looms, velocity_data_looms);
behaviour_data_flashes = cat(2, behaviour_data_flashes, freezing_data_flashes, velocity_data_flashes);

%% Load neural data
N_trials  = 20; % total number of trials
fs = 30000;  

% Extract event files
evt_folders = cell(num_mice, 1);

for mouse = mice_to_analyse
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        evt_folders{mouse} = fullfile(mouse_path, mouse_name, ...
            session, 'Neural Data', 'Triggers'); 
end

evt_data = cell(num_mice, 1);
for mouse = mice_to_analyse
    current_evt_folder = evt_folders{mouse};  
    current_evt_file = dir(fullfile(current_evt_folder, '*_extracted_events.mat'));
    tmp = load(fullfile(current_evt_file.folder, current_evt_file.name));
    evt_data{mouse} = tmp;
end

           
neural_data_folders = cell(num_mice, 1);

% Extract kilosort filepaths
for mouse = mice_to_analyse
    for behaviour = 1:num_behaviours
        current_behaviour = behaviours{behaviour};
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        neural_data_folders{mouse} = fullfile(mouse_path, mouse_name, ...
            session, 'Neural Data', 'Concatenated Data', 'kilosort4\');           
    end
end

% Load neural data from kilosort paths
spike_matrix = cell(num_mice, 1);
cluster_value_matrix = cell(num_mice, 1);

for mouse = mice_to_analyse
    clu = readNPY(fullfile(neural_data_folders{mouse}, 'spike_clusters.npy'));
    spk = readNPY(fullfile(neural_data_folders{mouse}, 'spike_times.npy'));
    spk = double(spk) / fs;  
    clu_val = unique(clu);  

    % Store both the spike times and the cluster labels
    spike_matrix{mouse} = spk;
    cluster_id_matrix{mouse} = clu;
    cluster_value_matrix{mouse} = double(clu_val);
end


%% Generate save location for data

% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

data_folders = cell(length(mice_to_analyse), 1);

% Get behaviour folders for specified session
for mouse = mice_to_analyse
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    data_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Combined Data');
end

%% Create .mat files for input to Poisson neural network

fps = 15;
frames_per_trial = 502;
trial_length = frames_per_trial / fps; % 33.47 s

bin_size = 0.5; % 0.5 s bins
num_bins = floor(trial_length / bin_size); % 66 bins

t_start = -10;  % Align to trial start (10s before event)
t_end = trial_length-abs(t_start);
period_length = abs(t_start) + abs(t_end);
num_behaviours = size(behaviour_data_flashes, 2);
neuron_counter = -1; % Start at -1 to preseve python-based numbering

for mouse = mice_to_analyse
    evts{1} = evt_data{mouse}.evt_loom;
    evts{2} = evt_data{mouse}.evt_flash;

    X = [];
    y = [];
    cell_ids = [];

    spk = spike_matrix{mouse};
    clu = cluster_id_matrix{mouse};
    clu_val = cluster_value_matrix{mouse};

    for neuron = 1:numel(clu_val)
        neuron_counter = neuron_counter + 1;

        for trial_type = 1:numel(evts)
            tsp = spk(clu == clu_val(neuron));

            [spike_counts, bin_edges, bin_centres] = raster_binned(tsp, evts{trial_type}, t_start, t_end, num_bins);
            spike_count = reshape(spike_counts', [], 1);

            % Select correct behavioural data
            if trial_type == 1
                behaviour_matrix = behaviour_data_looms(mouse, :);
            else
                behaviour_matrix = behaviour_data_flashes(mouse, :);
            end

            N_trials = size(behaviour_matrix{1}, 1);
            frames_per_trial = size(behaviour_matrix{1}, 2);

            % Initialise binned behaviour matrix
            binned_behaviour = zeros(num_behaviours, N_trials, num_bins);

            for behaviour = 1:num_behaviours
                for trial = 1:N_trials
                    for bin = 1:num_bins
                        start_idx = floor((bin-1)*bin_size*fps) + 1;
                        end_idx   = min(round(bin*bin_size*fps), frames_per_trial);
                        trial_data = behaviour_matrix{behaviour}(trial, :);
                        binned_behaviour(behaviour, trial, bin) = mean(trial_data(start_idx:end_idx));
                    end
                end
            end

            % Flatten and assemble features
            binned_grooming = reshape(squeeze(binned_behaviour(1, :, :))', [], 1);
            binned_rearing = reshape(squeeze(binned_behaviour(2, :, :))', [], 1);
            binned_darting = reshape(squeeze(binned_behaviour(3, :, :))', [], 1);
            binned_freezing = reshape(squeeze(binned_behaviour(4, :, :))', [], 1);
            binned_velocity = reshape(squeeze(binned_behaviour(5, :, :))', [], 1);
            
            binned_behaviour = cat(2, binned_grooming, binned_rearing, binned_darting, ...
                binned_freezing, binned_velocity);

            time_bins = repmat((1:num_bins)', N_trials, 1);
            trials_binned = repelem((1:N_trials)', num_bins);
            neuron_number = repelem(neuron_counter, N_trials * num_bins)';

            if trial_type == 1
                trial_id = 1;
            else
                trial_id = 0;
            end
            trial_identifier = repelem(trial_id, N_trials * num_bins)';

            % Create stimulus indicator columns
            stim_indicator_loom = zeros(size(binned_velocity, 1), 1);
            stim_indicator_flash = zeros(size(binned_velocity, 1), 1);
            
            % Stimulus occurs 152–202 frames (convert to bins)
            stim_start_s = 152 / fps;  
            stim_end_s = 202 / fps;    
            stim_start_bin = ceil(stim_start_s / bin_size);
            stim_end_bin = ceil(stim_end_s / bin_size);
            
            % Extract stimulus timings
            for trial = 1:N_trials
                % indices in flattened vector
                trial_offset = (trial-1) * num_bins;
                stim_bins = trial_offset + (stim_start_bin:stim_end_bin);
                if trial_type == 1   % Loom
                    stim_indicator_loom(stim_bins) = 1;
                else                 % Flash
                    stim_indicator_flash(stim_bins) = 1;
                end
            end

            tmp_X = [binned_behaviour, time_bins, trials_binned, trial_identifier, stim_indicator_loom, stim_indicator_flash];
            tmp_cell_ids = neuron_number;
            tmp_y = spike_count;

            X = [X; tmp_X];
            y = [y; tmp_y];
            cell_ids = [cell_ids; tmp_cell_ids];
            
        end

    end

    % Transpose for correct Poisson NN format
    X = X'; 
    y = y'; 
    cell_ids = cell_ids';

    % Save
    mouse_name = mouse_files(mouse).name;
    save_folder = data_folders{mouse};
    save_name = ['mouse', num2str(mouse), '_poissonNN_prepped_data_', lower(session)];
    save(fullfile(save_folder, save_name), 'X', 'y', 'cell_ids', '-v7.3');
    disp([save_name, ' has been saved!']);

end
