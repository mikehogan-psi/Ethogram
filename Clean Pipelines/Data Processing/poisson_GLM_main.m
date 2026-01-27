% POISSON_GLM_MAIN
%
% Preps neural and behavioural data for input into a GLM that predicts
% firing based on behavioural covariates and saves pseudo-R2 values and 
% input tables

%% Directory Setup

master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';
 
session = 'Extinction';

mice_to_analyse = [2 3 4 5 6 7 8 9];

received_stim_set_1 = [2 3 6 7];
received_stim_set_2 = [4 5 8 9];

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
n_trials  = 40; % total number of trials
fs = 30000;  % Ephys sampling rate (30kHz)

% Extract trigger file filepaths
evt_folders = cell(num_mice, 1);

for mouse = mice_to_analyse
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        evt_folders{mouse} = fullfile(mouse_path, mouse_name, ...
            session, 'Neural Data', 'Triggers'); 
end

% Load trigger files
evt_data = cell(num_mice, 2);
for mouse = mice_to_analyse
    current_evt_folder = evt_folders{mouse};  
    current_evt_file_1 = dir(fullfile(current_evt_folder, '*_p1_triggers.mat'));
    current_evt_file_2 = dir(fullfile(current_evt_folder, '*_p2_triggers.mat'));
    tmp1 = load(fullfile(current_evt_file_1.folder, current_evt_file_1.name));
    evt_data{mouse, 1} = tmp1;
    tmp2 = load(fullfile(current_evt_file_2.folder, current_evt_file_2.name));
    evt_data{mouse, 2} = tmp2;
end


% Extract trigger timings, convert to seconds, reshape into trials x frames
% and then add extra trigger at the end as edge for input to histc
evt_data_cat = cell(num_mice, 1);

for mouse = mice_to_analyse
    evt_data_cat{mouse} = cat(1, (evt_data{mouse, 1}.evt(2:end)), evt_data{mouse, 2}.evt(2:end));
    evt_data_cat{mouse}= double(evt_data_cat{mouse});
    evt_data_cat{mouse} = evt_data_cat{mouse}/fs;
    evt_data_cat{mouse} = reshape(evt_data_cat{mouse}, ...
        numel(evt_data_cat{mouse})/n_trials,n_trials)';
    current_mouse_evt = evt_data_cat{mouse}; % nTrials x nEdges (40 x 502)
    bin_duration = median(diff(current_mouse_evt, 1, 2), 2);  % (nTrials x 1) per-trial bin width
    % Append one extra edge per trial equal to last bin + avg bin duration
    % for that trial
    evt_data_cat{mouse} = [current_mouse_evt, current_mouse_evt(:, end) + bin_duration];  

end

% Extract kilosort filepaths
neural_data_folders = cell(num_mice, 1);

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

num_behaviours = size(behaviour_data_flashes, 2);

for mouse = mice_to_analyse

    mouse_name = mouse_files(mouse).name;
    mouse_name = lower(strrep(mouse_name, ' ', ''));
    save_name = [mouse_name, '_', lower(session), '_poisson_GLM_data_and_results.mat'];
    save_path = fullfile(data_folders{mouse}, save_name);
    
    if ~exist(save_path, 'file')
    
        X = [];
        y = [];
        cell_ids = [];
        % Select trigger timings for correct mouse
        evt = evt_data_cat{mouse};
        
        % Extract spike timings and cluster values for correct mouse
        spike_times = spike_matrix{mouse};
        clusters = cluster_id_matrix{mouse};
        cluster_values = cluster_value_matrix{mouse};
    
        for neuron = 1:numel(cluster_values)
           
    
            if ismember(mouse, received_stim_set_1)
                stim_set = 1;
            elseif ismember(mouse, received_stim_set_2)
                stim_set = 2;
            end
             
            spike_count = raster_full_trial(evt, spike_times, clusters, cluster_values, stim_set);
            
            n_trials = size(evt, 1);
    
            % Bin neural data
            spikes_binned = zeros(n_trials, num_bins);
    
            for trial = 1:n_trials
                for bin = 1:num_bins
                    start_idx = floor((bin-1)*bin_size*fps) + 1;
                    end_idx   = min(round(bin*bin_size*fps), frames_per_trial);
                    spikes_binned(trial, bin) = sum(spike_count{neuron}(trial, start_idx:end_idx));
                end
            end
    
            tmp_y = reshape(spikes_binned', [], 1);  % Flatten 
            
            binned_behaviour_all = zeros(num_behaviours, n_trials, num_bins);  % 5 x 40 x 66
    
            for trial_type = 1:2
             
                % Select correct behavioural data
                if trial_type == 1
                    behaviour_matrix = behaviour_data_looms(mouse, :);
                else
                    behaviour_matrix = behaviour_data_flashes(mouse, :);
                end
             
                frames_per_trial = size(behaviour_matrix{1}, 2);
        
                % Initialise binned behaviour matrix
                binned_behaviour = zeros(num_behaviours, n_trials/2, num_bins);
        
                for behaviour = 1:num_behaviours
                    for trial = 1:n_trials/2
                        for bin = 1:num_bins
                            start_idx = floor((bin-1)*bin_size*fps) + 1;
                            end_idx   = min(round(bin*bin_size*fps), frames_per_trial);
                            trial_data = behaviour_matrix{behaviour}(trial, :);
                            binned_behaviour(behaviour, trial, bin) = mean(trial_data(start_idx:end_idx));
                        end
                    end
                end
                
                if trial_type == 1
                    % Looms go first (trials 1..20)
                    binned_behaviour_all(:, 1:(n_trials/2), :) = binned_behaviour;
                else
                    % Flashes go second (trials 21..40)
                    binned_behaviour_all(:, (n_trials/2+1):n_trials, :) = binned_behaviour;
                end
    
                % Flatten and assemble features
                binned_grooming = reshape(squeeze(binned_behaviour_all(1, :, :))', [], 1);
                binned_rearing = reshape(squeeze(binned_behaviour_all(2, :, :))', [], 1);
                binned_darting = reshape(squeeze(binned_behaviour_all(3, :, :))', [], 1);
                binned_freezing = reshape(squeeze(binned_behaviour_all(4, :, :))', [], 1);
                binned_velocity = reshape(squeeze(binned_behaviour_all(5, :, :))', [], 1);
                
                binned_behaviour = cat(2, binned_grooming, binned_rearing, binned_darting, ...
                    binned_freezing, binned_velocity);
        
                time_bins = repmat((1:num_bins)', n_trials, 1);
                trials_binned = repelem((1:n_trials)', num_bins);
                neuron_number = repelem(neuron-1, n_trials * num_bins)';
        
                if trial_type == 1
                    trial_id = 1;
                else
                    trial_id = 0;
                end
                
                % Looms first (1), flashes second (0)
                trial_identifier = [ones((n_trials/2)*num_bins, 1); zeros((n_trials/2)*num_bins, 1)];
        
                % Create stimulus indicator columns
                stim_indicator_loom = zeros(size(binned_velocity, 1), 1);
                stim_indicator_flash = zeros(size(binned_velocity, 1), 1);
                
                % Stimulus occurs 152–202 frames (convert to seconds and then bins)
                stim_start_s = 152 / fps;  
                stim_end_s = 202 / fps;    
                stim_start_bin = ceil(stim_start_s / bin_size);
                stim_end_bin = ceil(stim_end_s / bin_size);
                
                % Extract stimulus timings
                for trial = 1:n_trials
                    % indices in flattened vector
                    trial_offset = (trial-1) * num_bins;
                    stim_bins = trial_offset + (stim_start_bin:stim_end_bin);
                    if trial <= n_trials/2   % Loom
                        stim_indicator_loom(stim_bins) = 1;
                    else                 % Flash
                        stim_indicator_flash(stim_bins) = 1;
                    end
                end
    
            end
    
            tmp_X = [binned_behaviour, time_bins, trials_binned, trial_identifier, stim_indicator_loom, stim_indicator_flash];
            tmp_cell_ids = neuron_number;
    
            X = [X; tmp_X];
            y = [y; tmp_y];
            cell_ids = [cell_ids; tmp_cell_ids];
    
        end
        
        % Create table for glm input (all rows, all cells)
        X_table = array2table(X, 'VariableNames', ...
            {'Grooming','Rearing','Darting','Freezing','Velocity', ...
             'TimeBin','Trial','TrialIdentifier','LoomON','FlashON'});
        
        % Scale only continuous predictors
        % X_table.Velocity = zscore(X_table.Velocity);
        % X_table.TimeBin  = zscore(X_table.TimeBin);
        % X_table.Trial    = zscore(X_table.Trial);
        
        % ---- Compute test pseudo-R2 per cell ----
        unique_cells = unique(cell_ids);
        nCells = numel(unique_cells);
        pseudoR2_test = nan(nCells,1);
        weights_per_cell = cell(nCells, 1);   % one entry per cell
        coef_names = {};                      % will fill once from first successful model
    
        
        opts = statset('glmfit');
        opts.MaxIter = 100;
        
        for c = 1:nCells
            this_cell = unique_cells(c);
        
            idx_all = (cell_ids == this_cell);
        
            X_cell = X_table(idx_all, :);
            y_cell = y(idx_all);
        
            % % Skip dead / almost-dead cells
            % if var(y_cell) == 0 || sum(y_cell) < 20
            %     pseudoR2_test(c) = NaN;
            %     continue
            % end
            
            % 80/20 train/test split
            N = numel(y_cell);
            cv = cvpartition(N, 'HoldOut', 0.2);
            train_idx = training(cv);
            test_idx  = test(cv);
        
            X_train = X_cell(train_idx, :);
            y_train = y_cell(train_idx);
        
            X_test  = X_cell(test_idx, :);
            y_test  = y_cell(test_idx);
        
                   % Fit GLM on train
            try
                mdl = fitglm(X_train, y_train, 'Distribution', 'poisson', 'Options', opts);
            catch
                pseudoR2_test(c) = NaN;
                weights_per_cell{c} = NaN;   % store NaN for alignment
                continue
            end
    
            % ---- Extract weights (coefficients) for this cell ----
            weights_per_cell{c} = mdl.Coefficients.Estimate;  % (Intercept + 10 predictors)
    
            % Save coefficient names once (same for all cells)
            if isempty(coef_names)
                coef_names = mdl.CoefficientNames;
            end
        
            % Deviance on TEST data
            lambda_test = predict(mdl, X_test);
            D_full_test = poissonDeviance(y_test, lambda_test);
    
        
            % Null deviance on TEST data (constant rate fit on TRAIN)
            mu_null = mean(y_train);
            mu_test_null = mu_null * ones(size(y_test));
            D_null_test = poissonDeviance(y_test, mu_test_null);
        
            pseudoR2_test(c) = 1 - D_full_test / D_null_test;


        end

        % Save all data in a table for later analysis
        X_table_full = X_table;
        X_table_full.SpikeCount = y;        % firing / spike counts
        X_table_full.CellID     = cell_ids; % neuron identifier per row
    
        % Save per-mouse results
        save(save_path, 'X_table_full', 'pseudoR2_test', ...
            'unique_cells', 'weights_per_cell', 'coef_names');
        disp([save_name ' saved to ' save_path '!'])
    
    else
        fprintf('%s already exists, skipping save\n', save_name);
    end

end

%%
function D = poissonDeviance(y, mu)
    % Compute Poisson deviance between observed y and mean mu
    % y, mu: column vectors of same length

    % Ensure column vectors
    y  = y(:);
    mu = mu(:);

    idx = y > 0;
    term_pos  = y(idx) .* log(y(idx) ./ mu(idx)) - (y(idx) - mu(idx));
    term_zero = mu(~idx);   % contribution when y=0

    D = 2 * (sum(term_pos) + sum(term_zero));
end
