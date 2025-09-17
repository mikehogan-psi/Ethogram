%% Directory Setup

% Define folder path that contains neronal spiking data 
% (i.e. firing rate matrices and response groups obtained from neural_data_analysis_pipeline.m)

% neuronal_data_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\processed_data\spiking_data\';

% Define folder path that contains kilosorted dats of this session
kilosort_dir = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\kilosort4\'; % 


triggers_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Triggers\';
  
%% Load data

% get important variables
N_trials  = 20; % total number of trials
fs = 30000;  % neuronal data sampling rate  
tpre = 10;        % Time before each event to include in the window (i.e. pre-stimulus time)
tpost = 23;       % Time after each event to include.               (i.e. post-stimulus time)

% neuronal data - converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

spk = double(spk)/fs;   % converts sample number (by dividing by sampling rate) into seconds
clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
N_neurons = numel(clu_val); % get number of detected clusters i.e neurons

% % load list of Loom-response-types 
%   load([neuronal_data_path 'mouse3_extinction_LOOMresp_neurons'], 'resp_duringLOOM_nr', 'resp_postLOOM_nr')

% load extracted loom and flash (i.e. stimuli onset) events
load([triggers_path 'mouse3_extinction_extracted_events'],'evt_flash', 'evt_loom')

% Load freezing data
load('C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\for riccardo\m3_flashes_freezing_extinction.mat')
load('C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\for riccardo\m3_looms_freezing_extinction.mat')


% % Create variable with all responsive neurons
% resp_LOOM_nr = [resp_duringLOOM_nr; resp_postLOOM_nr];
% resp_LOOM_nr = unique(resp_LOOM_nr);


% Specify which type of responsive neurons you want to look at and save
% name
% current_group = resp_LOOM_nr;
current_group = double(clu_val+1);
  %% Fit poisson model for during and post-stim periods for each neuron

period_length = 450;   
t_start = 10; %10 seconds before stimulus         
t_end = 20; % 20 seconds after stimulus  
num_bins = 50;
fps = 15;
evts = {evt_flash, evt_loom};
  
all_data_tables = cell(size(current_group,1), 2);

for trial_type = 1:length(evts)

    for n = 1:length(current_group) % loop through each neuron
        idx = current_group(n);
        disp(['Creating data table for neuron: ' num2str(idx-1)]);
    
        % Calculate binsize 
        bin_size = period_length / num_bins;  % get binsize dependent on length of investigated period (if longer, greater bin size required)
        bin_size_s = bin_size/15;             % gets binsize in seconds (15fps)
        disp(['Binsize is ' num2str(bin_size_s) ' seconds'])  % display binsize
        
        % Extract all spikes from current cluster/neuron
        tsp = spk(clu==clu_val(idx)); 
        
        % extract the spikes that are occuring in each timebin fo each trial
        [~, ~, t_event, fr, spikes]   = raster_NM(tsp,evts{trial_type},tpre,tpost,bin_size_s,false);
        
        % extract firing rates in time_bins that are occuring within neuronal response period
        t = t_event; % returned from raster_NM
        idx_during = find(t >= -t_start & t <= t_end);

        fr_period = fr(:, idx_during, :);
        
        spike_count = reshape(fr_period', [], 1); 
        
        % Select correct freezing data (starting with flash)
        if trial_type == 1
            freeze_matrix = mouse3_flashes_freezing;
        else
            freeze_matrix = mouse3_looms_freezing;
        end

        % Extract freezing within bin bounds
        frames_start = 1;
        frames_end = period_length;
        cropped_freezing = freeze_matrix(:, frames_start:frames_end);
       
        % Initialise binned freezing matrix
        binned_freezing = zeros(N_trials, num_bins);
        
        % Extract freeze data for each bin over every trial
        for trial = 1:size(cropped_freezing, 1)
            for bin = 1:num_bins
                start_idx = floor((bin-1)*bin_size) + 1;
                end_idx = round(bin*bin_size);
                binned_freezing(trial, bin) = mean(cropped_freezing(trial, start_idx:end_idx));
            end
        end
        
        % Flatten into column vector
        binned_freezing = reshape(binned_freezing', [], 1);

        
        % Get time bin, trial indices and neuron numbers
        time_bins = repmat((1:num_bins)', N_trials, 1);   % Time bin indices
        trials_binned = repelem((1:N_trials)', num_bins); % Trial numbers
        neuron_number = repelem(idx-1, N_trials * num_bins)'; % Neuron numbers
       
        % Get trial type for table      
        if trial_type == 1
            trial_id = 0; % flashes
        else
            trial_id = 1; % looms
        end
        trial_identifier = repelem(trial_id, N_trials * num_bins)';

        %Create data table
        data_binned = table(spike_count, binned_freezing, time_bins, trials_binned, trial_identifier, ...
                neuron_number, 'VariableNames', {'Firing', 'Freezing', 'TimeBin', 'Trial',...
                'TrialIdentifier', 'NeuronNumber'});
    
        all_data_tables{n, trial_type} = data_binned;
    end    

end

% Merge flash and loom data
merged_neurons = cell(size(all_data_tables, 1), 1); 

for i = 1:size(all_data_tables, 1)
    table1 = all_data_tables{i, 1};
    table2 = all_data_tables{i, 2};

    % Stack vertically 
    merged_neurons{i} = [table1; table2];  
end
