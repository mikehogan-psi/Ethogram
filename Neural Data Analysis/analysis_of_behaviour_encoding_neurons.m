%% Analysis of Behaviour encoding neurons
% this script .... (give summary)

%% Plott Behaviour onset PSTH

%% Directory Setup

% Define folder path that contains neronal spiking data 
% (i.e. firing rate matrices from raster_NM and response groups obtained from neural_data_analysis_pipeline.m)
  neuronal_data_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\processed_data\spiking_data\';

% Define folder path that contains kilosorted dats of this session
    kilosort_dir = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Concatenated data\kilosort4\';  

% Define folder path that contains loom and flash event trigger timestamps
% (extracted events from triggers via neural_data_analysis pipeline)
    triggers_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\Triggers\';   

% Define folder path that contains behavioural labels (1x 20800 flattened
% freeze vector)
    behaviour_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\';

% Define folder path results shall be saved
    savepath = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\';



%% define variables

 % define target behaviour that shall be analysed
      behaviour = 'freezing';

 % define session 
      sesh = 'extinction';
  %   sesh = 'renewal';

 % define mouse

      mouse_nr = '2';



%% Loading neuronal spike, trial and behavioural data

% converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

    spk = double(spk)/30000;% converts sample number (by dividing through sampling rate) into seconds
    clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
    Ncell = numel(clu_val); % get number of detected clusters 

% load extracted loom and flash (i.e. stimuli onset) events
% Format: e.g. mouse3_extinction_extracted_events
load([triggers_path 'mouse' mouse '_' sesh '_' 'extracted_events'],'evt_loom', 'evt_flash')

    Ntrial = length(evt_loom);

% load and format behavioural labels 

    % load part 1
      load([ behaviour_path behaviour '_labels_mouse' mouse_nr '_' sesh '_p1.mat'], 'predicted_labels');
            labels_p1 = predicted_labels;
    % load part 2        
      load([ behaviour_path behaviour '_labels_mouse' mouse_nr '_' sesh '_p2.mat'], 'predicted_labels');
            labels_p2 = predicted_labels;
    % concatinate part 1 and part 2
      behaviour_labels = [labels_p1;labels_p2]; % [20080x1] - contains behaviour label (0 = no behaviour, 1 = behaviour) for every frame
    % reshape into [trials x frames]
      %behaviour_labels = reshape(behaviour_labels, 502, 40)';  % Transpose to get [40 x 502]

%% Identify mean duration of events

% Flatten to column vector
labels = behaviour_labels(:);

% Identify start and end of each behavioural event (runs of 1s)
is_behaviour = labels == 1;
d = diff([0; is_behaviour; 0]);  % pad with zeros to catch edges

onsets = find(d == 1);  % start of a run of 1s
offsets = find(d == -1) - 1;  % end of a run of 1s

% Compute durations in frames
durations_frames = offsets - onsets + 1;

% Convert to seconds
fps = 15;
durations_seconds = durations_frames / fps;

% Mean duration
mean_duration = mean(durations_seconds);

fprintf('Mean behavioural event duration: %.2f seconds\n', mean_duration);


 load([ behaviour_path behaviour '_labels_mouse' mouse '_' sesh '.mat'], 'predicted_labels');
 behaviour_labels = predicted_labels;
 behaviour_labels = reshape(behaviour_labels', [], 1);
    


%% identify behaviour onsets
% extact the frame numbers of behaviour onsets
% (filter out onsets which where preceeded by at least 3 seconds of non-behaviour before 
% and behaviour remained concistent for 3 seconds after onset)


% Constants
fps = 15;

% Variables:
tpre = 4;        % Time in sec before each event to include in the window (i.e. pre-stimulus time)
tpost = 4;       % Time in sec after each event to include.               (i.e. post-stimulus time)

    pre_frames =  tpre * fps;   % recalculate in frames
    post_frames = tpost * fps;  % recalculate in frames

valid_onsets = [];

for i = 1:length(onsets)
    onset = onsets(i);
    
    % Check bounds
    if onset <= pre_frames || onset + post_buffer - 1 > length(behaviour_labels)
        continue;
    end
    
    % Check 3 sec (45 frames) before are all 0s
    if all(behaviour_labels(onset - pre_frames : onset - 1) == 0)
        
        % Check 3 sec (45 frames) after are all 1s
        if all(behaviour_labels(onset : onset + post_buffer - 1) == 1)
            valid_onsets(end + 1) = onset; 
        end
    end
end


% convert frame numbet of behaviour onsets into seconds
  valid_onset_times = valid_onsets / fps;


%% Generating Firing rate matrices (spikes per frame aligned with behaviour onset)

% raster options (all in seconds)
    bin_size = 0.2;   % Bin size (time resolution of histogram)


% initialize variables of spike data to be saved matrices to hold firing rates [trials x times x cell] 
    fr_behav   = [];  % firing rates (Hz) over loom  trials [trials x time-bins x cell] 
    mfr_behav   = [];  % Mean firing rate (Hz) across  loom  trials [cell x time-bins]
    t_behav     = [];  % Time vector of all histogram timebins (excluding the first and last bin) - same foor loom and flash trials
    sfr_behav   = [];  % Standard error of firing rate across loom  trials [cell x time-bins]

% give list of kilosort clusters which shall be plotted
cluster = neuron_cluster_assignments(:,1); % MAKE CLUSTER VALUES STARTING FROM 0!!!  
                                        
% loop though all clusters to compute firing rated and load matrices

for n = 1: length(cluster)
=======

cluster = double(clu_val)+1;

for n = 1:length(cluster)

    clu_nr = cluster(n);

    tsp = spk(clu==clu_nr); % extracts spikes from this cluster

    [mfr_behav(n,:), sfr_behav(n,:), t_behav, fr_behav(:,:,n),  ~]   = raster_NM(tsp,valid_onset_times,tpre,tpost,bin_size,true,false);

    [mfr_loom(n,:), sfr_loom(n,:), t_loom, fr_loom(:,:,n),  ~]   = raster_NM(tsp,valid_onset_times,tpre,tpost,bin_size,false);


    fig = gcf;  
    sgtitle(sprintf(['mouse' mouse_nr ': ' behaviour ' ' sesh ' neuron %d'], clu_nr));

    %savefig(fig_flash, fullfile(figures_path, sprintf('Cluster_%d_FLASH.fig', clu_val(n))));
    %saveas(fig_flash, fullfile(figures_path, sprintf('Cluster_%d_FLASH.png', clu_val(n))));

    ginput(1); close all;

    % disp('Press any key or click to show next neuron...');
    % waitforbuttonpress;
    % close(fig);  % Optional: close current figure
end

% % save variables 
  save([savepath 'mouse_' mouse '_' sesh '_spikes_per_frame_loom_data' ],'mfr_loom', 'sfr_loom', 't_loom', 'fr_loom', 'valid_onset_times');
  % save([savepath 'mouse_' mouse '_' sesh 'spikes_per_frame_flash_data' ],'mfr_flash', 'sfr_flash', 't_flash', 'fr_flash');
% 





%% Separating behaviour into loom and flash trials
%  make sure you know which stimset the mouse you are analysis received!

stim_set_1 = [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1,...
    0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1];

stim_set_1 = stim_set_1(:);  % Ensure it's a column vector

% Get logical indices
loom_idx = stim_set_1 == 1;
flash_idx = stim_set_1 == 0;

% Extract trials
loom_labels = behaviour_labels(loom_idx, :);   % [20 x 502]
flash_labels = behaviour_labels(flash_idx, :); % [20 x 502]

%% Generating Firing rate matrices (spikes per frame aligned with stimulus onset)

% pre_stim_period = 2:152;
% during_stim_period = 153:201;
% post_stim_period = 202:502;


% raster options (all in seconds)
    tpre = 153/15;        % Time before each event to include in the window (i.e. pre-stimulus time)
    tpost = 350/15;       % Time after each event to include.               (i.e. post-stimulus time)
    bin_size = 1/15;   % Bin size (time resolution of histogram)


% initialize variables of spike data to be saved matrices to hold firing rates [trials x times x cell] 
    fr_loom   = [];  % firing rates (Hz) over loom  trials [trials x time-bins x cell] 
    fr_flash  = [];  % firing rates (Hz) over flash trials [trials x time-bins x cell] 
    mfr_loom  = [];  % Mean firing rate (Hz) across  loom  trials [cell x time-bins]
    mfr_flash = [];  % Mean firing rate (Hz) across  flash trials [cell x time-bins]
    t_loom    = [];  % Time vector of all histogram timebins (excluding the first and last bin) - same foor loom and flash trials
    t_flash   = [];  % Time vector of all histogram timebins (excluding the first and last bin) 
    sfr_loom  = [];  % Standard error of firing rate across loom  trials [cell x time-bins]
    sfr_flash = [];  % Standard error of firing rate across  flash trials [cell x time-bins]


% loop though all clusters to compute firing rated and load matrices
for n = 1:Ncell

    tsp = spk(clu==clu_val(n)); % extracts spike from one cluster/neuron
 
    [mfr_loom(n,:), sfr_loom(n,:), t_loom, fr_loom(:,:,n),  ~]   = raster_NM(tsp,evt_loom,tpre,tpost,bin_size,false);
    [mfr_flash(n,:),sfr_flash(n,:),t_flash,fr_flash(:,:,n), ~]   = raster_NM(tsp,evt_flash,tpre,tpost,bin_size,false);

end


% save variables 
  save([savepath 'mouse_' mouse_nr '_' sesh '_spikes_per_frame_loom_data' ],'mfr_loom', 'sfr_loom', 't_loom', 'fr_loom');
  save([savepath 'mouse_' mouse_nr '_' sesh 'spikes_per_frame_flash_data' ],'mfr_flash', 'sfr_flash', 't_flash', 'fr_flash');

% % save variables 
%   save([savepath 'mouse_' mouse '_' sesh '_spikes_per_frame_loom_data' ],'mfr_loom', 'sfr_loom', 't_loom', 'fr_loom');
%   save([savepath 'mouse_' mouse '_' sesh 'spikes_per_frame_flash_data' ],'mfr_flash', 'sfr_flash', 't_flash', 'fr_flash');


%% Plotting rasterplots of Loom trials with freezing behavior and firing rate (colormap version)

% Inputs
n = 5;  % Neuron index to plot
[nTrials, nFrames] = size(loom_labels);
gapSize = 1;  % We'll insert 1 color-coded row between each trial

% Extract firing rate for selected neuron [trials x bins]
fr_data = squeeze(fr_loom(:,:,n));  % now [trials x bins]

% Normalize firing rate to [0,1] for colormap
fr_norm = fr_data - min(fr_data(:));
fr_norm = fr_norm ./ max(fr_data(:));

% Ensure behavior and firing have same time resolution
minFrames = min(size(loom_labels, 2), size(fr_norm, 2));
loom_labels = loom_labels(:, 1:minFrames);
fr_norm = fr_norm(:, 1:minFrames);

% Prepare stacked matrix with behavior + firing rate rows
spacedMatrix = [];  % Will hold [behavior; firing; behavior; firing; ...]
rowType = [];       % Track row type: 0 = behavior, 1 = firing

for t = 1:nTrials
    spacedMatrix = [spacedMatrix; loom_labels(t, :)];
    rowType = [rowType; 0];  % behavior
    
    if t < nTrials
        spacedMatrix = [spacedMatrix; fr_norm(t, :) + 2];  % offset firing rows
        rowType = [rowType; 1];  % firing
    end
end

% Initialize RGB image
RGBimage = zeros(size(spacedMatrix,1), size(spacedMatrix,2), 3);

% --- Behavior Rows ---
behaviorRows = rowType == 0;
freezeMask = spacedMatrix == 1 & behaviorRows;
moveMask   = spacedMatrix == 0 & behaviorRows;
RGBimage(:,:,3) = freezeMask;  % Blue
RGBimage(:,:,1) = moveMask;    % Red

% --- Firing Rows (blue to red colormap) ---
firingRows = rowType == 1;
firingMask = spacedMatrix > 2 & firingRows;

% Extract normalized firing values back from spacedMatrix
firingVals = spacedMatrix(firingMask) - 2;  % [0,1]

% Map to colormap
cmap = jet(256);  % jet goes blue->green->red
idx = round(firingVals * 255) + 1;
rgbVals = cmap(idx, :);  % RGB triplets for each value

% Assign RGB values to image
[rIdx, cIdx] = find(firingMask);  % row/col indices
for i = 1:length(rIdx)
    RGBimage(rIdx(i), cIdx(i), :) = rgbVals(i, :);
end

% --- Plot ---
figure;
image(RGBimage);
axis tight;
xlabel('Time Bin');
ylabel('Trial');
title(sprintf('Freezing and Firing Rate (Neuron %d, Loom Trials)', n));

% Y-axis ticks for each trial (behavior rows only)
trialTicks = find(rowType == 0);
set(gca, 'YTick', trialTicks);
set(gca, 'YTickLabel', 1:nTrials);
