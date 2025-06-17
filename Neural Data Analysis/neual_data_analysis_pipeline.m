%% Neural Analysis pipeline

%% Directoty setup

% Define folder path that contains neuropixel data (open ephyis output) of this session
ephys_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\Extinction'; 

% Define folder path that contains kilosorted dats of this session
kilosort_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\conactenated_data\kilosort4\'; % 

% create output folders where processed data will be saved
    filepath_out   = 'Z:\Abi\neuronal_data\mouse_2\processed_data\';  
    triggers_path =  [filepath_out 'concatinated_triggers\'];
    spikes_path    = [filepath_out 'spiking_data\'];
    figures_path   = [filepath_out 'spiking_data\figures\'];
     
    % Create folders if they do not exist
    folders = {triggers_path, spikes_path, figures_path};
    for i = 1:length(folders)
        if ~exist(folders{i}, 'dir')
            mkdir(folders{i});
        end
        addpath(folders{i});
    end

%% Setup variables

% define which mouse
   mouse = 'mouse2';

% define which session
    sesh = 'extinction';
  % sesh = 'renewal';

% define different session parts which are included in the continuous recording
    sesh_parts = {'hab', 'p1', 'p2'}; 


% raster options (all in seconds)
    tpre = 10;   % Time before each event to include in the window
    tpost = 20;  % Time after each event to include.
    bin_size = 0.2;    % Bin size (time resolution of histogram)

% general options
    fs = 30000;  % neuronal data sampling rate   
    Ntrial = 20; % number of trials in each session part
    onset = 153; % trigger nr at which stimulus starts (within each trial)


%% concatinating event triggers

  cont0 = 0; % sets start at 0 triggers
    
  for p = 1:length(sesh_parts) % loop through all seshion parts

    % get filepaths to open_ephys data   
    file_path = dir(fullfile(ephys_dir, ['*' sesh_parts{p} '*']));
    TTL_path  = [ephys_dir '\' file_path(1).name '\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\'];
    cont_path = [ephys_dir '\' file_path(1).name '\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\'];
    
    % convert python into matlab files 
    evt    = readNPY([TTL_path   'sample_numbers.npy']);  % extracts trigger counts (one at each video frame)
    states = readNPY([TTL_path   'states.npy']);          % 1 or -1 (TTL signal on/off)
    cont   = readNPY([cont_path  'sample_numbers.npy']);  % extracts neuronal data sample counts
    
    % adjust trigger events so that align with neuronal data samples
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);

    % save aligned triggers
    save([triggers_path mouse '_' sesh '_' sesh_parts{p}],'evt');
 
  end

%% extraxting relevant events for aligning spikes 

evt_flash = [];
evt_loom = [];

for p = 1:length(sesh_parts) % loop through all seshion parts

    load([triggers_path mouse '_' sesh '_' sesh_parts{p}],'evt'); 
    evt = double(evt)/fs;
    evt = evt(2:end);
    
    % extract random events from habituation (no stimuli)
    if p == 1
    evt_hav = sort_evts_hab(evt, tpre, tpost);
    
    % extract loom and flash events (i.e. triggers of stimulus onset)
    elseif p == 2 || p == 3    
    [evt_loom1,evt_flash1] = sort_evts(evt,Ntrial,onset,p);

    evt_flash = cat(1,evt_flash,evt_flash1);
    evt_loom =  cat(1,evt_loom ,evt_loom1);

    end 
   
end 

%% Loading neuronal spike data

% converting kilosort python output files into matlab format
clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

spk = double(spk)/fs;   % converts sample number (by dividing through sampling rate) into seconds
clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
Ncell = numel(clu_val); % get number of detected clusters 

%% Neuronal spike analysis for Loom and Flash trials 

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
  save([spikes_path mouse '_' sesh '_loom_data' ],'mfr_loom', 'sfr_loom', 't_loom', 'fr_loom');
  save([spikes_path mouse '_' sesh '_flash_data' ],'mfr_flash', 'sfr_flash', 't_flash', 'fr_flash');

%% Plott Rasters and PSTHs for any cluser/neuron

cluster = 1:5; % change

for n = cluster
    tsp = spk(clu==clu_val(n)); % extracts spike from this cluster

    % === LOOM ===
    [mfr_loom, ~, ~, fr_loom] = raster_NM(tsp, evt_loom, tpre, tpost, bin_size, true);
    fig_loom = gcf;  
    sgtitle(sprintf('Cluster %d — LOOM', clu_val(n)));
    savefig(fig_loom, fullfile(figures_path, sprintf('Cluster_%d_LOOM.fig', clu_val(n))));
    %saveas(fig_loom, fullfile(figures_path, sprintf('Cluster_%d_LOOM.png', clu_val(n))));
    
    % === FLASH ===
    [mfr_flash, ~, ~, fr_flash] = raster_NM(tsp, evt_flash, tpre, tpost, bin_size, true);
    fig_flash = gcf; 
    sgtitle(sprintf('Cluster %d — FLASH', clu_val(n)));
    savefig(fig_flash, fullfile(figures_path, sprintf('Cluster_%d_FLASH.fig', clu_val(n))));
    %saveas(fig_flash, fullfile(figures_path, sprintf('Cluster_%d_FLASH.png', clu_val(n))));
end

%% comparing mean firing rates before, during and after stimulus (LOOM trials)

% load flash and loom trial mean firing rate (across trials)
load([spikes_path mouse '_' sesh '_loom_data'], 'fr_loom', 't_loom');

n_cells  = size(fr_loom, 3); % extract number of clusters


% define furing and after stimulus time epochs 
idx_pre    = find(t_loom < 0);
idx_during = find(t_loom >= 0 & t_loom <= 3.3);
idx_post   = find(t_loom > 3.3);

% get Trial-by-Trial Spike Counts 
n_spk_pre    = squeeze(sum(fr_loom(:, idx_pre, :) * bin_size, 2));     % [nTrials x nCells]
n_spk_during = squeeze(sum(fr_loom(:, idx_during, :) * bin_size, 2));  % [nTrials x nCells]
n_spk_post   = squeeze(sum(fr_loom(:, idx_post, :) * bin_size, 2));    % [nTrials x nCells]

% Convert to Hz (to account for differnt lengths of epochs)
fr_pre    = n_spk_pre    / tpre;
fr_during = n_spk_during / 3.3;
fr_post   = n_spk_post   / tpost;

% Run signrank test on pre-stim vs stim and post-stim period firing rates 
pvals_pre_vs_during = zeros(n_cells, 1);
pvals_pre_vs_post   = zeros(n_cells, 1);

for c = 1:n_cells
    pvals_pre_vs_during(c) = signrank(fr_during(:,c), fr_pre(:,c));  
    pvals_pre_vs_post(c)   = signrank(fr_post(:,c),   fr_pre(:,c));
end

% Identify Responsive clusters
alpha = 0.05;

responsive_during = pvals_pre_vs_during < alpha; 
resp_duringLOOM_nr    = find(responsive_during == 1); 

responsive_post   = pvals_pre_vs_post < alpha;
resp_postLOOM_nr      = find(responsive_post == 1); 

responsive_during_and_post = (pvals_pre_vs_during < alpha) & (pvals_pre_vs_post < alpha);
resp_during_and_postLOOM_nr = find(responsive_during_and_post == 1);

% Report number of neurons that are signifcantly responsive
fprintf('Responsive neurons during loom: %d / %d\n', sum(responsive_during), n_cells);
fprintf('Responsive neurons after loom : %d / %d\n', sum(responsive_post), n_cells);
fprintf('Responsive neurons during and after loom : %d / %d\n', sum(responsive_during_and_post), n_cells);

% save responsive nerons
  save([spikes_path mouse '_' sesh '_LOOMresp_neurons' ],'resp_duringLOOM_nr', 'resp_postLOOM_nr', 'resp_during_and_postLOOM_nr');


%% comparing mean firing rates before, during and after stimulus (FLASH trials)

% load flash and loom trial mean firing rate (across trials)
load([spikes_path mouse '_' sesh '_flash_data'], 'fr_flash', 't_flash');

n_cells  = size(fr_flash, 3); % extract number of clusters

% define furing and after stimulus time epochs 
idx_pre    = find(t_flash < 0);
idx_during = find(t_flash >= 0 & t_flash <= 3.3);
idx_post   = find(t_flash > 3.3);

% get Trial-by-Trial Spike Counts 
n_spk_pre    = squeeze(sum(fr_flash(:, idx_pre, :) * bin_size, 2));     % [nTrials x nCells]
n_spk_during = squeeze(sum(fr_flash(:, idx_during, :) * bin_size, 2));  % [nTrials x nCells]
n_spk_post   = squeeze(sum(fr_flash(:, idx_post, :) * bin_size, 2));    % [nTrials x nCells]

% Convert to Hz (to account for differnt lengths of epochs)
fr_pre    = n_spk_pre    / tpre;
fr_during = n_spk_during / 3.3;
fr_post   = n_spk_post   / tpost;

% Run signrank test on pre-stim vs stim and post-stim period firing rates 
pvals_pre_vs_during = zeros(n_cells, 1);
pvals_pre_vs_post   = zeros(n_cells, 1);

for c = 1:n_cells
    pvals_pre_vs_during(c) = signrank(fr_during(:,c), fr_pre(:,c));  
    pvals_pre_vs_post(c)   = signrank(fr_post(:,c),   fr_pre(:,c));
end

% Identify Responsive clusters
alpha = 0.05;

responsive_during = pvals_pre_vs_during < alpha; 
resp_duringFLASH_nr    = find(responsive_during == 1); 

responsive_post   = pvals_pre_vs_post < alpha;
resp_postFLASH_nr      = find(responsive_post == 1); 

responsive_during_and_post = (pvals_pre_vs_during < alpha) & (pvals_pre_vs_post < alpha);
resp_during_and_postFLASH_nr = find(responsive_during_and_post == 1);

% Report number of neurons that are signifcantly responsive
fprintf('Responsive neurons during flash: %d / %d\n', sum(responsive_during), n_cells);
fprintf('Responsive neurons after flash : %d / %d\n', sum(responsive_post), n_cells);
fprintf('Responsive neurons during and after flash : %d / %d\n', sum(responsive_during_and_post), n_cells);

% save responsive nerons
  save([spikes_path mouse '_' sesh '_FLASHresp_neurons' ],'resp_duringFLASH_nr', 'resp_postFLASH_nr', 'resp_during_and_postFLASH_nr');

  %% do the same for renewal!!!
