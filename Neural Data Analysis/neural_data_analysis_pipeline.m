%% Neural Analysis pipeline

%% Directory setup

% Define folder path that contains neuropixel data (open ephyis output) of this session
ephys_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\Extinction'; 

% Define folder path that contains kilosorted dats of this session
kilosort_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\conactenated_data\kilosort4\'; % 

% create output folders where processed data will be saved
    filepath_out   = 'Z:\Abi\neuronal_data\mouse_2\processed_data\';  
    triggers_path =  [filepath_out 'concatenated_triggers\'];
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
    dt = 0.2;    % Bin size (time resolution of histogram)

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
clu_val = unique(clu);
Ncell = numel(clu_val); % get number of detected clusters (i.e. potential neuronal cells)

%% Neuronal spike analysis for Loom and Flash trials 


%        - mfr: Mean firing rate across trials (Hz).
%        - sfr: Standard error of firing rate across trials (Hz).
%        - t: Time vector (excluding the first and last bin).
%        - fr: Firing rate matrix (trials × time bins).
%        - frspk: Vector of spike times aligned to events (all trials concatenated).
%        - trialspk: Trial numbers for each spike in frspk.

% initialize variables of spike data to be saved matrices to hold firing rates [trials x times x cell] 
fr_loom   = [];  % firing rates over loom  trials [trials x times x cell] 
fr_flash  = [];  % firing rates over flash trials [trials x times x cell] 
mfr_loom  = [];
mfr_flash = [];
sfr_loom  = [];
sfr_flash = [];
trialspk_loom  = [];
trialspk_flash = [];

% loop though all clusters to compute firing rated and load matriced
% for n = 1:Ncell
%     tsp = spk(clu==clu_val(n));
%     %[mfr,sfr,t,fr,frspk,trialspk] = raster_NM(tsp,evt,tpre,tpost,dt,graph)
%     [mfr_loom(:,n), sfr_loom(:,n), t_loom, fr_loom(:,:,n),  trialspk_loom(:,:,n)]   = raster_NM(tsp,evt_loom,tpre,tpost,dt,false);
%     [mfr_flash(:,n),sfr_flash(:,n),t_flash,fr_flash(:,:,n), trialspk_flash(:,:,n)]  = raster_NM(tsp,evt_flash,tpre,tpost,dt,false);
% 
% end

for n = 1:Ncell
    tsp = spk(clu==clu_val(n));
    %[mfr,sfr,t,fr,frspk,trialspk] = raster_NM(tsp,evt,tpre,tpost,dt,graph)
    [~,~,t,fr_loom(:,:,n)] = raster_NM(tsp,evt_loom,tpre,tpost,dt,true);title('LOOM')
    [~,~,t,fr_flash(:,:,n)] = raster_NM(tsp,evt_flash,tpre,tpost,dt,true);title('FLASH')
    ginput(1); close all;
end

for n = 1:Ncell
    tsp = spk(clu==clu_val(n));
    %[mfr,sfr,t,fr,frspk,trialspk] = raster_NM(tsp,evt,tpre,tpost,dt,graph)
    [~,~,t,fr_loom(:,:,n)] = raster_NM(tsp,evt_loom,tpre,tpost,dt,false);
    [~,~,t,fr_flash(:,:,n)] = raster_NM(tsp,evt_flash,tpre,tpost,dt,false);
end

