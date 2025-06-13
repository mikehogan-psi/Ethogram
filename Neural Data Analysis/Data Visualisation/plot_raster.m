function[]= plot_raster()

addpath('C:\Users\G71044MH\OneDrive - The University of Manchester\Documents\GitHub\Ethogram\Neural Data Analysis\Data Visualisation');

%general options
fs = 30000; % sampling rate (always 30,000)
Ntrial = 20; % number of trials
onset = 153; % which trigger time '0' of raster will be set to
%raster options
tpre = 10; % time (in seconds) before time 0 shown
tpost = 20; % time after
dt = 0.2; % bin size (in seconds)

path2 = 'C:\Cohort 4 Temp Data Storage\Mouse2\Extinction\Neural Data\mouse2_extinction_p1_2025-06-06_13-24-00\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\';
filepath_out = 'C:\Cohort 4 Temp Data Storage\Mouse2\Extinction\Neural Data';

load([filepath_out 'mouse2_p1'],'evt'); evt = double(evt)/fs;
[evt_loom1,evt_flash1] = sort_evts(evt,Ntrial,onset,1);
load([filepath_out 'mouse2_p2'],'evt'); evt = double(evt)/fs;
[evt_loom2,evt_flash2] = sort_evts(evt,Ntrial,onset,2);
evt_flash = cat(1,evt_flash1,evt_flash2);
evt_loom = cat(1,evt_loom1,evt_loom2);

% provide path to spike_clusters.npy from kilosort4 file
clu = readNPY('C:\Cohort 4 Temp Data Storage\Mouse2\Extinction\Neural Data\mouse2_extinction_conactenated_data\kilosort4\spike_clusters.npy');
% provide path to spike_times.npy from kilosort4 file
spk = readNPY('C:\Cohort 4 Temp Data Storage\Mouse2\Extinction\Neural Data\mouse2_extinction_conactenated_data\kilosort4\spike_times.npy');
spk = double(spk)/fs;
clu_val = unique(clu);
Ncell = numel(clu_val);

fr_loom = []; %trials x bins x neurons 
fr_flash = []; 
for n = 1:Ncell
    tsp = spk(clu==clu_val(n));
    %[mfr,sfr,t,fr,frspk,trialspk] = raster_NM(tsp,evt,tpre,tpost,dt,graph)
    [~,~,t,fr_loom(:,:,n)] = raster_NM(tsp,evt_loom,tpre,tpost,dt,true);title('LOOM')
    [~,~,t,fr_flash(:,:,n)] = raster_NM(tsp,evt_flash,tpre,tpost,dt,true);title('FLASH')
    ginput(1); close all;
end

%note abi: use fr and t from raster_NM