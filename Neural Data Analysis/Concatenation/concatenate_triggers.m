
path1 = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\mouse3_extinction_habituation\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\';
path2 = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\mouse3_extinction_p1\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\';
path3 = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\mouse3_extinction_p2\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\';
% path4 = 'C:\Cohort 4 Temp Data Storage\Mouse1\Extinction\Neural Data\mouse1_extinction_p3\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\';

path1_continuous = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\mouse3_extinction_habituation\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path2_continuous = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\mouse3_extinction_p1\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path3_continuous = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\mouse3_extinction_p2\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
% path4_continuous = 'C:\Cohort 4 Temp Data Storage\Mouse1\Extinction\Neural Data\mouse1_extinction_p3\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';

filepath_out = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\Neural Data\';

cont0 = 0;

evt = readNPY([path1 'sample_numbers.npy']); 
states = readNPY([path1 'states.npy']);
cont = readNPY([path1_continuous 'sample_numbers.npy']);
evt = evt-cont(1)+cont0;
evt = evt(states==1);
cont0 = cont0+numel(cont);
save([filepath_out 'mouse2_hab'],'evt');

evt = readNPY([path2 'sample_numbers.npy']); 
states = readNPY([path2 'states.npy']);
cont = readNPY([path2_continuous 'sample_numbers.npy']);
evt = evt-cont(1)+cont0;
evt = evt(states==1);
cont0 = cont0+numel(cont);
save([filepath_out 'mouse2_p1'],'evt');

evt = readNPY([path3 'sample_numbers.npy']); 
states = readNPY([path3 'states.npy']);
cont = readNPY([path3_continuous 'sample_numbers.npy']);
evt = evt-cont(1)+cont0;
evt = evt(states==1);
cont0 = cont0+numel(cont);
save([filepath_out 'mouse2_p2'],'evt');

% evt = readNPY([path4 'timestamps.npy']); 
% states = readNPY([path4 'states.npy']);
% cont = readNPY([path4_continuous 'timestamps.npy']);
% evt = evt-cont(1)+cont0;
% evt = evt(states==1);
% cont0 = cont0+numel(cont);
% save([filepath_out 'mouse1_p3'],'evt');