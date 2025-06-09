% Concatinating triggers


path1 = 'C:\Cohort 4 Temp Data Storage\Renewal\Neural Data\mouse1_renewal_p1\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\';
path2 = 'C:\Cohort 4 Temp Data Storage\Renewal\Neural Data\mouse1_renewal_p2\Record Node 101\experiment1\recording1\events\Neuropix-PXI-100.ProbeA\TTL\';

path1_continuous = 'C:\Cohort 4 Temp Data Storage\Renewal\Neural Data\mouse1_renewal_p1\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path2_continuous = 'C:\Cohort 4 Temp Data Storage\Renewal\Neural Data\mouse1_renewal_p2\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';

filepath_out = 'C:\Cohort 4 Temp Data Storage\Renewal\Neural Data\'

cont0 = 0;

evt = readNPY([path1 'timestamps.npy']); 
states = readNPY([path1 'states.npy']);
cont = readNPY([path1_continuous 'timestamps.npy']);
evt = evt-cont(1)+cont0;
evt = evt(states==1);
cont0 = cont0+numel(cont);
save([filepath_out 'evt_rec1'],'evt');

evt = readNPY([path2 'timestamps.npy']); 
states = readNPY([path2 'states.npy']);
cont = readNPY([path2_continuous 'timestamps.npy']);
evt = evt-cont(1)+cont0;
evt = evt(states==1);
cont0 = cont0+numel(cont);
save([filepath_out 'evt_rec2'],'evt');


