
fid_write = fopen('mouse2_concatenated_data_full.dat', 'w');

path1 = 'C:\Cohort 4 Temp Data Storage\Mouse2\Neural Data\mouse2_extinction_habituation_2025-06-06_13-13-19\Record Node 101\experiment2\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path2 = 'C:\Cohort 4 Temp Data Storage\Mouse2\Neural Data\mouse2_extinction_p1_2025-06-06_13-24-00\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path3 = 'C:\Cohort 4 Temp Data Storage\Mouse2\Neural Data\mouse2_extinction_p2_2025-06-06_13-48-21\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';

fid_read = fopen([path1 'continuous.dat']);
A = fread(fid_read, '*int16');
fwrite(fid_write, A, 'int16');
fclose(fid_read);

fid_read = fopen([path2 'continuous.dat']);
A = fread(fid_read, '*int16');
fwrite(fid_write, A, 'int16');
fclose(fid_read);

fid_read = fopen([path3 'continuous.dat']);
A = fread(fid_read, '*int16');
fwrite(fid_write, A, 'int16');
fclose(fid_read);

fclose(fid_write);