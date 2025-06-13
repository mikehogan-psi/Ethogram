
fid_write = fopen('mouse2_concatenated_data_full.dat', 'w');

path1 = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort_4_06_05_25 (SC PAG Implanted Animals)\Renewal\Mouse 2\Neural Data\mouse2_renewal_habituation_2025-06-12_11-56-40\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path2 = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort_4_06_05_25 (SC PAG Implanted Animals)\Renewal\Mouse 2\Neural Data\mouse2_renewal_p1_2025-06-12_12-04-58\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path3 = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort_4_06_05_25 (SC PAG Implanted Animals)\Renewal\Mouse 2\Neural Data\mouse2_renewal_p2_2025-06-12_12-28-14\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path4 = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort_4_06_05_25 (SC PAG Implanted Animals)\Renewal\Mouse 2\Neural Data\mouse2_renewal_checkerboard_2025-06-12_12-52-28\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA';

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

fid_read = fopen([path4 'continuous.dat']);
A = fread(fid_read, '*int16');
fwrite(fid_write, A, 'int16');
fclose(fid_read);

fclose(fid_write);