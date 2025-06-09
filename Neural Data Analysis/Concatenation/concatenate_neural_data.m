
fid_write = fopen('test.dat', 'w');

path1 = 'C:\Cohort 4 Temp Data Storage\Renewal\Neural Data\mouse1_renewal_p1\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';
path2 = 'C:\Cohort 4 Temp Data Storage\Renewal\Neural Data\mouse1_renewal_p2\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA\';


fid_read = fopen([path1 'continuous.dat']);
A = fread(fid_read, '*int16');
fwrite(fid_write, A, 'int16');
fclose(fid_read);


fid_read = fopen([path2 'continuous.dat']);
A = fread(fid_read, '*int16');
fwrite(fid_write, A, 'int16');
fclose(fid_read);

fclose(fid_write);