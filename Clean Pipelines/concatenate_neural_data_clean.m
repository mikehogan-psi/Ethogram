function [] = concatenate_neural_data_clean(hab_path, p1_path, p2_path, check_path, save_path, save_name)

cd(save_path)

fid_write = fopen(save_name, 'w');

if isempty(check_path)
    
    fid_read = fopen(hab_path);
    A = fread(fid_read, '*int16');
    fwrite(fid_write, A, 'int16');
    fclose(fid_read);
    
    fid_read = fopen(p1_path);
    A = fread(fid_read, '*int16');
    fwrite(fid_write, A, 'int16');
    fclose(fid_read);
    
    fid_read = fopen(p2_path);
    A = fread(fid_read, '*int16');
    fwrite(fid_write, A, 'int16');
    fclose(fid_read);
    
    fclose(fid_write);

else

    fid_read = fopen(hab_path);
    A = fread(fid_read, '*int16');
    fwrite(fid_write, A, 'int16');
    fclose(fid_read);
    
    fid_read = fopen(p1_path);
    A = fread(fid_read, '*int16');
    fwrite(fid_write, A, 'int16');
    fclose(fid_read);
    
    fid_read = fopen(p2_path);
    A = fread(fid_read, '*int16');
    fwrite(fid_write, A, 'int16');
    fclose(fid_read);    

    fid_read = fopen(check_path);
    A = fread(fid_read, '*int16');
    fwrite(fid_write, A, 'int16');
    fclose(fid_read);

    fclose(fid_write);

end

end
