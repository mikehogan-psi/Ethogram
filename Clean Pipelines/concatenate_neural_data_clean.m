function [] = concatenate_neural_data_clean(hab_path, p1_path, p2_path, check_path, save_path, save_name)

fid_write = fopen(fullfile(save_path, save_name), 'w');

% Check if this is a extinction session, if so, skip checkerboard file
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

    clear A

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

    clear A

end

end
