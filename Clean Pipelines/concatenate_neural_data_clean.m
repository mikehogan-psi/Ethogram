function [] = concatenate_neural_data_clean(hab_path, p1_path, p2_path, check_path, save_path, save_name)
% CONCATENATE_NEURAL_DATA_CLEAN Concatenates multiple neural recording files.
%   This function concatenates binary neural data files (int16 format) into
%   a single output file. The order of concatenation is:
%       1. habituation file (hab_path)
%       2. part 1 file (p1_path)
%       3. part 2 file (p2_path)
%       4. optional checkerboard file (check_path) if provided (in renewal)
%  
%   Inputs:
%       hab_path   - full path to habituation binary file
%       p1_path    - full path to part 1 binary file
%       p2_path    - full path to part 2 binary file
%       check_path - full path to checkerboard binary file (empty if not used)
%       save_path  - directory where concatenated file will be saved
%       save_name  - name of the output file 
%
%   Output:
%       A concatenated binary file written to [save_path, save_name].
%
%   To call this function with all inputs preformatted, use
%   neural_data_preprocessing_pipeline.mat

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
