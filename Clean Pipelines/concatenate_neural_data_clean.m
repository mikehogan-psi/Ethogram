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

% Open output file for writing
fid_write = fopen(fullfile(save_path, save_name), 'w');
if fid_write == -1
    error('Could not open output file for writing: %s', fullfile(save_path, save_name));
end

% List all input files
file_list = {hab_path, p1_path, p2_path};
if ~isempty(check_path)
    file_list{end+1} = check_path;
end

% Define chunk size (number of int16 elements to read at a time)
chunk_size = 1e8; %200MB at a time

% Loop through input files
for f = 1:numel(file_list)
    file_info = dir(file_list{f});
    total_bytes = file_info.bytes;
    processed_bytes = 0;
    
    fprintf('Concatenating %s (%.2f GB)...\n', file_list{f}, total_bytes/1e9);
    
    fid_read = fopen(file_list{f}, 'r');
    if fid_read == -1
        fclose(fid_write);
        error('Could not open input file for reading: %s', file_list{f});
    end
    
    while true
        A = fread(fid_read, chunk_size, '*int16'); % read a chunk
        if isempty(A)
            break;
        end
        fwrite(fid_write, A, 'int16');             % write chunk
        
        % Update progress
        processed_bytes = processed_bytes + numel(A)*2; % int16 = 2 bytes
        percent_done = (processed_bytes / total_bytes) * 100;
        fprintf('\r  Progress: %.1f%%', percent_done);
    end
    fprintf('\n'); % newline after each file finishes
    fclose(fid_read);
end

fclose(fid_write);
fprintf('Concatenation complete. File saved to %s\n', fullfile(save_path, save_name));

end
