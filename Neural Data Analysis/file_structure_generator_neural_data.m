% Define root directory
root_dir = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% Define folder names
mice = {'Mouse 1', 'Mouse 2', 'Mouse 3'};
sessions = {'Acquisition', 'Extinction', 'Retention'};
data_types = {'Behavioural Data', 'Neural Data'};
processed_unprocessed_b = {'SSM Fitted Data', 'Triangulated Data', 'Raw DLC Data', 'Video Data'};
processed_unprocessed_n = {'Raw Data', 'Concatenated Data', 'Triggers'};

% Create folder structure
for i = 1:length(mice)
    for j = 1:length(sessions)
        for k = 1:length(data_types)
            % Create main data type folder
            base_path = fullfile(root_dir, mice{i}, sessions{j}, data_types{k});
            if ~exist(base_path, 'dir')
                mkdir(base_path);
            end
            
            % If it's Behavioural Data, add relevant subfolders
            if strcmp(data_types{k}, 'Behavioural Data')
                for m = 1:length(processed_unprocessed_b)
                    sub_path_b = fullfile(base_path, processed_unprocessed_b{m});
                    if ~exist(sub_path_b, 'dir')
                        mkdir(sub_path_b);
                    end
                end
            end

            % If it's Neural Data, add relevant subfolders
            if strcmp(data_types{k}, 'Neural Data')
                for n = 1:length(processed_unprocessed_n)
                    sub_path_n = fullfile(base_path, processed_unprocessed_n{n});
                    if ~exist(sub_path_n, 'dir')
                        mkdir(sub_path_n);
                    end
                end
            end
        end
    end
end
