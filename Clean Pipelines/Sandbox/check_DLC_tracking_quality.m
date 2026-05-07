% !!!Specify which session is to be analysed!!!
% Acquisition, Extinction, or Renewal
session = 'Extinction';

% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

mouse_to_look_at = 2;
part_to_check = 'p1';

%% Get directory for DLC datafiles for each mouse for specified session
% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

% Initialise 
dlc_session_folders = cell(length(mouse_files), 1);

% Get DLC folders for specified session
for mouse = 1:length(mouse_files)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    dlc_session_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Raw DLC Data');
end

%% Sort into mouse order
% Extract mouse number from folder names
mouse_numbers = cellfun(@(x) sscanf(regexp(x, 'Mouse\s*(\d+)', 'match', 'once'), 'Mouse%d'), dlc_session_folders);

% Sort by mouse number and get index for sorted order
[~, sort_idx] = sort(mouse_numbers, 'ascend');

% Apply sorting to folder list
dlc_session_folders = dlc_session_folders(sort_idx);

% Optional: display sorted folder names
disp('Sorted Folder Names:');
disp(dlc_session_folders);
%%
% Master loop goes over each mouse and accesses its DLC data
for mouse = mouse_to_look_at
    current_dlc_folder_path = dlc_session_folders{mouse};
    dlc_file_list = dir(fullfile(current_dlc_folder_path, 'camera*_mouse*_p*.csv'));
    dlc_file_list = dlc_file_list(~[dlc_file_list.isdir]);
    base_names = cell(length(dlc_file_list), 1);
    dlc_file_names = cell(length(dlc_file_list), 1);
    
    % Extract file names from list of DLC files
    for n = 1:length(dlc_file_list)
        dlc_file_names{n} =  dlc_file_list(n).name;
    end
end
%%
part_mask = contains(dlc_file_names, part_to_check);
part_files = dlc_file_names(part_mask);

TH = 0.9;
L_all = [];

for c = 1:length(part_files)
    file_path = fullfile(current_dlc_folder_path, part_files{c});
    
    % Read CSV
    Tcsv = readtable(file_path);

    % Convert to array
    A = table2array(Tcsv);

    % DLC columns are usually [x y likelihood] repeated for each point
    likelihood_cols = 4:3:size(A,2);
    L_cam = A(:, likelihood_cols);   % frames x points

    % Stack into frames x points x cameras
    L_all(:,:,c) = L_cam;
end

%%
mean_L_frame = squeeze(median(L_all, [2 3], 'omitnan'));

figure;
plot(mean_L_frame, 'k')
hold on
yline(TH, '--r', 'TH')
xlabel('Frame')
ylabel('Median likelihood')
ylim([0 1])
title(sprintf('Mouse %d, %s: overall DLC quality', mouse, part_to_check))
%%
mean_L_cam = squeeze(median(L_all, 2, 'omitnan'));  % frames x cameras

figure;
plot(mean_L_cam, 'LineWidth', 1)
xlabel('Frame')
ylabel('Median likelihood')
title(sprintf('Mouse %d, %s: mean likelihood per camera', mouse, part_to_check))
legend(compose('Cam %d', 1:size(L_all,3)), 'Location', 'eastoutside')
% %%
% mean_L_point = squeeze(mean(L_all, 3, 'omitnan'));  % frames x points
% 
% figure;
% imagesc(mean_L_point')
% axis xy
% xlabel('Frame')
% ylabel('Point')
% title(sprintf('Mouse %d, %s: mean likelihood per point', mouse, part_to_check))
% colorbar
% caxis([0 1])
% 
% %%
% valid_cam_count = squeeze(sum(L_all > TH, 3));   % frames x points
% 
% figure;
% imagesc(valid_cam_count')
% axis xy
% xlabel('Frame')
% ylabel('Point')
% title(sprintf('Mouse %d, %s: cameras above threshold', mouse, part_to_check))
% colorbar