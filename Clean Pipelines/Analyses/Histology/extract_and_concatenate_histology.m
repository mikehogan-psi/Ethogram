%% EXTRACT_AND_CONCATENATE_HISTOLOGY
%
% Calls histology and cell depth data extraction functions to load data and
% then joins them to give each cell a corresponding brain region based on
% its assigned channel and where each channel is has been found in the brain
%
%% Provide directories and session info
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';
histology_folder = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Histology\Final Output';
session = 'Renewal';

%% Extract histology data and cell position data
[hist_data] = extract_histology_data(histology_folder);
[depth_data] = extract_cell_depth_data(master_directory, session);

%% Generate save paths

mouse_dirs = dir(fullfile(master_directory, 'Mouse*'));

num_mice = length(mouse_dirs);

mouse_nos = zeros(num_mice, 1);
expression = '\d+';

for i = 1:num_mice
    mouse_no = regexp(mouse_dirs(i).name, expression, 'match', 'once');
    mouse_nos(i) = str2double(mouse_no);
end

[mouse_nos, sort_idx] = sort(mouse_nos);
mouse_dirs = mouse_dirs(sort_idx);

depth_paths = cell(num_mice, 1);

for i = 1:num_mice
    current_depth_path = fullfile(mouse_dirs(i).folder, mouse_dirs(i).name, session,...
        'Neural Data', 'Cell Depths');
    depth_paths{i} = current_depth_path;
end


%% Extract mouse numbers from hist and depth data to match them up
hist_numbers = zeros(length(hist_data), 1);
depth_numbers = zeros(length(hist_data), 1);

for i = 1:length(hist_data)
    current_hist = hist_data{i};
    if isempty(current_hist)
        hist_numbers(i) = NaN;
        continue
    else
    hist_numbers(i) = current_hist.mouse(1);
    end
end

for i = 1:length(depth_data)
    current_depth = depth_data{i};
    if isempty(current_depth)
        depth_numbers(i) = NaN;
        continue
    else
    depth_numbers(i) = current_depth.mouse(1);
    end
end

[true_false, locations] = ismember(depth_numbers, hist_numbers);

depth_idx_match = find(true_false);
hist_idx_match = locations(true_false);

%% Load tables and append 'region' column from hist data to depth data

appended_data = cell(length(hist_data), 1);

for current_mouse = 1:length(hist_idx_match)
    % Use matching idxs to get data from same mouse
    hist_idx = hist_idx_match(current_mouse);
    depth_idx = depth_idx_match(current_mouse);
    current_hist_data = hist_data{hist_idx};
    current_depth_data = depth_data{depth_idx};
    % Extract channel assigned to each cell
    cells_to_analyse = current_depth_data.Channel;
    
    for i = 1:length(cells_to_analyse)
        % Extract channel ID of current cell
        current_channel_id = cells_to_analyse(i);
        % Find idx of this channel ID in hist data
        region_idx = (current_hist_data.channel_id == current_channel_id);
        % Use idx to extract region associated with channel ID from hist
        % data
        current_region = current_hist_data.acronym(region_idx);
        current_depth_data.brain_region(i) = current_region;
    end

    appended_data{current_mouse} = current_depth_data;

end

%% Save data

for i = 1:length(appended_data)
    to_save = appended_data{i};
    mouse_no = to_save.mouse(1);
    save_folder = depth_paths{mouse_no};
    save_name = (['mouse' num2str(mouse_no), '_' lower(session), '_cell_brain_positions.mat']);
    save_path = fullfile(save_folder, save_name);
    if exist(save_path, 'file')
        warning([save_name, ' already exists, skipping save']);
    else
    save(save_path, 'to_save');
    disp([save_name, ' saved to ', save_folder]);
    end
end

