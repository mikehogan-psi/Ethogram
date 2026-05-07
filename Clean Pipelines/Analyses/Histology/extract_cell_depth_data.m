%% EXTRACT_CELL_DEPTH_DATA
% Loads cell_position_extractor output and stores it in a cell array of tables
%
% INPUT
% master_directory: file path to where all data is stored
% session: 'Extinction' or 'Renewal'
%
% OUTPUT
% depth_data: cell array containing tables of cell depth data with mouse
% number appended

function[depth_data] = extract_cell_depth_data(master_directory, session)

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

depth_data = cell(num_mice, 1);

for mouse = 1:num_mice
    current_mouse_path = fullfile(mouse_dirs(mouse).folder, mouse_dirs(mouse).name);
    current_depth_folder = fullfile(current_mouse_path, session, 'Neural Data', 'Cell Depths');
    depth_dir = dir(fullfile(current_depth_folder, '*cell_depths.mat'));
    if isempty(depth_dir)
        warning(['Mouse ', num2str(mouse), ' has no cell depth data'])
        depth_data{mouse} = [];
    else
        current_depth_path = fullfile(depth_dir.folder, depth_dir.name);
        depth_data_temp = load(current_depth_path, 'depth_table');
        depth_data{mouse} = depth_data_temp.depth_table;
        depth_data{mouse}.mouse = repmat(mouse_nos(mouse), size(depth_data{mouse}, 1), 1);
    end
end


end
