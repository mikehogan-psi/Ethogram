%% Load 
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';
session = 'Renewal';

%%
if strcmp(session, 'Extinction')
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ext\cluster_assignments_extinction.mat")
    cluster_assignments = cluster_assignments_extinction;
elseif strcmp(session, 'Renewal')
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ren\cluster_assignments_renewal_to_extinction.mat")
    cluster_assignments = cluster_assignments_renewal;
end

%% Load cluster data and spiking data
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

spiking_paths = cell(numel(mouse_dirs), 1);

for mouse = 1:num_mice
    current_mouse_folder = mouse_dirs(mouse).folder;
    current_mouse_name = mouse_dirs(mouse).name;
    current_mouse_path = fullfile(current_mouse_folder, current_mouse_name,...
        session, 'Combined Data');
    spiking_paths{mouse} = current_mouse_path;
end

combined_data_matrices = cell(numel(mouse_dirs), 1);

for mouse = 1:num_mice
    current_file_to_load = dir(fullfile(spiking_paths{mouse}, '*GLM_data*'));
    if isempty(current_file_to_load)
        combined_data_matrices{mouse} = [];
        warning(['Mouse ', num2str(mouse), ' has no GLM data'])
    else 
        current_data = load(fullfile(current_file_to_load.folder, current_file_to_load.name));
        combined_data_matrices{mouse} = current_data;
    end
end
%% Load hist data
hist_paths = cell(numel(mouse_dirs), 1);

for mouse = 1:num_mice
    current_mouse_folder = mouse_dirs(mouse).folder;
    current_mouse_name = mouse_dirs(mouse).name;
    current_mouse_path = fullfile(current_mouse_folder, current_mouse_name,...
        session, 'Neural Data', 'Cell Depths');
    hist_paths{mouse} = current_mouse_path;
end

combined_hist_data = cell(numel(mouse_dirs), 1);

for mouse = 1:num_mice
    current_file_to_load = dir(fullfile(hist_paths{mouse}, '*brain_positions*'));
    if isempty(current_file_to_load)
        combined_hist_data{mouse} = [];
        warning(['Mouse ', num2str(mouse), ' has no hist data'])
    else 
        current_data = load(fullfile(current_file_to_load.folder, current_file_to_load.name));
        current_data = current_data.to_save;
        combined_hist_data{mouse} = current_data;
    end
end
    
%%

hist_and_spiking_matrices = cell(num_mice, 1);

for mouse = 1:length(combined_data_matrices)
    if isempty(combined_data_matrices{mouse}) || isempty(combined_hist_data{mouse})
        hist_and_spiking_matrices{mouse} = [];
        warning(['Mouse ' num2str(mouse), ' is missing data'])
        continue
    else
        current_cells_idx = cluster_assignments.MouseID == mouse;
        current_cells = cluster_assignments(current_cells_idx, :);
        current_spiking = combined_data_matrices{mouse}.X_table_full;
        keep_idx = ismember(current_spiking.CellID, current_cells.CellID);
        good_spiking = current_spiking(keep_idx, :);
        current_hist = combined_hist_data{mouse};
        good_spiking.ClusterID = repelem(NaN, height(good_spiking), 1);
        good_spiking.Brain_Region = strings(height(good_spiking), 1);


        for i = 1:length(current_cells.CellID)
            current_cellID = current_cells.CellID(i);
            cell_region_idx = current_hist.CellID == current_cellID;
            cell_cluster_idx = current_cells.CellID == current_cellID;
            current_region = current_hist.brain_region(cell_region_idx);
            current_cluster = current_cells.ClusterID(cell_cluster_idx);

            current_good_spiking_idx = good_spiking.CellID == current_cellID;
            good_spiking.ClusterID(current_good_spiking_idx) = current_cluster(1);

            good_spiking.Brain_Region(current_good_spiking_idx) = string(current_region{1});
        end
        
        hist_and_spiking_matrices{mouse} = good_spiking;

    end
end

%%