%% Setup

% !!!Specify which session is to be analysed!!!
% Acquisition, Extinction, or Renewal
session = 'Renewal';

% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!! Define likelihood threshold (only DLC data which exceed this
% confidence value will be used for triangulation) !!!
TH = 0.9; 

% !!! Provide filepath for camera projection matrix for triangulation !!!
load('C:\Users\G71044MH\OneDrive - The University of Manchester\Documents\GitHub\3D_camera_calibration\p_matrices\Pcal_hogan_9cameras.mat', 'P')

% !!! Provide filepath for SSM for estimating data !!!
SSM_model_path = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Models\SSMs\SSM_3D_implant_mouse1_head_fixed.mat';
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

%% Triangulate DLC data and save in appropriate subfolder

% Master loop goes over each mouse and accesses its DLC data
for mouse = 1:length(dlc_session_folders)
    current_dlc_folder_path = dlc_session_folders{mouse};
    dlc_file_list = dir(fullfile(current_dlc_folder_path, 'camera*_mouse*_p*.csv'));
    dlc_file_list = dlc_file_list(~[dlc_file_list.isdir]);
    base_names = cell(length(dlc_file_list), 1);
    dlc_file_names = cell(length(dlc_file_list), 1);
    
    % Extract file names from list of DLC files
    for n = 1:length(dlc_file_list)
        dlc_file_names{n} =  dlc_file_list(n).name;
    end
    
    % Loop through each filename and extract base name using regexp
    for i = 1:length(dlc_file_names)
    % Match pattern: mouseX_session_pY
        tokens = regexp(dlc_file_names{i}, ['mouse\d+_' session '_p\d+'] , 'match',...
            'ignorecase');
        if ~isempty(tokens)
          base_names{i} = tokens{1};
        end
    end

    % Remove double entries 
    base_names = unique(base_names(~cellfun('isempty', base_names)));

    % Loop through session parts (1 and 2)
    for part = 1 : length(base_names)
    
        % Extract filenames for each camera for each part
        camera_files_struct = dir(fullfile(current_dlc_folder_path, ['camera*' base_names{part} '*.csv'])); 
        % Ensure cameras are in the correct order (matching order of P matrix)
        
        % Load camera file paths into cell array
        camera_files = cell(1, length(camera_files_struct)); % initiate cell array with length corresponding to number of cameras used
    
        for camera = 1:length(camera_files_struct)
            camera_files{camera} = fullfile(current_dlc_folder_path, camera_files_struct(camera).name);
            % -> each cell will contain file paths for CSV files recorded by one camera   
        end              

        
        % Get mouse names for file creation
        mouse_name = regexp(current_dlc_folder_path, 'Mouse\s*\d+', 'match', 'once');

        triangulated_data_folder = fullfile(master_directory, mouse_name,...
        session, 'Behavioural Data', 'Triangulated Data'); 
        % Skip if save directory cannot be found and display name of
        % directory that needs creating
        if ~exist(triangulated_data_folder, 'dir')
            warning(['Save directory: ' triangulated_data_folder ' does not exist, skipping this file'])
        continue
        end

        triangulated_data_save_path = fullfile(triangulated_data_folder,...
            [base_names{part} '_triangulated.mat']);
        % Check if triangulated data already exists, skip this file if it
        % does
        if exist(triangulated_data_save_path, 'file')
        warning([base_names{part} '_triangulated.mat already exists, skipping triangulation.']);
        else
        
        % Run load_data function 
        disp(['Loading DLC data for ' base_names{part} '...']);
        [x,L] = load_data_og(camera_files); % Extacts 2D coordinates and their likelihood values and loads them into x and L variables
        
        disp(['Triangulating data for ' base_names{part} '...'])
        % Run triangulation function 
        [X, W] = triangulate_simple_og(x, L, P, TH); 
        
        % Save data in appropriate subfolder and with appropriate name
        save(triangulated_data_save_path,"W","X","x","L","TH","camera_files"); 
        disp([base_names{part} '_triangulated.mat' ' saved to ' triangulated_data_folder])
        end
        
    end

end

% call 'plot3d_video(X, false)' if you want to check triangulation has worked
%% Fit SSM to triangulated data

% Initialise 
triangulated_session_folders = cell(length(mouse_files), 1);

% Get triangulated data folders for each mouse for specified session
for mouse = 1:length(mouse_files)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    triangulated_session_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Triangulated Data');

    current_triangulated_folder_path = triangulated_session_folders{mouse};
    
    triangulated_file_list = dir(fullfile(current_triangulated_folder_path, 'mouse*_triangulated.mat'));
    triangulated_file_list = triangulated_file_list(~[triangulated_file_list.isdir]);
    
    % Initialise
    triangulated_file_paths = cell(length(triangulated_file_list), 1);
    
    % Extract file paths from list of triangulated files
    for n = 1:length(triangulated_file_list)
        triangulated_file_paths{n} =  fullfile(triangulated_file_list(n).folder, triangulated_file_list(n).name);
    end
    
    % Loop through each filename and extract base name using regex
    base_names = cell(length(triangulated_file_paths), 1);
    for i = 1:length(triangulated_file_paths)
    % Match pattern: mouseX_session_pY
        tokens = regexp(triangulated_file_paths{i}, ['mouse\d+_' session '_p\d+'] , 'match',...
            'ignorecase');
        if ~isempty(tokens)
          base_names{i} = tokens{1};
        end
    end
    
    % Setting save folder
    SSM_data_save_folder = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'SSM Fitted Data' );
    if ~exist(SSM_data_save_folder, 'dir')
        warning(['Save directory: ' SSM_data_save_folder ' does not exist,' ...
            ' skipping this file' ])
        continue
    end
    
    % Fitting SSM for each part
    for part = 1:length(triangulated_file_paths)
        SSM_save_path = fullfile(SSM_data_save_folder, [base_names{part},...
            '_SSM_fitted.mat']);
        if exist(SSM_save_path, 'file')
            warning([base_names{part}, '_SSM_fitted.mat already exists,' ...
            ' skipping this file'])
            continue
        end
        disp(['Fitting SSM for ' base_names{part}, '...'])
        Fit_SSM_3D_new(triangulated_file_paths{part}, SSM_save_path, SSM_model_path)
    end

end

% call 'plot3d_video(Xfit, false)' if you want to check fitting has worked

%% Triangulate DLC data for habituation and save in appropriate subfolder

% Master loop goes over each mouse and accesses its DLC data
for mouse = 1:length(dlc_session_folders)
    current_dlc_folder_path = dlc_session_folders{mouse};
    dlc_file_list = dir(fullfile(current_dlc_folder_path, 'camera*_mouse*_habituation*.csv'));
    dlc_file_list = dlc_file_list(~[dlc_file_list.isdir]);
    base_names = cell(length(dlc_file_list), 1);
    dlc_file_names = cell(length(dlc_file_list), 1);
    
    % Extract file names from list of DLC files
    for n = 1:length(dlc_file_list)
        dlc_file_names{n} =  dlc_file_list(n).name;
    end
    
    % Loop through each filename and extract base name using regexp
    for i = 1:length(dlc_file_names)
    % Match pattern: mouseX_session_pY
        tokens = regexp(dlc_file_names{i}, ['mouse\d+_' session '_habituation'] , 'match',...
            'ignorecase');
        if ~isempty(tokens)
          base_names{i} = tokens{1};
        end
    end

  
    
    % Extract filenames for each camera for each part
    camera_files_struct = dir(fullfile(current_dlc_folder_path, ['camera*' base_names{1} '*.csv'])); 
    % Ensure cameras are in the correct order (matching order of P matrix)
    
    % Load camera file paths into cell array
    camera_files = cell(1, length(camera_files_struct)); % initiate cell array with length corresponding to number of cameras used

    for camera = 1:length(camera_files_struct)
        camera_files{camera} = fullfile(current_dlc_folder_path, camera_files_struct(camera).name);
        % -> each cell will contain file paths for CSV files recorded by one camera   
    end              

    
    % Get mouse names for file creation
    mouse_name = regexp(current_dlc_folder_path, 'Mouse\s*\d+', 'match', 'once');

    triangulated_data_folder = fullfile(master_directory, mouse_name,...
    session, 'Behavioural Data', 'Triangulated Data'); 
    % Skip if save directory cannot be found and display name of
    % directory that needs creating
    if ~exist(triangulated_data_folder, 'dir')
        warning(['Save directory: ' triangulated_data_folder ' does not exist, skipping this file'])
    continue
    end

    triangulated_data_save_path = fullfile(triangulated_data_folder,...
        [base_names{1} '_triangulated.mat']);
    % Check if triangulated data already exists, skip this file if it
    % does
    if exist(triangulated_data_save_path, 'file')
    warning([base_names{1} '_triangulated.mat already exists, skipping triangulation.']);
    else
    
    % Run load_data function 
    disp(['Loading DLC data for ' base_names{1} '...']);
    [x,L] = load_data_og(camera_files); % Extacts 2D coordinates and their likelihood values and loads them into x and L variables
    
    disp(['Triangulating data for ' base_names{1} '...'])
    % Run triangulation function 
    [X, W] = triangulate_simple_og(x, L, P, TH); 
    
    % Save data in appropriate subfolder and with appropriate name
    save(triangulated_data_save_path,"W","X","x","L","TH","camera_files"); 
    disp([base_names{1} '_triangulated.mat' ' saved to ' triangulated_data_folder])
    end
        
end

%% Fit SSM to triangulated data for habituation

% Initialise 
triangulated_session_folders = cell(length(mouse_files), 1);

% Get triangulated data folders for each mouse for specified session
for mouse = 1:length(mouse_files)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    triangulated_session_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Triangulated Data');

    current_triangulated_folder_path = triangulated_session_folders{mouse};
    
    triangulated_file_list = dir(fullfile(current_triangulated_folder_path, 'mouse*_habituation_triangulated.mat'));
    triangulated_file_list = triangulated_file_list(~[triangulated_file_list.isdir]);
    
    % Initialise
    triangulated_file_paths = cell(length(triangulated_file_list), 1);
    
    % Extract file paths from list of triangulated files
    for n = 1:length(triangulated_file_list)
        triangulated_file_paths{n} =  fullfile(triangulated_file_list(n).folder, triangulated_file_list(n).name);
    end
    
    % Loop through each filename and extract base name using regex
    base_names = cell(length(triangulated_file_paths), 1);
    for i = 1:length(triangulated_file_paths)
    % Match pattern: mouseX_session_habituation
        tokens = regexp(triangulated_file_paths{i}, ['mouse\d+_' session '_habituation'] , 'match',...
            'ignorecase');
        if ~isempty(tokens)
          base_names{i} = tokens{1};
        end
    end
    
    % Setting save folder
    SSM_data_save_folder = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'SSM Fitted Data' );
    if ~exist(SSM_data_save_folder, 'dir')
        warning(['Save directory: ' SSM_data_save_folder ' does not exist,' ...
            ' skipping this file' ])
        continue
    end
    
    % Fitting SSM 
    SSM_save_path = fullfile(SSM_data_save_folder, [base_names{1},...
        '_SSM_fitted.mat']);
    % if exist(SSM_save_path, 'file')
    %     warning([base_names{1}, '_SSM_fitted.mat already exists,' ...
    %     ' skipping this file'])
    %     continue
    % end
    disp(['Fitting SSM for ' base_names{1}, '...'])
    Fit_SSM_3D_habituation(triangulated_file_paths{1}, SSM_save_path, SSM_model_path)


end