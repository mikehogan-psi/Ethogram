%% 3D triangulation of 2D-DEEPLABCUT-data from multiple camara angles
%  - this script reconstructs the 3D positions of mouse body-marker positions 

%  MAIN INPUT:  1. deeplabcut data (2D coordinates) of markers recorded from 4 different camera angles
%               2. previously computet camera projection matrix [P] (specific for these camera angles)
%               !!! deeplabcut data files MUST be in the same carmera orderas the order as the [P] is organised !!! 
%                           (data from camera_files{n} stems from same camera used to calibrate P{n})

%  PROCESSING STEPS: 1. Load data - extract DEEPLAB cut data from csv. files and reshapes it into correct matlab format
%                    2. triangulate data - performs triangulation to reconstruct the 3D coordinates

%  OUTPUT: 3D-data of markers and which information about cameras recorded them
%          the final files that are saved will contain the following variables:
%           - camera_files: filepaths to the camera files that data for triangulation was extracted from
%           - x : deeplabcut_data - for each frame (1st dim), x- y-positions of each body marker (2nd dim), each camera (3rd dim)
%           - L : likelihoods (confidence values) for each each frame (1rst dim), each body marker (2nd dim), each camera (3rd dim) 
%           - TH: likelihood treshold
%           - X : triangulated 3D data - for each body marker (1st dim), x- y- z-positions (2nd dim), each frame (3rd dim)
%           - W : number of cameras which detected body marker (1st dim) in a frame (2nd dim)


%% General Setup

% define session 
 sesh = 'extinction';
  % sesh = 'renewal';

% % load camera projection matrix (P)
%   load('C:\Users\Abi Hogan\Documents\GitHub\3D_camera_calibration\p_matrices\Pcal_hogan_9cameras.mat', 'P');
load('C:\Users\G71044MH\OneDrive - The University of Manchester\Documents\GitHub\3D_camera_calibration\p_matrices\Pcal_hogan_9cameras.mat', 'P')
 
% define folder path to folder that contains dlc data from all cameras
  % dlc_folder_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\data_all_mice\Extinction\data_from_DLC_iteration_3\DLC_data\';
  dlc_folder_path = 'C:\Cohort 4 Temp Data Storage\Mouse3\Extinction\DLC Data\';

% define folder path to where triangulated data files shall be saved
  % save_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\data_all_mice\Extinction\data_from_DLC_iteration_3\triangulated_data\camera_8_included\';
  save_path = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort_4_06_05_25 (SC PAG Implanted Animals)\Extinction\Mouse 3\Triangulated Data\';  

  % define likelihood treshold (only markers which exeed this threshold will be used)
  TH = 0.9; % 

%% A: test for single mouse

% define mouse file name
  base_name = ['mouse3_' sesh '_p1'];
  

% extract filenames from different cameras
  camera_files_struct = dir(fullfile(dlc_folder_path, ['camera*' base_name '*']));

% load camera file paths into cell array
  camera_files = cell(1, length(camera_files_struct));

  for n = 1:length(camera_files_struct)
      camera_files{n} = [dlc_folder_path camera_files_struct(n).name];    
  end 

% run load_data function 
  [x,L] = load_data_og(camera_files); % extacts 2D coordinates and their likelihood values and loads them into x and L variables

% run triangualtion function 
  [X, W] = triangulate_simple_og(x, L, P, TH); 

% plot data in 3D coordinate system 
  plot3d_video(X,false)

%% B - for all mice

% get a list of all filenames in the folder
  file_list = dir([dlc_folder_path '/camera*_mouse*']);
  file_names = cell(1,length(file_list));

  for n = 1: length(file_list)
      file_names{n} =  file_list(n).name;
  end

% Initialize a cell array to hold base names (e.g. mouse4_extinction_p1)
base_names = cell(size(file_names));

% Loop through each filename and extract base name using regex
for i = 1:length(file_names)
    % Match pattern: mouseX_extinction_pY

    tokens = regexp(file_names{i}, ['mouse\d+_' sesh '_p\d+'] , 'match');
    if ~isempty(tokens)
        base_names{i} = tokens{1};
    end
end

% Remove double entries 
base_names = unique(base_names(~cellfun('isempty', base_names)));

% loop through each mouse file and triagnulate data
for n = 1 : length(base_names)

    % extract filenames from different cameras
      camera_files_struct = dir(fullfile(dlc_folder_path, ['camera*' base_names{n} '*'])); 
      % !!! make sure cameras are in the correct orders (matching order of P matrix) !!!
    
    % load camera file paths into cell array
      camera_files = cell(1, length(camera_files_struct)); % initiate cell array with length corresponding to number of cameras used

      for i = 1:length(camera_files_struct)
          camera_files{i} = [dlc_folder_path camera_files_struct(i).name];
          % -> each cell will contain file paths for CSV files recorded by one camera   
      end 
    
    % run load_data function 
      [x,L] = load_data_og(camera_files); % extacts 2D coordinates and their likelihood values and loads them into x and L variables
        
    % run triangualtion function 
      [X, W] = triangulate_simple_og(x, L, P, TH); 

    % save data  
      save([save_path base_names{n} '_3D_triangulated'],"W","X","x","L","TH","camera_files"); 

end

