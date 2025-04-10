%% Load in freezing data from tracked animals
load() %add path to validated_freeze_matrix
num_frames = 502;
num_trials = 40;
num_mice = 4; % Change as needed

%% Load in eye tracking data

% Specify the folder containing your .csv files
folderPathEyeTracking = 'D:\simulated data';

% Get a list of all .csv files in the specified folder
fileListEyeTracking = dir(fullfile(folderPathEyeTracking, '*.csv'));

% Initialise a cell array to store each T matrix
eye_data = cell(1, length(fileListEyeTracking));
file_names = (cell(length(fileListEyeTracking), 1));

% iterate over each file, load it, and extract eye landmarks
for i = 1:length(fileListEyeTracking)
    filePathEyeTracking = fullfile(folderPathEyeTracking, fileListEyeTracking(i).name);
     % load the .csv file
    data = readmatrix(filePathEyeTracking); 
    % remove unecessary columns
    data(:, 1) = [];
    likelihood_ind = 3:3:size(data, 2);
    data(:, likelihood_ind) = [];
    eye_data{i} = data';
    % stores body_anterior x and y positional data in T_matrices array
    file_names_eye_tracking{i} = fileListEyeTracking(i).name;
    % stores file names for sorting later
end
%% Sorting eye tracking data into mouse/part order

mouse_numbers_eye_tracking = cellfun(@(x) sscanf(x, 'mouse%d'), file_names_eye_tracking);
% extracting the mouse number from the file names

[~, sort_idx] = sort(mouse_numbers_eye_tracking);
% sort by mouse number amd get index for sorted order

eye_data = eye_data(sort_idx);
% sort matrix
disp('Sorted Eye Tracking File Names:');
disp(file_names_eye_tracking(sort_idx)');
%% Putting parts of eyetracking videos together for each mouse

eye_data_paired = cell(1, length(eye_data)/2);
%initialising new cell array half the size of original

for video_ix = 1:2:length(eye_data)
    %iterating over pairs of cells
    eye_data_paired{(video_ix+1)/2} = cat(2, eye_data{video_ix}, eye_data{video_ix+1});
    %stores concatenated matrices in correct indices of new matrix
end

eye_data = eye_data_paired;

%% Splitting into pupil points and eye edge points

pupil_points = cell(size(eye_data));

% Load pupil points and sort into xy
for mouse_idx = 1:num_mice
    
    pupil_raw = eye_data{mouse_idx}(13:20, :);
    
    nFrames = size(pupil_raw, 2);
    
    x_cords = pupil_raw(1:2:end, :);
    y_cords = pupil_raw(2:2:end, :);

    temp = zeros(nFrames, 4, 2);
    temp(:, : , 1) = x_cords';
    temp(:, :, 2) = y_cords';

    pupil_points{mouse_idx} = temp;  

end

% Load eye edge points and sort into xy
eye_edge_points = cell(size(eye_data));

for mouse_idx = 1:num_mice
    
    edge_raw = eye_data{mouse_idx}(1:12, :);
    
    nFrames = size(edge_raw, 2);
    
    x_cords = edge_raw(1:2:end, :);
    y_cords = edge_raw(2:2:end, :);

    temp = zeros(nFrames, 6, 2);
    temp(:, : , 1) = x_cords';
    temp(:, :, 2) = y_cords';

    eye_edge_points{mouse_idx} = temp;
    
end
%%
eye_areas = cell(num_mice, 1);
pupil_areas = cell(num_mice, 1);
eye_area_max = cell(num_mice, 1);
pupil_to_eye_ratio = cell(num_mice, 1);

for mouse_idx = 1:num_mice
   
    eye = eye_edge_points{mouse_idx};
    pupil = pupil_points{mouse_idx};

    nFrames = size(eye, 1);
    eye_area = zeros(nFrames, 1);
    pupil_area = zeros(nFrames, 1);

    for f = 1:nFrames
        ex = eye(f, :, 1);
        ey = eye(f, :, 2);
        px = pupil(f, :, 1);
        py = pupil(f, :, 2);

        eye_area(f) = 0.5 * abs( sum(ex .* circshift(ey, -1) - ey .* circshift(ex, -1)) );
        pupil_area(f) = 0.5 * abs( sum(px .* circshift(py, -1) - py .* circshift(px, -1)) );
    end
    % Extract pupil area, eye area and pupil size (relative to maximum eye
    % size)
    eye_areas{mouse_idx} = eye_area;
    eye_area_max{mouse_idx} = max(eye_area);
    pupil_areas{mouse_idx} = pupil_area;
    pupil_to_eye_ratio{mouse_idx} = pupil_area./eye_area_max{mouse_idx};

end

%% Index eye data using freezing data to remove movement as a variable in pupil size

pupil_to_eye_ratio_double = zeros(num_trials,num_frames,num_mice);


% Converting eye data from cell array to double for compatability with
% freezing data
for mouse_idx = 1:num_mice
    pupil_to_eye_ratio_double(:, :, mouse_idx) = reshape(pupil_to_eye_ratio{mouse_idx}, num_trials, num_frames);
end

pupil_ratio_no_movement =  zeros(num_trials, num_frames, num_mice);

for mouse_idx = 1:num_mice
    pupil_ratio = pupil_to_eye_ratio_double(:, : , mouse_idx);
    freeze_idx = squeeze(logical(validated_freeze_matrix(:, :, mouse_idx)));
    pupil_ratio(~freeze_idx) = NaN; % set non-freezing indices to NaN
    pupil_ratio_no_movement(:, :, mouse_idx) = pupil_ratio;
end

%% Plot change in pupil size over time