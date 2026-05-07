function[imgs, worldPoints, boardSize] = load_video_set_params(video_path)

%video_path = 'C:\Users\mqbphrs3\OneDrive - The University of Manchester\code_cam_cal\cam calibration 240226';
%video_path = 'C:\Users\mqbphrs3\OneDrive - The University of Manchester\code_cam_cal\test_250226';
%video_path = 'C:\Users\mqbphrs3\OneDrive - The University of Manchester\code_cam_cal\calib280226';

% Videos location
% imgs{c} is a cell array of images for camera c
% e.g., imgs{1}{k}, imgs{2}{k}, ..

% Parameters
squareSize = 20; % in millimeters (change to your checker size)
boardSize  = [7, 11]; % [numCols, numRows] internal corners
worldPoints = generateCheckerboardPoints(boardSize, squareSize);

%load images
f = dir([video_path '\*.avi']);
numCams = numel(f);
imgs = cell(1,numCams);
for n = 1:numCams
    obj = VideoReader([video_path '\' f(n).name]);
    numFrames = obj.NumFrames;
    for m = 1:numFrames
        imgs{n}{m} = read(obj, m);
    end
end