%% For extracting names from 'video files'
% Define the base directory where your videos are stored
baseDir = 'D:\PhD 2nd Year\DeepLabCut Models\Fear Conditioning With Implant 9 Cam-Mike-2025-05-16\videos';

% Get all files recursively
allFiles = dir(fullfile(baseDir, '**', '*.*'));

% Keep only .avi files (case-insensitive)
aviFiles = allFiles(endsWith(lower({allFiles.name}), '.avi'));

% Crop parameters
cropStr = ': crop: 0, 1280, 0, 1024';

% Loop and print full paths formatted for YAML
for i = 1:length(aviFiles)
    fullPath = fullfile(aviFiles(i).folder, aviFiles(i).name);
    fprintf('  ? %s\n  %s\n', fullPath, cropStr);
end

%% For extracting names from 'labeled data'
% Path to your DLC project
projectPath = 'D:\PhD 2nd Year\DeepLabCut Models\Fear Conditioning With Implant 9 Cam-Mike-2025-05-16';
labeledDataPath = fullfile(projectPath, 'labeled-data');

% Get all subfolders inside labeled-data
d = dir(labeledDataPath);
subfolders = {d([d.isdir] & ~ismember({d.name},{'.','..'})).name};

% Crop parameters
cropStr = ': crop: 0, 1280, 0, 1024';  % adjust if needed

% Print header
fprintf('video_sets:\n');

% Loop and print YAML-style block
for i = 1:numel(subfolders)
    videoName = subfolders{i};
    
    % Dummy/fake path with correct extension (adjust extension if your videos were .avi)
    fakePath = fullfile('C:\dummy\', [videoName '.avi']);
    
    % Print in your preferred style
    fprintf('  ? %s\n  %s\n', fakePath, cropStr);
end
