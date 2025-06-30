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
