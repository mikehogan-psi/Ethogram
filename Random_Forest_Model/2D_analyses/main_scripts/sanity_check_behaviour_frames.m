%% Show all behaviour labelled frames (sanity check)

% Define folders
labelFolder = 'Z:\Abi\behavioral_analysis\2D_behavioural analysis\habituation_analysis\darting_analyses\RFM_darting_17\predicted_labels_data';
videoFolder = 'Z:\Abi\behavioral_analysis\2D_behavioural analysis\data_all_mice\video_data\extinction_habituation';

% Get list of all video files
videoFiles = dir(fullfile(videoFolder, '*.avi'));

for k = 1:length(videoFiles)
    videoName = videoFiles(k).name;
    baseNameShort = erase(videoName, '.avi');
    labelFileName = ['predicted_darting_labels_' baseNameShort '.mat'];
    
    % Load predicted labels
    labelFile = fullfile(labelFolder, labelFileName);
    if ~isfile(labelFile)
        warning('Label file not found for %s', videoName);
        continue;
    end
    
    data = load(labelFile);
    if ~isfield(data, 'predicted_labels_TH')
        warning('Variable "predicted_labels_TH" not found in %s', labelFile);
        continue;
    end
    predicted_labels_TH = data.predicted_labels_TH;
    
    % Get darting frame indices
    dartingIndices = find(predicted_labels_TH == 1);
    if isempty(dartingIndices)
        fprintf('No darting frames in %s\n', videoName);
        continue;
    end

    % Load video
    videoPath = fullfile(videoFolder, videoName);
    vid = VideoReader(videoPath);
    
    % Loop through only darting frame indices
    for i = 1:length(dartingIndices)
        frameIdx = dartingIndices(i);
        frameTime = (frameIdx - 1) / vid.FrameRate;
        vid.CurrentTime = frameTime;
        
        if hasFrame(vid)
            frame = readFrame(vid);
            imshow(frame);
            title(sprintf('%s - Darting Frame %d', baseNameShort, frameIdx));
            pause(0.01);  % minimal delay
        end
    end
end
