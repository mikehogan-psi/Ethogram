
predicted_labels_path = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Extracted Behaviours\Darting\mouse3_extinction_darting.mat";
video_file_path = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 3\Extinction\Behavioural Data\Video Data\camera6_mouse3_extinction_p1_2025-07-11-132427-0000.avi";


% Load predicted labels
data = load(predicted_labels_path);  
predicted_labels = data.predicted_labels;
predicted_labels = predicted_labels(1:end/2, :);

% Open video
vidObj = VideoReader(video_file_path);

% Find frames predicted as the behaviour
frames_to_check = find(predicted_labels == 1);  % change to 0 to check non-behaviour

% Loop through selected frames
for i = 1:length(frames_to_check)
    frame_idx = frames_to_check(i);
    
    % Read and display frame
    frame = read(vidObj, frame_idx);
    imshow(frame);
    title(sprintf('Frame %d - Predicted: %d', frame_idx, predicted_labels(frame_idx)));
    
    % Pause to inspect (adjust duration or press a key to continue)
    pause(0.05);  % pause 0.1 seconds per frame
end
