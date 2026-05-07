% REVIEW_LABELLING
% Display frames classified  as behaviour or classified as not behaviour
%
% INPUTS
% vid_path: full path of video chosen to review frames 
% labels_path: full path of *labels.mat file 
% view behaviour: set to true to view behaviour, set to false to view
% non-behaviour

function[] = review_labelling(vid_path, labels_path, view_behaviour)

load(labels_path, 'labels')

vid_obj = VideoReader(vid_path);

vid_frames = 1:vid_obj.NumFrames;

if view_behaviour == true
    labels = logical(labels');
elseif view_behaviour == false
    temp = logical(labels');
    labels = ~temp; % Invert the labels for viewing non-behaviour frames
end

frames_to_watch = vid_frames(labels);

for i = frames_to_watch
    frame = read(vid_obj, i);
    imshow(frame);
    pause(0.05)
end

end