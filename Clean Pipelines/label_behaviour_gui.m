function [labels, extracted_features] = label_behaviour_gui(video_paths, Xfit, b, R,...
    T, missing, behaviour, output_folder, base_name)
%LABEL_BEHAVIOUR_GUI: Opens a GUI to label frames of your behaviour of choice
%and extract the relevant features for training a model to predict this
%behaviour. 
% 
% To use:
%   -Use left arrow/right arrow scrub through the video
%   -Use L to label a frame as behaviour
%   -Use D to unlabel a frame as behaviour
%   -Use the 'Save Labels' button to save labels and their associated
%   features to your output folder
%
% Inputs:
%   -video_paths: a cell array containing paths the 4 camera angles to use
%   as a reference when labelling
%   -behaviour: a char array of the behaviour being labelled, e.g., 'Grooming'
%   -output_folder: path to the location where labels/features will be
%   saved
%   -base_name: char array of the mouse number, e.g., 'mouse5'
%   -All other inputs are obtained from the relevant SSM file
%   
%   To run the function with all inputs pre-formatted, use label_and_train_RFM.mat
%
% Outputs:
%   -labels: a binary array the length of the video in frames, with 1
%   corresponding to behaviour and 0 to not behaviour
%   -extracted_features: a frame x feature array containing values for
%   each feature for every frame

num_cams = length(video_paths);
obj = cell(num_cams, 1);
for i = 1:num_cams
    obj{i} = VideoReader(video_paths{i});
end

num_frames = obj{1}.NumFrames;

% Extract features for all frames
extracted_features = extract_behavioural_features(Xfit, b, R, T, missing, behaviour);

% Initialise labels
labels = zeros(num_frames, 1);

% GUI setup
fig = figure('Name', ['Label ', behaviour], 'NumberTitle', 'off', ...
             'KeyPressFcn', @key_press_callback, 'Units', 'normalized', ...
             'Position', [0.05, 0.05, 0.9, 0.85]);

% Video axes
video_positions = {
    [0.19, 0.54, 0.40, 0.42];  % Cam 1
    [0.47, 0.54, 0.40, 0.42];  % Cam 4
    [0.19, 0.12, 0.40, 0.42];  % Cam 6
    [0.47, 0.12, 0.40, 0.42];  % Cam 9
};

axes_handles = gobjects(num_cams,1);
image_handles = gobjects(num_cams,1);
for cam = 1:num_cams
    ax = axes(fig, 'Units', 'normalized', 'Position', video_positions{cam});
    axes_handles(cam) = ax;
    frame = read(obj{cam}, 1);
    image_handles(cam) = imshow(frame, 'Parent', ax);
    axis(ax,'off');
end

% Shared title
super_title = annotation(fig,'textbox',[0.33,0.96,0.4,0.04], ...
                         'String','Frame 1 - Label: 0', ...
                         'HorizontalAlignment','center','FontWeight','bold','EdgeColor','none');

% Slider
slider = uicontrol(fig, 'Style', 'slider', 'Min',1,'Max',num_frames, ...
                   'Value',1,'SliderStep',[1/(num_frames-1),0.1], ...
                   'Units','normalized','Position',[0.17,0.1,0.7,0.02], ...
                   'Callback', @slider_callback);

% Buttons
uicontrol(fig,'Style','pushbutton','String','Play', ...
          'Units','normalized','Position',[0.24,0.03,0.1,0.05], ...
          'Callback',@play_video);
uicontrol(fig,'Style','pushbutton','String','Pause', ...
          'Units','normalized','Position',[0.36,0.03,0.1,0.05], ...
          'Callback',@pause_video);
uicontrol(fig,'Style','pushbutton','String','Label Frame', ...
          'Units','normalized','Position',[0.48,0.03,0.1,0.05], ...
          'Callback',@label_frame);
uicontrol(fig,'Style','pushbutton','String','De-label Frame', ...
          'Units','normalized','Position',[0.60,0.03,0.1,0.05], ...
          'Callback',@delabel_frame);
uicontrol(fig,'Style','pushbutton','String','Save Labels', ...
          'Units','normalized','Position',[0.72,0.03,0.1,0.05], ...
          'Callback',@save_labels);

% --- Nested Functions ---
    function key_press_callback(~, event)
        frame_idx = round(slider.Value);
        key = lower(event.Key); % ensure case-insensitive

        switch key
            case 'rightarrow'
                if frame_idx < num_frames
                    slider.Value = frame_idx + 1;
                    update_frame();
                end
            case 'leftarrow'
                if frame_idx > 1
                    slider.Value = frame_idx - 1;
                    update_frame();
                end
            case 'l'
                labels(frame_idx) = 1;
                if frame_idx < num_frames
                    slider.Value = frame_idx + 1; % advance
                end
                disp(['Frame ', num2str(frame_idx), ' labelled as ', lower(behaviour)])
                update_frame();
            case 'd'
                labels(frame_idx) = 0;
                if frame_idx < num_frames
                    slider.Value = frame_idx + 1; % advance
                end
                disp(['Frame ', num2str(frame_idx), ' unlabelled as ', lower(behaviour)])
                update_frame();
        end
    end

    function label_frame(~,~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 1;
        if frame_idx < num_frames
            slider.Value = frame_idx + 1; % advance
        end
        update_frame();
    end

    function delabel_frame(~,~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 0;
        if frame_idx < num_frames
            slider.Value = frame_idx + 1; % advance
        end
        update_frame();
    end


    function save_labels(~,~)
        if ~exist(output_folder,'dir')
            mkdir(output_folder);
        end
        save_path = fullfile(output_folder, [base_name, '_', lower(behaviour), '_labels.mat']);
        save(save_path,'labels','extracted_features');
        disp(['Saved labels and features to ', save_path]);
    end

    function update_frame()
        f = round(slider.Value);
        for cam = 1:num_cams
            frame = read(obj{cam}, f);
            set(image_handles(cam), 'CData', frame);
        end
        super_title.String = sprintf('Frame %d - Label: %d', f, labels(f));
    end

    function slider_callback(~,~)
        frame_idx = round(slider.Value);     % snap to nearest integer
        slider.Value = frame_idx;            % force slider to that value
        update_frame();                      % refresh display
    end

    isPlaying = false; % shared flag for play/pause
    
    % Play button
    function play_video(~,~)
        isPlaying = true;
        frame_idx = round(slider.Value);
        while frame_idx <= num_frames && isPlaying
            slider.Value = frame_idx;
            update_frame();
            pause(1/15); % adjust to FPS
            frame_idx = frame_idx + 1;
    
            % Stop if figure is closed
            if ~isvalid(fig)
                break
            end
        end
    end
    
    % Pause button
    function pause_video(~,~)
        isPlaying = false; % this will stop the loop in play_video
    end

end
