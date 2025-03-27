function from_preselect_feature_extractor(video_file_path, preselected_file_path , SSM_file_path, manual_labels_path, frame_features_strings, window_size)
    % Load the video
    obj = VideoReader(video_file_path);

    % Load data
    data = load(SSM_file_path);
    b = data.b;
    A = data.A;
    T = data.T;
    missing = data.missing; 

    % Load preselected_labels data
    preselected_frames = load(preselected_file_path, "prior_freezing_selection");
    preselected_frames = preselected_frames.prior_freezing_selection;

    % Get indices of preselected frames
    selected_indices = find(preselected_frames);
    num_selected_frames = numel(selected_indices);

    % Initialise feature matrix 
    features = [];  % Store features for all selected frames as rows

    % Initialise labels
    labels = zeros(num_selected_frames, 1);  % Start with all frames as non-behaviour (0)
    all_frames_labels = zeros(num_selected_frames, 1); % Track labels before balancing

    % Create GUI
    fig = figure('Name', 'Video Player with Labelling', 'NumberTitle', 'off', ...
                 'KeyPressFcn', @key_press_callback);
    ax = axes(fig, 'Position', [0.1, 0.2, 0.8, 0.7]);
    im = imshow(read(obj, selected_indices(1)), 'Parent', ax);
    title(ax, sprintf('Frame %d - Label: %d', selected_indices(1), labels(1)));

    % Create slider
    slider = uicontrol(fig, 'Style', 'slider', 'Min', 1, 'Max', num_selected_frames, ...
                       'Value', 1, 'SliderStep', [1/(num_selected_frames-1), 0.1], ...
                       'Units', 'normalized', ...
                       'Position', [0.17, 0.1, 0.7, 0.02], ...
                       'Callback', @slider_callback);

    % Create buttons
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play', ...
              'Units', 'normalized', 'Position', [0.23, 0.03, 0.1, 0.05], ...
              'Callback', @play_video);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Pause', ...
              'Units', 'normalized', 'Position', [0.35, 0.03, 0.1, 0.05], ...
              'Callback', @pause_video);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Label Frame', ...
              'Units', 'normalized', 'Position', [0.47, 0.03, 0.1, 0.05], ...
              'Callback', @label_frame);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'De-label Frame', ...
              'Units', 'normalized', 'Position', [0.59, 0.03, 0.1, 0.05], ...
              'Callback', @delabel_frame);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save Labels', ...
              'Units', 'normalized', 'Position', [0.71, 0.03, 0.1, 0.05], ...
              'Callback', @save_labels);

    % Playback control
    playing = false;

    % Slider callback
    function slider_callback(~, ~)
        frame_idx = round(slider.Value);
        update_frame(frame_idx);
    end

    % Update displayed frame
    function update_frame(frame_idx)
        if frame_idx >= 1 && frame_idx <= num_selected_frames
            frame = read(obj, selected_indices(frame_idx));
            im.CData = frame;
            current_label = labels(frame_idx);
            title(ax, sprintf('Frame %d - Label: %d', selected_indices(frame_idx), current_label));
            slider.Value = frame_idx;
        end
    end

    % Play video
    function play_video(~, ~)
        playing = true;
        while playing && slider.Value < num_selected_frames
            slider.Value = slider.Value + 1;
            update_frame(round(slider.Value));
            pause(1 / obj.FrameRate);
        end
    end

    % Pause video
    function pause_video(~, ~)
        playing = false;
    end

    % Label current frame as behavior (1)
    function label_frame(~, ~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 1;
        all_frames_labels(frame_idx) = 1;
        disp(['Frame ', num2str(selected_indices(frame_idx)), ' labelled as behaviour.']);
    end

    % De-label current frame (set to 0)
    function delabel_frame(~, ~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 0;
        all_frames_labels(frame_idx) = 0;
        disp(['Frame ', num2str(selected_indices(frame_idx)), ' de-labelled.']);
        update_frame(frame_idx);
    end

    % Save labels
    function save_labels(~, ~)
        [~, video_name, ~] = fileparts(video_file_path);
        save_path = fullfile(manual_labels_path, [video_name, '_labels.mat']);
        save(save_path, 'labels', 'all_frames_labels');
        disp(['Labels saved to ', save_path]);
    end

    % Keyboard shortcut callback
    function key_press_callback(~, event)
        if strcmp(event.Key, 'rightarrow') && slider.Value < num_selected_frames
            slider.Value = slider.Value + 1;
            update_frame(round(slider.Value));
        elseif strcmp(event.Key, 'leftarrow') && slider.Value > 1
            slider.Value = slider.Value - 1;
            update_frame(round(slider.Value));
        elseif strcmp(event.Key, 'l')
            label_frame();
        elseif strcmp(event.Key, 'd')
            delabel_frame();
        end
    end

    % Display the first frame
    update_frame(1);
end
