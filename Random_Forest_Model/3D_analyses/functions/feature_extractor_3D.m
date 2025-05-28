function feature_extractor_3D(video_path, mat_file_path, output_folder, frame_features_strings, window_size)

    % Load the 4 videos
    num_cams = 4;
    obj = cell(num_cams, 1);
    for i = 1:num_cams
        obj{i} = VideoReader(video_path{i});
    end

    % Check frame count consistency
    num_frames = obj{1}.NumFrames;
    for i = 2:num_cams
        if obj{i}.NumFrames ~= num_frames
            error('All videos must have the same number of frames');
        end
    end

    % Load pose data
    data = load(mat_file_path);
    b = data.b;
    R = rotm2eul(data.R, 'ZYX')';
    T = data.T;
    missing = data.missing;
    X = data.X;

    % Initialise
    features = [];
    labels = zeros(num_frames, 1);
    all_frames_labels = zeros(num_frames, 1);
    playing = false;

    % Create GUI with tiled layout for multiple cameras
   fig = figure('Name', 'Video Player with Labelling', ...
             'NumberTitle', 'off', ...
             'KeyPressFcn', @key_press_callback, ...
             'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.85]);

    % Create a panel for video display (top 85% of window)
    % Create GUI window
    fig = figure('Name', 'Video Player with Labelling', ...
                 'NumberTitle', 'off', ...
                 'KeyPressFcn', @key_press_callback, ...
                 'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.85]);
    
    % Create a panel for video display (top 85% of window)
    video_panel = uipanel(fig, 'Units', 'normalized', ...
                          'Position', [0, 0.15, 1, 0.85]);
    
    % Tiled layout inside video panel
    t = tiledlayout(video_panel, 2, 2, ...
                    'Padding', 'none', 'TileSpacing', 'none');  % Maximize video area
    
    axes_handles = gobjects(num_cams, 1);
    image_handles = gobjects(num_cams, 1);
    
    for cam = 1:num_cams
        ax = nexttile(t);
        axes_handles(cam) = ax;
        frame = read(obj{cam}, 1);
        if cam <= 2
            frame = flipud(frame);  % Flip upside down
        end
        image_handles(cam) = imshow(frame, 'Parent', ax);
        axis(ax, 'off');  % Turn off axis lines and ticks
    end

    % Shared title with label info
    super_title = sgtitle(sprintf('Frame 1 - Label: %d', labels(1)));

    % Slider
    slider = uicontrol(fig, 'Style', 'slider', 'Min', 1, 'Max', num_frames, ...
                       'Value', 1, 'SliderStep', [1/(num_frames-1), 0.1], ...
                       'Units', 'normalized', 'Position', [0.17, 0.1, 0.7, 0.02], ...
                       'Callback', @slider_callback);

    % Buttons
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

    % Display first frame
    update_frame(1);

    % Slider callback
    function slider_callback(~, ~)
        update_frame(round(slider.Value));
    end

    % Frame update
    function update_frame(frame_idx)
    if frame_idx >= 1 && frame_idx <= num_frames
        for cam = 1:num_cams
            frame = read(obj{cam}, frame_idx);
            if cam <= 2
                frame = flipud(frame);  % Flip camera 1 and 2 vertically
            end
            image_handles(cam).CData = frame;
            %title(axes_handles(cam), sprintf('Camera %d - Frame %d', cam, frame_idx));
        end
        super_title.String = sprintf('Frame %d - Label: %d', frame_idx, labels(frame_idx));
        slider.Value = frame_idx;
    end
 end


    % Playback controls
    function play_video(~, ~)
        playing = true;
        while playing && slider.Value < num_frames
            slider.Value = slider.Value + 1;
            update_frame(round(slider.Value));
            pause(1 / obj{1}.FrameRate);
        end
    end

    function pause_video(~, ~)
        playing = false;
    end

    % Labelling functions
    function label_frame(~, ~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 1;
        all_frames_labels(frame_idx) = 1;
        disp(['Frame ', num2str(frame_idx), ' labelled as behaviour.']);
        frame_features = calculate_features(frame_idx);
        features = [features; frame_features];
        if frame_idx < num_frames
            slider.Value = frame_idx + 1;
            update_frame(round(slider.Value));
        end
    end

    function delabel_frame(~, ~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 0;
        all_frames_labels(frame_idx) = 0;
        disp(['Frame ', num2str(frame_idx), ' de-labelled.']);
        update_frame(frame_idx);
    end

    function save_labels(~, ~)
        balance_dataset();
        if ~endsWith(output_folder, filesep)
            output_folder = [output_folder, filesep];
        end
        [~, base_name, ~] = fileparts(video_path{1});
        save_path = [output_folder, base_name, '_labels.mat'];
        try
            save(save_path, 'labels', 'features', 'all_frames_labels');
            disp(['Saved to ', save_path]);
        catch ME
            disp(['Error saving: ', ME.message]);
        end
    end

    % Feature extraction for single frame (same as before)
    function frame_features = calculate_features(frame_idx)
        start_idx = max(1, frame_idx - (window_size / 2));
        end_idx = min(num_frames, frame_idx + (window_size / 2));
        b_window = b(:, start_idx:end_idx);
        T_window = T(:, start_idx:end_idx);
        R_window = R(:, start_idx:end_idx);
        miss_window = missing(start_idx:end_idx);
        miss_num = sum(miss_window);
        if miss_num >= 1 && miss_num <= 6
            valid_cols = ~any(isnan([b_window; T_window; R_window]), 1);
            b_window = b_window(:, valid_cols);
            T_window = T_window(:, valid_cols);
            R_window = R_window(:, valid_cols);
        end
        frame_features = zeros(length(frame_features_strings), 1);
        for i = 1:length(frame_features_strings)
            frame_features(i) = eval(frame_features_strings{i});
        end
    end

    function balance_dataset()
        behaviour_indices = find(labels == 1);
        non_behaviour_indices = find(labels == 0);
        num_behaviour = numel(behaviour_indices);
        num_non_behaviour_to_sample = min(5 * num_behaviour, numel(non_behaviour_indices));
        selected_non_behaviour_indices = datasample(non_behaviour_indices, num_non_behaviour_to_sample, 'Replace', false);
        balanced_indices = [behaviour_indices; selected_non_behaviour_indices];
        balanced_features = arrayfun(@calculate_features, balanced_indices, 'UniformOutput', false);
        features = cell2mat(balanced_features');
        features = features';
        labels = labels(balanced_indices);
        valid_rows = all(~isnan(features), 2);
        features = features(valid_rows, :);
        labels = labels(valid_rows);
    end

    function key_press_callback(~, event)
        if strcmp(event.Key, 'rightarrow') && slider.Value < num_frames
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
end
