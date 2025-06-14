function validate_predictions(base_name, video_file_path, predicted_labels_path, output_folder)
    % Load the video
    obj = VideoReader(video_file_path);

    % Load predicted labels
    data = load(predicted_labels_path);
    predicted_labels = data.predicted_labels;

    % Ensure the number of labels matches the number of frames
    num_frames = obj.NumFrames;
    if length(predicted_labels) ~= num_frames
        error('The number of predicted labels does not match the number of video frames.');
    end

    % Initialize verified labels with predicted labels
    verified_labels = predicted_labels;

    % Create GUI
    fig = figure('Name', 'Video Player with Prediction Validation', 'NumberTitle', 'off', ...
                 'KeyPressFcn', @key_press_callback);
    ax = axes(fig, 'Position', [0.1, 0.2, 0.8, 0.7]);
    im = imshow(read(obj, 1), 'Parent', ax);
    title_text = title(ax, sprintf('Frame 1 - Prediction: %d, Verified: %d', predicted_labels(1), verified_labels(1)));

    % Create slider
    slider = uicontrol(fig, 'Style', 'slider', 'Min', 1, 'Max', num_frames, ...
                       'Value', 1, 'SliderStep', [1/(num_frames-1), 0.1], ...
                       'Units', 'normalized', ...
                       'Position', [0.17, 0.1, 0.7, 0.02], ...
                       'Callback', @slider_callback);

    % Create buttons
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play', ...
              'Units', 'normalized', 'Position', [0.27, 0.03, 0.1, 0.05], ...
              'Callback', @play_video);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Pause', ...
              'Units', 'normalized', 'Position', [0.39, 0.03, 0.1, 0.05], ...
              'Callback', @pause_video);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Correct Label', ...
              'Units', 'normalized', 'Position', [0.51, 0.03, 0.1, 0.05], ...
              'Callback', @correct_label);
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save Labels', ...
              'Units', 'normalized', 'Position', [0.63, 0.03, 0.1, 0.05], ...
              'Callback', @save_labels);

    % Playback control
    playing = false;

    % --- Glider bar for showing labeled frames ---
    glider_bar_ax = axes(fig, 'Position', [0.17, 0.13, 0.7, 0.01]);
    glider_bar_patches = [];  % Will hold patch handles
    set(glider_bar_ax, 'XLim', [1, num_frames], ...
                       'YLim', [0, 1], ...
                       'XTick', [], 'YTick', [], ...
                       'Color', 'w', 'XColor', 'none', 'YColor', 'none', ...
                       'Box', 'on');
    update_glider_bar();  % Initialize glider bar

    % Slider callback
    function slider_callback(~, ~)
        frame_idx = round(slider.Value);
        update_frame(frame_idx);
    end

    % Update displayed frame
    function update_frame(frame_idx)
        if frame_idx >= 1 && frame_idx <= num_frames
            frame = read(obj, frame_idx);
            im.CData = frame;

            % Update title
            set(title_text, 'String', ...
                sprintf('Frame %d - Prediction: %d, Verified: %d', ...
                        frame_idx, predicted_labels(frame_idx), verified_labels(frame_idx)));

            % Update background color
            if verified_labels(frame_idx) == 1
                set(fig, 'Color', [1, 0.8, 0.8]); % Light red background
            else
                set(fig, 'Color', [1, 1, 1]); % White background
            end

            slider.Value = frame_idx;
        end
    end

    % Play video
    function play_video(~, ~)
        playing = true;
        while playing && slider.Value < num_frames
            slider.Value = slider.Value + 1;
            update_frame(round(slider.Value));
            pause(1 / obj.FrameRate);
        end
    end

    % Pause video
    function pause_video(~, ~)
        playing = false;
    end

    % Correct label for current frame
    function correct_label(~, ~)
        frame_idx = round(slider.Value);
        verified_labels(frame_idx) = 1 - verified_labels(frame_idx);  % Toggle label
        fprintf('Frame %d label corrected to %d.\n', frame_idx, verified_labels(frame_idx));
        update_frame(frame_idx);
        update_glider_bar();
    end

    % Save verified labels to file
    function save_labels(~, ~)
        if ~endsWith(output_folder, filesep)
            output_folder = [output_folder, filesep];
        end
        save_path = [output_folder, base_name, '_verified_labels.mat'];

        try
            save(save_path, 'verified_labels');
            disp(['Verified labels saved to ', save_path]);
        catch ME
            disp(['Error saving labels: ', ME.message]);
        end
    end

    % Handle key press shortcuts
    function key_press_callback(~, event)
        if strcmp(event.Key, 'rightarrow') && slider.Value < num_frames
            slider.Value = slider.Value + 1;
            update_frame(round(slider.Value));
        elseif strcmp(event.Key, 'leftarrow') && slider.Value > 1
            slider.Value = slider.Value - 1;
            update_frame(round(slider.Value));
        elseif strcmp(event.Key, 'c')
            correct_label();
        end
    end

    % Glider bar update function
    function update_glider_bar()
        if ~isempty(glider_bar_patches)
            delete(glider_bar_patches(ishandle(glider_bar_patches)));
        end

        labeled_frames = find(verified_labels == 1);
        glider_bar_patches = gobjects(size(labeled_frames));

        for i = 1:length(labeled_frames)
            x = labeled_frames(i);
            glider_bar_patches(i) = patch(glider_bar_ax, ...
                [x-0.5, x+0.5, x+0.5, x-0.5], ...
                [0, 0, 1, 1], ...
                'r', 'EdgeColor', 'none');
        end
    end

    % Display the first frame
    update_frame(1);
end
