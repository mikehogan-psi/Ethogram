function feature_extractor_3D(video_path, mat_file_path, output_folder, frame_features_strings, window_size)
    % Load the video
    obj = VideoReader(video_path);

    % Load data
    data = load(mat_file_path);
    b = data.b;
    R = data.R;
    T = data.T;
    missing = data.missing; 

    % calculate Euler angles from R (rotational matrix) -> [yaw, pitch, roll] / [yaw, roll, pitch] for each frame
    R = rotm2eul(R, 'ZYX');
    R = R';

    % Initialise feature matrix 
    num_frames = obj.NumFrames;
    features = [];  % Store features for all frames as rows

    % Initialise labels
    labels = zeros(num_frames, 1);  % Start with all frames as non-behaviour (0)
    all_frames_labels = zeros(num_frames, 1); % Track labels before balancing

    % Create GUI
    fig = figure('Name', 'Video Player with Labelling', 'NumberTitle', 'off', ...
                 'KeyPressFcn', @key_press_callback);
    ax = axes(fig, 'Position', [0.1, 0.2, 0.8, 0.7]);
    im = imshow(read(obj, 1), 'Parent', ax);
    title_text = title(ax, sprintf('Frame 1 - Label: %d, Verified: %d', labels(1), 0));

    % Create slider
    slider = uicontrol(fig, 'Style', 'slider', 'Min', 1, 'Max', num_frames, ...
                       'Value', 1, 'SliderStep', [1/(num_frames-1), 0.1], ...
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
        if frame_idx >= 1 && frame_idx <= num_frames
            frame = read(obj, frame_idx);
            im.CData = frame;
            
            % Get the label for the current frame (default 0, or 1 if labeled as behavior)
            current_label = labels(frame_idx);
            
            % Update the title to display both frame number and its current label
            title(ax, sprintf('Frame %d - Label: %d', frame_idx, current_label));
            
            % Update the slider value
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

    % Label current frame as behavior (1)
    function label_frame(~, ~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 1;  % Mark as behaviour
        all_frames_labels(frame_idx) = 1; % Track original labeling
        disp(['Frame ', num2str(frame_idx), ' labelled as behaviour.']);

        % Compute and store features for the current frame
        frame_features = calculate_features(frame_idx);
        features = [features; frame_features];

        % Move to next frame
        if frame_idx < num_frames
            slider.Value = frame_idx + 1;
            update_frame(round(slider.Value));
        end
    end

    % De-label current frame (set to 0)
    function delabel_frame(~, ~)
        frame_idx = round(slider.Value);
        labels(frame_idx) = 0; % Remove behavior label
        all_frames_labels(frame_idx) = 0; % Track original labeling
        disp(['Frame ', num2str(frame_idx), ' de-labelled.']);
        
        % Update the displayed frame after de-labeling
        update_frame(frame_idx);
    end

    % Save labels and features to .mat file
    function save_labels(~, ~)
        % Balance dataset
        balance_dataset();

        % Ensure the output folder path ends with a separator
        if ~endsWith(output_folder, filesep)
            output_folder = [output_folder, filesep];
        end

        % Create the full path for saving the file
        [~, video_name, ~] = fileparts(video_path);
        save_path = [output_folder, video_name, '_labels.mat'];

        try
            save(save_path, 'labels', 'features', 'all_frames_labels');
            disp(['Labels, features, and all_frames_labels saved to ', save_path]);
        catch ME
            disp(['Error saving labels: ', ME.message]);
        end
    end

    % Feature extraction function 
    function frame_features = calculate_features(frame_idx)
        
        % Define window size that seems relevant to behaviour being analysed
        start_idx = max(1, frame_idx - (window_size / 2));
        end_idx = min(num_frames, frame_idx + (window_size / 2));

        % Extract data for the window
        b_window = b(:, start_idx:end_idx);
        T_window = T(:, start_idx:end_idx);
        R_window = R(:,start_idx:end_idx);

        % find how many frames are missing datapoints
        miss_window = missing(start_idx:end_idx);
        % Count number of frames qith missing datapoints
        miss_num = sum(miss_window);

        % If 1-6 frames with missing datapoints, exclude these  
        % (If 7 or more columns have NaNs, we proceed without excluding them)
        if miss_num >= 1 && miss_num <= 6
            valid_cols = ~any(isnan([b_window; T_window; R_window]), 1);
            b_window = b_window(:, valid_cols);
            T_window = T_window(:, valid_cols);
            R_window = R_window(:, valid_cols);
        end   
        

        % compute all features for this frame        
       

        frame_features = zeros(length(frame_features_strings), 1); % Preallocate for efficiency

        for i = 1:length(frame_features_strings)
            frame_features(i) = eval(frame_features_strings{i});
        end


    end
         
      % Balance dataset by adding more non-behaviour frames 
        function balance_dataset()
        behaviour_indices = find(labels == 1);
        non_behaviour_indices = find(labels == 0);

        % select non-behaviour frames (want eg. 5x more than behaviour frames)
        num_behaviour = numel(behaviour_indices);
        num_non_behaviour_to_sample = min(5 * num_behaviour, numel(non_behaviour_indices)); % to have a 5-to-1 non-behaviour/behaviour ratio
                                                                                            % OR if not enough non-behaviour frames to reach ratio -> just takes all non-behaviour frames
        selected_non_behaviour_indices = datasample(non_behaviour_indices, num_non_behaviour_to_sample, 'Replace', false);

        % Combine indices
        balanced_indices = [behaviour_indices; selected_non_behaviour_indices];
   
        % Extract features for the balanced dataset
        balanced_features = arrayfun(@calculate_features, balanced_indices, 'UniformOutput', false);
        features = cell2mat(balanced_features');
        features = features';
        labels = labels(balanced_indices);

        valid_rows = all(~isnan(features), 2);
        features = features(valid_rows, :);
        labels = labels(valid_rows);
    end

     % Keyboard shortcut callback
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

    % Display the first frame
    update_frame(1);

   
end
