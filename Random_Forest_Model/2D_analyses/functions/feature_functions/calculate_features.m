 % Feature extraction function
    function frame_features = calculate_features(frame_idx)
        num_frames = 10040;
        window_size = 30;  % Define window size (/fps to get time over which features are being computed)
        start_idx = max(1, frame_idx - floor(window_size / 2));
        end_idx = min(num_frames, frame_idx + floor(window_size / 2));

        % Extract data for the window
        b_window = b(:, start_idx:end_idx);
        T_window = T(:, start_idx:end_idx);
        A_window = A(start_idx:end_idx);

        % Calculate velocity and angular velocity, these can be changed to
        % capture the desired behaviour      
        angular_velocity = [0, abs(diff(A_window))];  % Prepend 0 for alignment
        velocity = [0, sqrt(sum(diff(T_window, 1, 2).^2))];  % Euclidean velocity

        % Compute features for the current window
        frame_features = [
            mean(b_window(3, :));  % Mean of eigenpose 3 (looking up/down)
            var(b_window(3, :));   % Variance of eigenpose 3
            mean(diff(b_window(3, :)));  % Temporal derivative of eigenpose 3
            mean(b_window(2, :));  % Mean of eigenpose 2 (elongation/hunching)
            mean(velocity);        % Mean velocity
            mean(angular_velocity) % Mean angular velocity
        ];
     end