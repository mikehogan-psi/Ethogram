%%
load('D:\PhD 2nd Year\T_matrices.mat');

num_frames = 502;
num_trials = 40;
num_mice = 30;

% Initialise storage for all positional data
all_positional_data = zeros(2, num_frames, num_trials, num_mice);

for mouse_idx = 1:num_mice
    % Extract x and y position data for the current mouse
    x_posi = T_matrices{mouse_idx}(1, :);
    y_posi = T_matrices{mouse_idx}(2, :);
    
    % Temporary storage for trials of the current mouse
    current_positional_data = zeros(2, num_frames, num_trials);
    
    for trial_idx = 1:num_trials
        % Extract data for the current trial
        start_idx = (trial_idx - 1) * num_frames + 1;
        end_idx = trial_idx * num_frames;
        
        % Assign positional data to temporary storage
        current_positional_data_x = x_posi(start_idx:end_idx);
        current_positional_data_y = y_posi(start_idx:end_idx);
        
        % Reshape and assign to the correct position in
        % current_positional_data
        current_positional_data(1, :, trial_idx) = current_positional_data_x; % x positions
        current_positional_data(2, :, trial_idx) = current_positional_data_y; % y positions
    end
    
    % Assign data for the current mouse to the main storage
    all_positional_data(:, :, :, mouse_idx) = current_positional_data;
end
%% Aligning positional data between mice with corner labels as reference data
data = readmatrix('D:\PhD 2nd Year\DeepLabCut Models\corner_labels.xlsx');
data = data(:, 4:end);

data_x = data(:, 1:2:end);
data_y = data(:, 2:2:end);

corner_matrix = zeros(4,2,num_mice);

for i = 1:num_mice
    corner_matrix(:, 1, i) = data_x(i, :)';
    corner_matrix(:, 2, i) = data_y(i, :)';
end

% Get dimensions
[num_dims, num_frames, num_trials, num_mice] = size(all_positional_data);

% Loop through each mouse to align their positional data
for mouse_idx = 1:num_mice
    % Extract this mouse's corner coordinates (4x2)
    mouse_corners = squeeze(corner_matrix(:, :, mouse_idx));

    % Compute Procrustes transformation using this mouse's own reference
    [~, ~, transform] = procrustes(mouse_corners, mouse_corners); % Identity transform

    % Loop through trials
    for trial_idx = 1:num_trials
        % Extract x and y positions for this mouse & trial (2x502)
        pos_data = squeeze(all_positional_data(:, :, trial_idx, mouse_idx));

        % Apply Procrustes transformation
        transformed_pos = transform.b * pos_data' * transform.T + transform.c(1, :);

        % Store back in all_positional_data
        all_positional_data(:, :, trial_idx, mouse_idx) = transformed_pos';
    end
end

all_positional_data = all_positional_data * 0.048;
% Converting positional data into centimeters

%% Splitting data

psilocybin_positional_data = all_positional_data(:, :, :, 2:2:20);
vehicle_positional_data = all_positional_data(:, :, :, 1:2:19);

%0 = flash, 1 = loom
stim_set_1 = [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1,...
    0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1];

stim_set_2 = [0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0,...
    1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0 ,0 ,1 ,0];

% defining pre and post stimulus frames
pre_stim_period = 2:152;
post_stim_period = 202:502;

% setting indexes for looms and flashes
stim_set_1_looms_idx = logical(stim_set_1);
stim_set_1_flashes_idx = ~stim_set_1_looms_idx;
stim_set_2_looms_idx = logical(stim_set_2);
stim_set_2_flashes_idx = ~stim_set_2_looms_idx;

psi_loom_ss1 = psilocybin_positional_data(:, :, stim_set_1_looms_idx, 1:2:9);
psi_flash_ss1 = psilocybin_positional_data(:, :, stim_set_1_flashes_idx, 1:2:9);
psi_loom_ss2 = psilocybin_positional_data(:, :, stim_set_2_looms_idx, 2:2:10);
psi_flash_ss2 = psilocybin_positional_data(:, :, stim_set_2_flashes_idx, 2:2:10);

veh_loom_ss1 = vehicle_positional_data(:, :, stim_set_1_looms_idx, 1:2:9);
veh_flash_ss1 = vehicle_positional_data(:, :, stim_set_1_flashes_idx, 1:2:9);
veh_loom_ss2 = vehicle_positional_data(:, :,  stim_set_2_looms_idx, 2:2:10);
veh_flash_ss2 = vehicle_positional_data(:, :, stim_set_2_flashes_idx, 2:2:10);

% concatenating all psi and veh freezing data
psi_loom_data = cat(4, psi_loom_ss1, psi_loom_ss2);
psi_flash_data = cat(4, psi_flash_ss1, psi_flash_ss2);
veh_loom_data = cat(4, veh_loom_ss1, veh_loom_ss2);
veh_flash_data = cat(4, veh_flash_ss1, veh_flash_ss2);

%extracting pre and post-stim freezing data
psi_post_stim_loom = psi_loom_data(:, post_stim_period, :, :);
psi_post_stim_flash = psi_flash_data(:, post_stim_period, :, :);
veh_post_stim_loom = veh_loom_data(:, post_stim_period, :, :);
veh_post_stim_flash = veh_flash_data(:, post_stim_period, :, :);

psi_pre_stim_loom = psi_loom_data(:, pre_stim_period, :, :);
psi_pre_stim_flash = psi_flash_data(:, pre_stim_period, :, :);
veh_pre_stim_loom = veh_loom_data(:, pre_stim_period, :, :);
veh_pre_stim_flash = veh_flash_data(:, pre_stim_period, :, :);

psi_pre_stim_data = cat(3, psi_pre_stim_loom, psi_pre_stim_flash);
veh_pre_stim_data = cat(3, veh_pre_stim_loom, veh_pre_stim_flash);

%% Plots traces of each mouse in arena (shows alignment has worked)
figure;
hold on;
axis equal;
title('Animal Positions in Arena (Aligned Data)');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');

% Loop through each mouse to plot their positional data
for mouse_idx = 1:2:num_mice
    % Loop through each trial for the current mouse
    for trial_idx = 1:num_trials
        % Extract x and y positions for the current mouse and trial
        pos_data = squeeze(all_positional_data(:, :, trial_idx, mouse_idx)); % 2x502

        % Plot the positions for this trial (x, y coordinates for each frame)
        plot(pos_data(1, :), pos_data(2, :), 'o-', 'MarkerSize', 3, 'LineWidth', 1);
    end
end

% Optionally, you can add specific colours for each stimulus (loom vs flash)
% Add more customisation options, like colours or labels, based on stimuli if desired

% Display the plot
hold off;

% Create a new figure for visualization
figure;
hold on;
axis equal;
title('Animal Positions in Arena (Aligned Data)');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');

% Loop through each mouse to plot their positional data
for mouse_idx = 2:2:num_mice
    % Loop through each trial for the current mouse
    for trial_idx = 1:num_trials
        % Extract x and y positions for the current mouse and trial
        pos_data = squeeze(all_positional_data(:, :, trial_idx, mouse_idx)); % 2x502

        % Plot the positions for this trial (x, y coordinates for each frame)
        plot(pos_data(1, :), pos_data(2, :), 'o-', 'MarkerSize', 3, 'LineWidth', 1);
    end
end

% Optionally, you can add specific colours for each stimulus (loom vs flash)
% Add more customisation options, like colours or labels, based on stimuli if desired

% Display the plot
hold off;
%% Plotting heatmaps of time spent in parts of the arena

% Create a figure for the heatmaps
% Create a figure for the heatmaps
figure;
subplot(1, 2, 1); % Create two subplots: one for vehicle and one for psilocybin

% Vehicle Post-Stimulus Heatmap
% Aggregating all positions of vehicle (pre-stim)
vehicle_x = reshape(veh_pre_stim_data(1, :, :, :), [], 1);  % Extract x positions
vehicle_y = reshape(veh_pre_stim_data(2, :, :, :), [], 1);  % Extract y positions

% Plot the heatmap
% Define the bin size for the heatmap (e.g., a larger number of bins for finer resolution)
nbins = [30, 30];  % Adjust this to increase/decrease bin size
[counts, xedges, yedges] = histcounts2(vehicle_x, vehicle_y, nbins);

% Plot using imagesc for a regular heatmap
imagesc(xedges, yedges, counts'); % Transpose counts to match x/y axes
set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
colorbar;
title('Vehicle Pre-Stim');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
clim([0, 1500]);
xlim([10, 60]);
ylim([5, 40]);

% Psilocybin Post-Stimulus Heatmap
% Aggregating all positions of psilocybin (pre-stim)
psilocybin_x = reshape(psi_pre_stim_data(1, :, :, :), [], 1);  % Extract x positions
psilocybin_y = reshape(psi_pre_stim_data(2, :, :, :), [], 1);  % Extract y positions

% Plot the heatmap
subplot(1, 2, 2); % Second subplot for psilocybin
[counts, xedges, yedges] = histcounts2(psilocybin_x, psilocybin_y, nbins);

% Plot using imagesc for a regular heatmap
imagesc(xedges, yedges, counts');
set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
colorbar;
title('Psilocybin Pre-Stim');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
clim([0, 1500]);
xlim([10, 60]);
ylim([5, 40]);

% Adjusting the layout
sgtitle('Heatmaps of Positions Pre-Stim');




%%
