%%
% load('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\positional_analysis\T_matrices.mat');
load('D:\PhD 2nd Year\T_matrices.mat')

% define constants
num_frames = 502;
num_trials = 40;
num_mice = 30;
num_dims = 2;

% Initialise storage for all positional data
all_positional_data = zeros(num_dims, num_frames, num_trials, num_mice);

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
        
        % Reshape and assign to the correct position in current_positional_data
        current_positional_data(1, :, trial_idx) = current_positional_data_x; % x positions
        current_positional_data(2, :, trial_idx) = current_positional_data_y; % y positions
    end
    
    % Assign data for the current mouse to the main storage
    all_positional_data(:, :, :, mouse_idx) = current_positional_data;
end
%% Aligning positional data between mice with corner labels as reference data

% data = readmatrix('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\positional_analysis\corner_labels.xlsx');
data = readmatrix('D:\PhD 2nd Year\DeepLabCut Models\corner_labels.xlsx');
data = data(:, 4:end);

data_x = data(:, 1:2:end);
data_y = data(:, 2:2:end);

corner_matrix = zeros(4,2,num_mice);

for i = 1:num_mice
    corner_matrix(:, 1, i) = data_x(i, :)';
    corner_matrix(:, 2, i) = data_y(i, :)';
end

ref_corner_matrix = corner_matrix(:,:,1);

% Loop through each mouse to align their positional data
for mouse_idx = 1:num_mice

    % Extract this mouse's arena corner coordinates (4x2)
    mouse_corners = squeeze(corner_matrix(:, :, mouse_idx));

    % Compute Procrustes transformation using this mouse's own reference
    [~, ~, transform] = procrustes(ref_corner_matrix, mouse_corners); % Identity transform

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


%% removing all positional data that is outside of arena

% Define arena boundaries (in cm) using reference mouse (mouse 1)
x_min = min(ref_corner_matrix(:,1)) * 0.048; 
x_max = max(ref_corner_matrix(:,1)) * 0.048; 
y_min = min(ref_corner_matrix(:,2)) * 0.048; 
y_max = max(ref_corner_matrix(:,2)) * 0.048;

% Loop through all mice and trials to remove outliers
for mouse_idx = 1:num_mice
    for trial_idx = 1:num_trials
        % Extract x and y positions for this mouse & trial
        x_data = squeeze(all_positional_data(1, :, trial_idx, mouse_idx));
        y_data = squeeze(all_positional_data(2, :, trial_idx, mouse_idx));

        % Find outliers
        outlier_idx = (x_data < x_min) | (x_data > x_max) | ...
                      (y_data < y_min) | (y_data > y_max);

        % Replace outliers with NaN
        all_positional_data(1, outlier_idx, trial_idx, mouse_idx) = NaN;
        all_positional_data(2, outlier_idx, trial_idx, mouse_idx) = NaN;
    end
end


%% Splitting data

psilocybin_positional_data = all_positional_data(:, :, :, 2:2:end);
vehicle_positional_data = all_positional_data(:, :, :, 1:2:end);

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

psi_loom_ss1 = psilocybin_positional_data(:, :, stim_set_1_looms_idx, 1:2:end);
psi_flash_ss1 = psilocybin_positional_data(:, :, stim_set_1_flashes_idx, 1:2:end);
psi_loom_ss2 = psilocybin_positional_data(:, :, stim_set_2_looms_idx, 2:2:end);
psi_flash_ss2 = psilocybin_positional_data(:, :, stim_set_2_flashes_idx, 2:2:end);

veh_loom_ss1 = vehicle_positional_data(:, :, stim_set_1_looms_idx, 1:2:end);
veh_flash_ss1 = vehicle_positional_data(:, :, stim_set_1_flashes_idx, 1:2:end);
veh_loom_ss2 = vehicle_positional_data(:, :,  stim_set_2_looms_idx, 2:2:end);
veh_flash_ss2 = vehicle_positional_data(:, :, stim_set_2_flashes_idx, 2:2:end);

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
nbins = [15, 15];  % Adjust this to increase/decrease bin size
[veh_counts, xedges, yedges] = histcounts2(vehicle_x, vehicle_y, nbins);

% Plot using imagesc for a regular heatmap
imagesc(xedges, yedges, veh_counts'); % Transpose counts to match x/y axes
set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
axis equal;
title('Vehicle Pre-Stim');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
clim([0,  7358]); % max(veh_counts(:))]???
xlim([x_min, x_max]);
ylim([y_min, y_max]);

% Get current tick positions
xt = xticks(); % Get x-axis tick positions
yt = yticks(); % Get y-axis tick positions

% Shift tick labels so they display from 0 onwards
xticklabels(xt - xt(1)); % Shift x labels
yticklabels(yt - yt(1)); % Shift y labels



% Psilocybin Post-Stimulus Heatmap
% Aggregating all positions of psilocybin (pre-stim)
psilocybin_x = reshape(psi_pre_stim_data(1, :, :, :), [], 1);  % Extract x positions
psilocybin_y = reshape(psi_pre_stim_data(2, :, :, :), [], 1);  % Extract y positions

% Plot the heatmap
subplot(1, 2, 2); % Second subplot for psilocybin
[psi_counts, xedges, yedges] = histcounts2(psilocybin_x, psilocybin_y, nbins);

% Plot using imagesc for a regular heatmap
imagesc(xedges, yedges, psi_counts');
set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
axis equal;
title('Psilocybin Pre-Stim');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
clim([0, 7358]); % max(psi_counts(:))???
xlim([x_min, x_max]);
ylim([y_min, y_max]);

% Get current tick positions
xt = xticks(); % Get x-axis tick positions
yt = yticks(); % Get y-axis tick positions

% Shift tick labels so they display from 0 onwards
xticklabels(xt - xt(1)); % Shift x labels
yticklabels(yt - yt(1)); % Shift y labels


% Create a single colorbar
c = colorbar('Position', [0.92, 0.29, 0.02, 0.44]); % Adjust [left, bottom, width, height]

% Convert colorbar tick labels from frames to seconds (assuming 15 fps)
colorbar_ticks = c.Ticks;  % Get the current colorbar ticks (which are frame counts)
colorbar_ticks_in_seconds = colorbar_ticks / 15;  % Convert frames to seconds

% Set the new tick labels to be in seconds
c.TickLabels = arrayfun(@(x) sprintf('%.0f s', x), colorbar_ticks_in_seconds, 'UniformOutput', false);

% Adjusting the layout
sgtitle('Heatmaps of Positions Pre-Stim');


%% create difference heatmap (Vehicle - Psilocybin)

difference_counts = veh_counts - psi_counts; % calculate difference in bin counts

% Plot using imagesc for a regular heatmap
figure
imagesc(xedges, yedges, difference_counts'); % Transpose counts to match x/y axes
set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
axis equal;
title('Difference Heatmap Veh - Psi (Pre-Stim)');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
clim([min(difference_counts(:)), max(difference_counts(:))]);
xlim([x_min, x_max]);
ylim([y_min, y_max]);

% Get current tick positions
xt = xticks(); % Get x-axis tick positions
yt = yticks(); % Get y-axis tick positions

% Shift tick labels so they display from 0 onwards
xticklabels(xt - xt(1)); % Shift x labels
yticklabels(yt - yt(1)); % Shift y labels

% Get current tick positions
xt = xticks(); % Get x-axis tick positions
yt = yticks(); % Get y-axis tick positions

% Shift tick labels so they display from 0 onwards
xticklabels(xt - xt(1)); % Shift x labels
yticklabels(yt - yt(1)); % Shift y labels

% Create a single colorbar
c = colorbar; %('Position', [0.92, 0.29, 0.02, 0.44]); % Adjust [left, bottom, width, height]

% Convert colorbar tick labels from frames to seconds (assuming 15 fps)
colorbar_ticks = c.Ticks;  % Get the current colorbar ticks (which are frame counts)
colorbar_ticks_in_seconds = colorbar_ticks / 15;  % Convert frames to seconds

% Set the new tick labels to be in seconds
c.TickLabels = arrayfun(@(x) sprintf('%.0f s', x), colorbar_ticks_in_seconds, 'UniformOutput', false);


