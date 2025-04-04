%%
 load('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\positional_analysis\T_matrices.mat');
% load('D:\PhD 2nd Year\T_matrices.mat')

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

 data = readmatrix('C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\positional_analysis\corner_labels.xlsx');
%data = readmatrix('D:\PhD 2nd Year\DeepLabCut Models\corner_labels.xlsx');
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

% Shift all positions so (x_min, y_min) becomes (0,0)
all_positional_data(1, :, :, :) = all_positional_data(1, :, :, :) - x_min;
all_positional_data(2, :, :, :) = all_positional_data(2, :, :, :) - y_min;

% defining new courner coodinates
  bot_left = [0, 0];
  top_left = [0,(y_max - y_min)];
  bot_right = [(x_max - x_min), 0];
  top_right = [(x_max - x_min), (y_max - y_min)];


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

%% extracting data (whole trial-period) of flash tials precursored by another flash trial (flash - FLASH)

% Identify indices of flash tials preceded by a flash (fF)
    ss1_fF_indx = [7, 8, 16, 17, 22, 23, 30, 34];
    ss2_fF_indx = [1, 12, 19, 20, 27, 28, 22, 28,]; 

% get freezing data from fF trials in veh and psi mice

psi_ss1_fF = psilocybin_positional_data(:, :, ss1_fF_indx, 1:2:end); % gets all the fF trials for psi mice that received stimset 1
psi_ss2_fF = psilocybin_positional_data(:, :, ss2_fF_indx, 2:2:end); % gets all the fF trials for psi mice that received stimset 2
veh_ss1_fF = vehicle_positional_data(:, :, ss1_fF_indx, 1:2:end);  % gets all the fF trials for psi mice that received stimset 1
veh_ss2_fF = vehicle_positional_data(:, :, ss2_fF_indx, 2:2:end);  % gets all the fF trials for psi mice that received stimset 2

% concatinate data frim ss1 and ss2 mice together
psi_fF = cat(4, psi_ss1_fF, psi_ss2_fF);
veh_fF = cat(4, veh_ss1_fF, veh_ss2_fF);


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
%% Plotting heatmaps of time spent in parts of the arena (Psi and Veh individually)

% choose which data you want analyse

  % vehicle_data = veh_fF;    % flash-FLASH trials veh
  % psilocybin_data = psi_fF; % flash-FLASH trials Psi

  vehicle_data = veh_post_stim_loom;    % loom trials (post-stim) veh
  psilocybin_data = psi_post_stim_loom; % loom trials (post-stim) Psi

% Create a figure for the heatmaps
% Create a figure for the heatmaps
figure;
subplot(1, 2, 1); % Create two subplots: one for vehicle and one for psilocybin

% Adjusting the layout
% sgtitle('Heatmaps of Positions in flash-FLASH trials ');
sgtitle('Positions in loom trials (post-stim)');

color_lim = 2400; %define bin count value that sets colourlimit

%%%% Vehicle Post-Stimulus Heatmap
% Aggregating all positions of vehicle (pre-stim)
vehicle_x = reshape(vehicle_data(1, :, :, :), [], 1);  % Extract x positions
vehicle_y = reshape(vehicle_data(2, :, :, :), [], 1);  % Extract y positions

% Define the bin size for the heatmap (e.g., a larger number of bins for finer resolution)
nbins = [15, 15];  % Adjust this to increase/decrease bin size
[veh_counts, xedges, yedges] = histcounts2(vehicle_x, vehicle_y, nbins, 'XBinLimits', [0, top_right(1)], 'YBinLimits', [0, top_right(2)]);
x_centers = (xedges(1:end-1) + xedges(2:end)) / 2;
y_centers = (yedges(1:end-1) + yedges(2:end)) / 2;
imagesc(x_centers, y_centers, veh_counts');

set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
axis equal;
title('Vehicle Pre-Stim');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
clim([0,  color_lim]); % max(veh_counts(:))]???
xlim([0, top_right(1)]);
ylim([0, top_right(2)]);

% Psilocybin Post-Stimulus Heatmap
% Aggregating all positions of psilocybin (pre-stim)
psilocybin_x = reshape(psilocybin_data(1, :, :, :), [], 1);  % Extract x positions
psilocybin_y = reshape(psilocybin_data(2, :, :, :), [], 1);  % Extract y positions

% Plot the heatmap
subplot(1, 2, 2); % Second subplot for psilocybin
[psi_counts, xedges, yedges] = histcounts2(psilocybin_x, psilocybin_y, nbins, 'XBinLimits', [0, top_right(1)], 'YBinLimits', [0, top_right(2)]);
x_centers = (xedges(1:end-1) + xedges(2:end)) / 2;
y_centers = (yedges(1:end-1) + yedges(2:end)) / 2;
imagesc(x_centers, y_centers, psi_counts');

set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
axis equal;
title('Psilocybin Pre-Stim');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
clim([0, color_lim]); % max(psi_counts(:))???
xlim([0, top_right(1)]);
ylim([0, top_right(2)]);


% Create a single colorbar
c = colorbar('Position', [0.92, 0.29, 0.02, 0.44]); % Adjust [left, bottom, width, height]

% colormap(parula); % Default MATLAB colormap    
%colormap(jet); % Rainbow-like colors
%colormap(hot); % Red-yellow gradient
%colormap(gray); % Grayscale
%colormap(cool); % Cyan-magenta gradient
% colormap(hsv); % Hue-saturation-value color wheel
colormap(turbo); % High-contrast, perceptually uniform
% make on green-blue-yellow

% Convert colorbar tick labels from frames to seconds (assuming 15 fps)
colorbar_ticks = c.Ticks;  % Get the current colorbar ticks (which are frame counts)
colorbar_ticks_in_seconds = colorbar_ticks / 15;  % Convert frames to seconds

% Set the new tick labels to be in seconds
c.TickLabels = arrayfun(@(x) sprintf('%.0f s', x), colorbar_ticks_in_seconds, 'UniformOutput', false);


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
clim([-500, 500]); % ([min(difference_counts(:)), max(difference_counts(:))]);
xlim([0, top_right(1)]);
ylim([0, top_right(2)]);


% Create a single colorbar
c = colorbar; %('Position', [0.92, 0.29, 0.02, 0.44]); % Adjust [left, bottom, width, height]

% Convert colorbar tick labels from frames to seconds (assuming 15 fps)
colorbar_ticks = c.Ticks;  % Get the current colorbar ticks (which are frame counts)
colorbar_ticks_in_seconds = colorbar_ticks / 15;  % Convert frames to seconds

% colormap(parula); % Default MATLAB colormap
%colormap(jet); % Rainbow-like colors
%colormap(hot); % Red-yellow gradient
%colormap(gray); % Grayscale
%colormap(cool); % Cyan-magenta gradient
% colormap(hsv); % Hue-saturation-value color wheel
%colormap(turbo); % High-contrast, perceptually uniform
% 

% Set the new tick labels to be in seconds
c.TickLabels = arrayfun(@(x) sprintf('%.0f s', x), colorbar_ticks_in_seconds, 'UniformOutput', false);


% %% Calculate binsize in cm to define courner bins
% 
% % Compute bin width and height in cm
% bin_width = (x_max - x_min) / nbins(1);  % Width of each bin in cm
% bin_height = (y_max - y_min) / nbins(2); % Height of each bin in cm
% 
% % Display results
% fprintf('Bin Width: %.2f cm\n', bin_width);
% fprintf('Bin Height: %.2f cm\n', bin_height);
% 
% % calculate how many columns and rows of bin-matrix fall into courners
% mouse_cm = 8; % average mouse size = 6-10 cm
% row_nr = round(mouse_cm / bin_height); 
% column_nr = round(mouse_cm / bin_width); 


%% determining time spent along wall and in center in flash-FLASH trials (psi vs vehicle)

% define boundaries of center of arena (in cm)
  edge_widht = 4; % width of the area around all the walls that shall be counted as thigmotactic area (mouse is touching walls)
  x_cen_min = edge_widht;
  y_cen_min = edge_widht;
  x_cen_max = top_right(1) - edge_widht;
  y_cen_max = top_right(2) - edge_widht;

% define which frames mouse is spending time in periphery of center of arena 

  thigmotaxis_matrix_psi = zeros(502,8,15); % later all frames where mouse is in center filled with 1s
  thigmotaxis_matrix_veh = zeros(502,8,15); 

% loop trough all frames in all trials and all mice 
  for m = 1: size(thigmotaxis_matrix_psi,3)

    for t = 1: size(thigmotaxis_matrix_psi,2)

      for f = 1: size(thigmotaxis_matrix_psi,1)
          
          % for psilocybin data
          if x_cen_min < psi_fF(1,f,t,m) && psi_fF(1,f,t,m) < x_cen_max && ...  % if mouse coordinates fall within center boundaries...
             y_cen_min < psi_fF(2,f,t,m) && psi_fF(2,f,t,m) < y_cen_max     
        
             thigmotaxis_matrix_psi(f,t,m) = 1;                                 % ... set to 1
          end
          
          % for vehicle data
          if x_cen_min < veh_fF(1,f,t,m) && veh_fF(1,f,t,m) < x_cen_max && ...  % if mouse coordinates fall within center boundaries...
             y_cen_min < veh_fF(2,f,t,m) && veh_fF(2,f,t,m) < y_cen_max     
        
             thigmotaxis_matrix_veh(f,t,m) = 1;                                  % ... set to 1
          end
      end
    end
  end

% calculate total time spent in periphery and center per mouse
  
 % Initialize counters for total frames spent in center and periphery
    psi_center = zeros(1, size(thigmotaxis_matrix_psi, 3));
    psi_periphery = zeros(1, size(thigmotaxis_matrix_psi, 3));
    
    veh_center = zeros(1, size(thigmotaxis_matrix_veh, 3));
    veh_periphery = zeros(1, size(thigmotaxis_matrix_veh, 3));

% Loop through each mouse
for m = 1:size(thigmotaxis_matrix_psi, 3)
    % Count frames in center and periphery for psilocybin
    psi_center(m) = sum(thigmotaxis_matrix_psi(:,:,m) == 1, 'all'); % Center (1 indicates center)
    psi_periphery(m) = sum(thigmotaxis_matrix_psi(:,:,m) == 0, 'all'); % Periphery (0 indicates periphery)
    
    % Count frames in center and periphery for vehicle
    veh_center(m) = sum(thigmotaxis_matrix_veh(:,:,m) == 1, 'all'); % Center
    veh_periphery(m) = sum(thigmotaxis_matrix_veh(:,:,m) == 0, 'all'); % Periphery
end

% convert to seconds (15 fps)
  psi_center = psi_center / 15;
  psi_periphery = psi_periphery / 15;
  veh_center    = veh_center / 15;
  veh_periphery =  veh_periphery / 15;

%% determining time spent along wall and in center in loom trials (psi vs vehicle)

% define boundaries of center of arena (in cm)
  edge_widht = 6; % width of the area around all the walls that shall be counted as thigmotactic area (mouse is touching walls)
  x_cen_min = edge_widht;
  y_cen_min = edge_widht;
  x_cen_max = top_right(1) - edge_widht;
  y_cen_max = top_right(2) - edge_widht;

% define which frames mouse is spending time in periphery of center of arena 

  thigmotaxis_matrix_psi = zeros(502,8,15); % later all frames where mouse is in center filled with 1s
  thigmotaxis_matrix_veh = zeros(502,8,15); 

% loop trough all frames in all trials and all mice 
  for m = 1: size(thigmotaxis_matrix_psi,3)

    for t = 1: size(thigmotaxis_matrix_psi,2)

      for f = 1: size(thigmotaxis_matrix_psi,1)
          
          % for psilocybin data
          if x_cen_min < psi_fF(1,f,t,m) && psi_fF(1,f,t,m) < x_cen_max && ...  % if mouse coordinates fall within center boundaries...
             y_cen_min < psi_fF(2,f,t,m) && psi_fF(2,f,t,m) < y_cen_max     
        
             thigmotaxis_matrix_psi(f,t,m) = 1;                                 % ... set to 1
          end
          
          % for vehicle data
          if x_cen_min < veh_fF(1,f,t,m) && veh_fF(1,f,t,m) < x_cen_max && ...  % if mouse coordinates fall within center boundaries...
             y_cen_min < veh_fF(2,f,t,m) && veh_fF(2,f,t,m) < y_cen_max     
        
             thigmotaxis_matrix_veh(f,t,m) = 1;                                  % ... set to 1
          end
      end
    end
  end

% calculate total time spent in periphery and center per mouse
  
 % Initialize counters for total frames spent in center and periphery
    psi_center = zeros(1, size(thigmotaxis_matrix_psi, 3));
    psi_periphery = zeros(1, size(thigmotaxis_matrix_psi, 3));
    
    veh_center = zeros(1, size(thigmotaxis_matrix_veh, 3));
    veh_periphery = zeros(1, size(thigmotaxis_matrix_veh, 3));

% Loop through each mouse
for m = 1:size(thigmotaxis_matrix_psi, 3)
    % Count frames in center and periphery for psilocybin
    psi_center(m) = sum(thigmotaxis_matrix_psi(:,:,m) == 1, 'all'); % Center (1 indicates center)
    psi_periphery(m) = sum(thigmotaxis_matrix_psi(:,:,m) == 0, 'all'); % Periphery (0 indicates periphery)
    
    % Count frames in center and periphery for vehicle
    veh_center(m) = sum(thigmotaxis_matrix_veh(:,:,m) == 1, 'all'); % Center
    veh_periphery(m) = sum(thigmotaxis_matrix_veh(:,:,m) == 0, 'all'); % Periphery
end

% convert to seconds (15 fps)
  psi_center = psi_center / 15;
  psi_periphery = psi_periphery / 15;
  veh_center    = veh_center / 15;
  veh_periphery =  veh_periphery / 15;

%% Statitistical analysis of thigmotaxis behaviour

% test if data normally distributed (h = 0 yes, h = 1 no)
    % For time spent in center
    [h_psi_c, p_psi_c] = lillietest(psi_center); % 
    [h_veh_c, p_veh_c] = lillietest(veh_center); % 
    % For time spent in periphery
    [h_psi_p, p_psi_p] = lillietest(psi_periphery); % 
    [h_veh_p, p_veh_p] = lillietest(veh_periphery); % 

% % IF NOT NROMALLY DISTRIBUTED: Mann–Whitney U test 
%     % Compare time spent in center
%     [p_center, h_center, stats_center] = ranksum(psi_center, veh_center);
%     % Compare time spent in periphery
%     [p_periphery, h_periphery, stats_periphery] = ranksum(psi_periphery, veh_periphery);
%     % Display results
%         fprintf('Mann–Whitney U test for time spent in center: p-value = %.4f\n', p_center);
%         fprintf('Mann–Whitney U test for time spent in periphery: p-value = %.4f\n', p_periphery);


  
 % IF NOT NROMALLY DISTRIBUTED: independent t-tests
    % Independent t-test for time spent in center
    [center_h, center_p] = ttest2(psi_center, veh_center); % Independent t-test for center timee   
    % Independent t-test for time spent in periphery
    [periphery_h, periphery_p] = ttest2(psi_periphery, veh_periphery); % Independent t-test for periphery time   
    % Display results
    fprintf('Independent t-test for time spent in center: p-value = %.4f\n', center_p);
    fprintf('Independent t-test for time spent in periphery: p-value = %.4f\n', periphery_p);


% plot resuts     
% Create a figure for the plots
figure;

% Plot time spent in center
subplot(1, 2, 1); % Create subplot (1 row, 2 columns, first plot)
hold on;

% Bar plot with error bars
bar([mean(psi_center), mean(veh_center)], 'FaceColor', [0.2, 0.6, 1]); % Psilocybin and Vehicle
errorbar([1, 2], [mean(psi_center), mean(veh_center)], ...
         [std(psi_center)/sqrt(length(psi_center)), std(veh_center)/sqrt(length(veh_center))], ...
         'k', 'LineStyle', 'none'); % Standard error bars

% Set plot properties
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Psilocybin', 'Vehicle'});
ylabel('Time Spent in Center (s)');
title('Time Spent in Center');

hold off;

% Plot time spent in periphery
subplot(1, 2, 2); % Second plot
hold on;

% Bar plot with error bars
bar([mean(psi_periphery), mean(veh_periphery)], 'FaceColor', [0.9, 0.4, 0.4]); % Psilocybin and Vehicle
errorbar([1, 2], [mean(psi_periphery), mean(veh_periphery)], ...
         [std(psi_periphery)/sqrt(length(psi_periphery)), std(veh_periphery)/sqrt(length(veh_periphery))], ...
         'k', 'LineStyle', 'none'); % Standard error bars

% Set plot properties
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Psilocybin', 'Vehicle'});
ylabel('Time Spent in Periphery (s)');
title('Time Spent in Periphery');

hold off;

