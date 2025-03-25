% Specify the folder containing your .mat files
folderPath = 'C:\PhD 2nd Year\DLC Renewal Data';

% Get a list of all .csv files in the specified folder
fileList = dir(fullfile(folderPath, '*.csv'));

% Initialise a cell array to store each T matrix
T_matrices = cell(1, length(fileList));
file_names = (cell(length(fileList), 1));

% iterate over each file, load it, and extract the T variable
for i = 1:length(fileList)
    filePath = fullfile(folderPath, fileList(i).name);
     % load the .mat file
    data = readmatrix(filePath); 
    T_matrices{i} = data(:, [14, 15])'; 
    % stores body_anterior x and y positional data in T_matrices array
    file_names{i} = fileList(i).name;
    % stores file names for sorting later
end


%% Sorting velocity matrix into mouse/part order

mouse_numbers = cellfun(@(x) sscanf(x, 'mouse%d'), file_names);
% extracting the mouse number from the file names

[~, sort_idx] = sort(mouse_numbers);
% sort by mouse number amd get index for sorted order

T_matrices = T_matrices(sort_idx);
% sort matrix
disp('Sorted File Names:');
disp(file_names(sort_idx));
%% Putting parts of videos together for each mouse
T_matrices_paired = cell(1, length(T_matrices)/2);
%initialising new cell array half the size of original

for video_ix = 1:2:length(T_matrices)
    %iterating over pairs of cells
    T_matrices_paired{(video_ix+1)/2} = cat(2, T_matrices{video_ix}, T_matrices{video_ix+1});
    %stores concatenated matrices in correct indices of new matrix
end

T_matrices = T_matrices_paired;
%% Getting xy positional data and calculating Euclidian distance between frames

num_frames = 502; % 502 frames per trial
num_trials = 40; % 40 trials
num_mice = 29; % change as needed
all_velocity_data = zeros(num_trials, num_frames, num_mice);

for mouse_idx = 1:length(T_matrices)
    x_pos = T_matrices{mouse_idx}(1,:);
    y_pos = T_matrices{mouse_idx}(2,:);

    dist = zeros(40, 502);
    
    for trial_ix = 1:size(dist,1)
        dx = [0, diff(x_pos(((trial_ix-1)*num_frames + 1):(trial_ix*num_frames)))];
        dy = [0, diff(y_pos(((trial_ix-1)*num_frames + 1):(trial_ix*num_frames)))];
        % calculating dx and dy for all frames for each trial (e.g., 1:502,
        % 503:1004 etc. and adds a 0 at start to maintain correct trial length)
        dist(trial_ix, :) = sqrt(dx.^2 + dy.^2);
        % calcuting Euclidian distance and adding values to columns of each
        % trial row
    end

all_velocity_data(:, :, mouse_idx) = dist;

end

%% Defining freezing behaviour

freezing_velocity_threshold = 1.1; % pixels per frame, change as needed (*0.048*15 to get cm/s value)
freezing_duration_threshold = 1*15; % convert seconds to frames (1 second × 15 FPS)

freeze_counter_matrix = all_velocity_data < freezing_velocity_threshold;

validated_freeze_matrix = zeros(size(freeze_counter_matrix));

for mouse_ix = 1:size(validated_freeze_matrix, 3)
    %iterating over each video
    for trials = 1:size(freeze_counter_matrix,1)
        % iterating over each trial
        [labeled_freezing, num_segments] = bwlabel(freeze_counter_matrix(trials, :, mouse_ix));
        % labeled_freezing is a vector same size as freeze_counter_matrix, but
        % has a unique value for each 'blob' of consecutive trials
        %num_segments is as scalar representing how many 'blobs' there are
        segment_props = regionprops(labeled_freezing, 'Area');
        % extracts length of connected regions in labeled_freezing and stores
        % it in segment_props.Area
        for segment = 1:num_segments
            %iterates over number of consecutive segments
            if segment_props(segment).Area >=freezing_duration_threshold
                validated_freeze_matrix(trials, labeled_freezing == segment, mouse_ix) = 1; % checks if current segment of freezing meets duration threshold, if it does, sets current segment's values to 1 
            end
        end
    end
end
%% Splitting into treatment groups, stim types and trial periods

psilocybin_freezing_data = validated_freeze_matrix(:, :, 1:2:end);
vehicle_freezing_data = validated_freeze_matrix(:, :, 2:2:end);

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

psi_loom_ss1 = psilocybin_freezing_data(stim_set_1_looms_idx, :, 1:2:end);
psi_flash_ss1 = psilocybin_freezing_data(stim_set_1_flashes_idx, :, 1:2:end);
psi_loom_ss2 = psilocybin_freezing_data(stim_set_2_looms_idx, :, 2:2:end);
psi_flash_ss2 = psilocybin_freezing_data(stim_set_2_flashes_idx, :, 2:2:end);

veh_loom_ss1 = vehicle_freezing_data(stim_set_1_looms_idx, :, 2:2:end);
veh_flash_ss1 = vehicle_freezing_data(stim_set_1_flashes_idx, :, 2:2:end);
veh_loom_ss2 = vehicle_freezing_data(stim_set_2_looms_idx, :, 1:2:end);
veh_flash_ss2 = vehicle_freezing_data(stim_set_2_flashes_idx, :, 1:2:end);

% concatenating all psi and veh freezing data
psi_loom_data = interleave(psi_loom_ss1, psi_loom_ss2);
psi_flash_data = interleave(psi_flash_ss1, psi_flash_ss2);
veh_loom_data = interleave_renewal(veh_loom_ss1, veh_loom_ss2);
veh_flash_data = interleave_renewal(veh_flash_ss1, veh_flash_ss2);

%extracting pre and post-stim freezing data
psi_post_stim_loom = psi_loom_data(:, post_stim_period, :);
psi_post_stim_flash = psi_flash_data(:, post_stim_period, :);
veh_post_stim_loom = veh_loom_data(:, post_stim_period, :);
veh_post_stim_flash = veh_flash_data(:, post_stim_period, :);

psi_pre_stim_loom = psi_loom_data(:, pre_stim_period, :);
psi_pre_stim_flash = psi_flash_data(:, pre_stim_period, :);
veh_pre_stim_loom = veh_loom_data(:, pre_stim_period, :);
veh_pre_stim_flash = veh_flash_data(:, pre_stim_period, :);
%% plotting mean freezing after stimulus per trial

fps = 15;
% defining frame rate

% Calculating mean freezing time (in seconds) and SEM
psi_loom_trial_means = mean(sum(psi_post_stim_loom, 2) / fps, 3);
psi_flash_trial_means = mean(sum(psi_post_stim_flash, 2) / fps, 3);
veh_loom_trial_means = mean(sum(veh_post_stim_loom, 2) / fps, 3);
veh_flash_trial_means = mean(sum(veh_post_stim_flash, 2) / fps, 3);

psi_loom_trial_sems = std(sum(psi_post_stim_loom, 2) / fps, [], 3) / sqrt(size(psi_post_stim_loom, 3));
psi_flash_trial_sems = std(sum(psi_post_stim_flash, 2) / fps, [], 3) / sqrt(size(psi_post_stim_flash, 3));
veh_loom_trial_sems = std(sum(veh_post_stim_loom, 2) / fps, [], 3) / sqrt(size(veh_post_stim_loom, 3));
veh_flash_trial_sems = std(sum(veh_post_stim_flash, 2) / fps, [], 3) / sqrt(size(veh_post_stim_flash, 3));

trials = 1:20;
% trial numbers

% Plot Psilocybin Group
figure;
hold on;
errorbar(trials, psi_loom_trial_means, psi_loom_trial_sems, '-o', 'DisplayName', 'Loom');
errorbar(trials, psi_flash_trial_means, psi_flash_trial_sems, '-x', 'DisplayName', 'Flash');
hold off;
title('Freezing Behaviour Across Repeated Trials: Psilocybin');
legend('Location', 'Best');
grid on;

% Plot Vehicle Group
figure;
hold on;
errorbar(trials, veh_loom_trial_means, veh_loom_trial_sems, '-o', 'DisplayName', 'Loom');
errorbar(trials, veh_flash_trial_means, veh_flash_trial_sems, '-x', 'DisplayName', 'Flash');
hold off;
title('Freezing Behaviour Across Repeated Trials: Vehicle');
legend('Location', 'Best');
grid on;

% Set common axis labels
for i = 1:2
    figure(i);
    xlabel('Trial Number');
    ylabel('Mean Time Spent Freezing Post-Stimulus (s)');
    % Set consistent axis limits (modify based on data range)
    xlim([1, 20]);
    ylim([0, 14])
end


%% plotting group variation in freezing

%Plotting median of mean freezing probability between all mice over all
%loom trials (post-stim period only)
mean_freezing_psi_loom_post = squeeze(mean(psi_post_stim_loom, 1));
mean_freezing_veh_loom_post = squeeze(mean(veh_post_stim_loom, 1));

figure;
plot(median(mean_freezing_psi_loom_post, 2))
hold on 
plot(median(mean_freezing_veh_loom_post, 2))
xlabel('Frame');
ylabel('Median Freezing Probability')
legend('Psilocybin', 'Vehicle');
title('Freezing Per Frame, All Mice (Post-Loom Period Only)')
xlim([1,301]);
ylim([0, 1]);
hold off


%Plotting median of mean freezing probability between all mice over all
%loom trials (post-stim period only)
mean_freezing_psi_flash_post = squeeze(mean(psi_post_stim_flash, 1));
mean_freezing_veh_flash_post = squeeze(mean(veh_post_stim_flash, 1));

figure;
plot(median(mean_freezing_psi_flash_post, 2))
hold on 
plot(median(mean_freezing_veh_flash_post, 2))
xlabel('Frame');
ylabel('Median Freezing Probability')
legend('Psilocybin', 'Vehicle');
title('Freezing Per Frame, All Mice (Post-Flash Period Only)')
xlim([1,301]);
ylim([0, 1]);
hold off


%% Plotting individual mean freezing probability between all mice over all loom trials (post-stim period only)
figure;
% Using imagesc to plot the data
imagesc([mean_freezing_veh_loom_post'; mean_freezing_psi_loom_post']); 
colorbar; % Add colour bar

% Add titles and labels
title('Freezing Per Frame Per Mouse (Post-Loom Period Only)'); % Graph title
xlabel('Frame'); % x-axis label
ylabel('Mice'); % y-axis label

% Label for the colour bar
c = colorbar;
c.Label.String = 'Mean Freezing Probability';

% Adjust the axis for a better fit
axis tight; % Fit the axis tightly around the data

% Modify y-axis ticks and labels to show group information
yticks(1:30); % Set y-axis ticks for each mouse
yticklabels({'M1 (V)', 'M3 (V)', 'M5 (V)', 'M7 (V)', 'M9 (V)', 'M11 (V)', 'M13 (V)', 'M15(V)', 'M17 (V)', 'M19 (V)', 'M21 (V)', 'M23 (V)', 'M25 (V)', 'M27 (V)', 'M29 (V)', ...
             'M2 (P)', 'M4 (P)', 'M6 (P)', 'M8 (P)', 'M10 (P)', 'M12 (P)', 'M14 (P)', 'M16 (P)', 'M18 (P)', 'M20 (P)', 'M22 (P)', 'M24 (P)', 'M26 (P)', 'M28 (P)', 'M30 (P)'});



%% Plotting mean freezing probability between all mice over all loom trials
% (all frames)
mean_freezing_psi = squeeze(mean(psi_loom_data, 1));
mean_freezing_veh = squeeze(mean(veh_loom_data, 1));

figure;
% Using imagesc to plot the data
imagesc([mean_freezing_veh'; mean_freezing_psi']); 
colorbar; % Add colour bar

% Adding stimulus markers
hold on;
xline(152, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus Start'); % Stimulus start
xline(202, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus End');   % Stimulus end
hold off;

% Add titles and labels
title('Freezing Per Frame Per Mouse (Whole Loom Trial)'); % Graph title
xlabel('Frame'); % x-axis label
ylabel('Mice (Psilocybin Group and Vehicle Group)'); % y-axis label

% Colour bar labelling
c = colorbar;
c.Label.String = 'Mean Freezing Probability';
clim([0,1])

% Adjust the axis for better appearance
axis tight; % Fit the axis tightly around the data

yticks(1:30); % Set y-axis ticks for each mouse
yticklabels({'M3 (V)', 'M5 (V)', 'M7 (V)', 'M9 (V)', '11 (V)', 'M13 (V)', 'M15 (V)', 'M17(V)', 'M19 (V)', 'M21 (V)', 'M23 (V)', 'M25 (V)', 'M27 (V)', 'M29 (V)', ...
             'M2 (P)', 'M4 (P)', 'M6 (P)', 'M8 (P)', 'M10 (P)', 'M12 (P)', 'M14 (P)', 'M16 (P)', 'M18 (P)', 'M20 (P)', 'M22 (P)', 'M24 (P)', 'M26 (P)', 'M28 (P)', 'M30 (P)'});

%% Plotting mean freezing probability between all mice over all loom trials (split into two graphs)
% (all frames)
mean_freezing_psi = squeeze(mean(psi_loom_data, 1));
mean_freezing_veh = squeeze(mean(veh_loom_data, 1));

figure;
% Using imagesc to plot the data
imagesc(mean_freezing_veh'); 
colorbar; % Add colour bar

% Adding stimulus markers
hold on;
xline(152, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus Start'); % Stimulus start
xline(202, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus End');   % Stimulus end
hold off;

% Add titles and labels
title('Freezing Per Frame Per Mouse (Whole Loom Trial: Vehicle)'); % Graph title
xlabel('Frame'); % x-axis label
ylabel('Mice'); % y-axis label

% Colour bar labelling
c = colorbar;
c.Label.String = 'Mean Freezing Probability';

% Adjust the axis for better appearance
axis tight; % Fit the axis tightly around the data

figure;
% Using imagesc to plot the data
imagesc(mean_freezing_psi'); 
colorbar; % Add colour bar

% Adding stimulus markers
hold on;
xline(152, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus Start'); % Stimulus start
xline(202, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus End');   % Stimulus end
hold off;

% Add titles and labels
title('Freezing Per Frame Per Mouse (Whole Loom Trial: Psilocybin)'); % Graph title
xlabel('Frame'); % x-axis label
ylabel('Mice'); % y-axis label

% Colour bar labelling
c = colorbar;
c.Label.String = 'Mean Freezing Probability';

% Adjust the axis for better appearance
axis tight; % Fit the axis tightly around the data

%% Comparing early and late trials and pre and post stimulus periods (% of time spent freezing)


%
psi_loom_avg = reshape(mean(psi_loom_data,2), 20, 1, 10);
psi_flash_avg = reshape(mean(psi_flash_data,2), 20, 1, 10);
veh_loom_avg = reshape(mean(veh_loom_data,2), 20, 1 ,10);
veh_flash_avg = reshape(mean(veh_flash_data,2), 20, 1, 10);

response_matrix = cat(2, psi_loom_avg, psi_flash_avg,...
    veh_loom_avg, veh_flash_avg);

early_responses = mean(response_matrix(1:10, :, :), 1); % Trials 1–10
late_responses = mean(response_matrix(11:20, :, :), 1); % Trials 11–20

% Aggregate across mice
early_mean = squeeze(mean(early_responses, 3)); % [4 conditions]
late_mean = squeeze(mean(late_responses, 3));

% Create a grouped bar plot
figure;
hBar = bar([early_mean; late_mean]', 'grouped');  % Transpose to match the structure

% Set colors for the bars
hBar(1).FaceColor = 'b';  % Blue for early trials
hBar(2).FaceColor = 'g';  % Green for late trials

% Add axis labels
xlabel('Condition (Psilocybin Loom, Psilocybin Flash, Vehicle Loom, Vehicle Flash)');
ylabel('Mean Response (Freezing Probability)');

% Add a legend
legend({'Early Trials', 'Late Trials'}, 'Location', 'northeast');

% Set the title
title('Early vs Late Responses Across Conditions');

%% Plotting proportion of time spent freezing per trial, with habituation freezing also shown
fps = 15; % defining frame rate
total_duration = 20; % total duration of the observation period in seconds

% Calculating mean proportion of freezing time and SEM
psi_loom_trial_means = mean(sum(psi_post_stim_loom, 2) / fps / total_duration, 3);
psi_flash_trial_means = mean(sum(psi_post_stim_flash, 2) / fps / total_duration, 3);
veh_loom_trial_means = mean(sum(veh_post_stim_loom, 2) / fps / total_duration, 3);
veh_flash_trial_means = mean(sum(veh_post_stim_flash, 2) / fps / total_duration, 3);

psi_loom_trial_sems = std(sum(psi_post_stim_loom, 2) / fps / total_duration, [], 3) / sqrt(size(psi_post_stim_loom, 3));
psi_flash_trial_sems = std(sum(psi_post_stim_flash, 2) / fps / total_duration, [], 3) / sqrt(size(psi_post_stim_flash, 3));
veh_loom_trial_sems = std(sum(veh_post_stim_loom, 2) / fps / total_duration, [], 3) / sqrt(size(veh_post_stim_loom, 3));
veh_flash_trial_sems = std(sum(veh_post_stim_flash, 2) / fps / total_duration, [], 3) / sqrt(size(veh_post_stim_flash, 3));

trials = 1:20; % trial numbers

% Plot Psilocybin Group
figure;
hold on;
errorbar(trials, psi_loom_trial_means, psi_loom_trial_sems, '-o', 'DisplayName', 'Loom');
errorbar(trials, psi_flash_trial_means, psi_flash_trial_sems, '-x', 'DisplayName', 'Flash');
hab_line_psi = yline(mean_psilocybin_habituation, '--', 'DisplayName', 'Habituation', ...
    'LabelHorizontalAlignment', 'left', 'Color', 'k');
hold off;
title('Proportion of Freezing Behaviour Across Trials: Psilocybin');
legend('Location', 'Best');
grid on;

% Plot Vehicle Group
figure;
hold on;
errorbar(trials, veh_loom_trial_means, veh_loom_trial_sems, '-o', 'DisplayName', 'Loom');
errorbar(trials, veh_flash_trial_means, veh_flash_trial_sems, '-x', 'DisplayName', 'Flash');
hab_line_veh = yline(mean_vehicle_habituation, '--', 'DisplayName', 'Habituation', ...
    'LabelHorizontalAlignment', 'left', 'Color', 'k');
hold off;
title('Proportion of Freezing Behaviour Across Trials: Vehicle');
legend('Location', 'Best');
grid on;

% Set common axis labels
for i = 1:2
    figure(i);
    xlabel('Trial Number');
    ylabel('Proportion of Time Spent Freezing Post-Stimulus');
    % Set consistent axis limits (modify based on data range)
    xlim([1, 20]);
    ylim([0, 1]);
end
%% Comparing stim sets of same treatment group to investigate order effect
fps = 15; % defining frame rate

% Calculating mean freezing time (in seconds) for each trial
% You might want to average over the trials and stimulus sets
psi_loom_ss1_mean = mean(sum(psi_loom_ss1, 2) / fps, 3); % Psi Loom Set 1
psi_loom_ss2_mean = mean(sum(psi_loom_ss2, 2) / fps, 3); % Psi Loom Set 2
veh_loom_ss1_mean = mean(sum(veh_loom_ss1, 2) / fps, 3); % Vehicle Loom Set 1
veh_loom_ss2_mean = mean(sum(veh_loom_ss2, 2) / fps, 3); % Vehicle Loom Set 2

% Check if you are summing over the correct dimension (frames)
% Typically, frame-based calculations would sum over the 2nd dimension (time axis), not trials.

% Calculate SEM for each set
psi_loom_ss1_sem = std(sum(psi_loom_ss1, 2) / fps, [], 3) / sqrt(size(psi_loom_ss1, 3));
psi_loom_ss2_sem = std(sum(psi_loom_ss2, 2) / fps, [], 3) / sqrt(size(psi_loom_ss2, 3));
veh_loom_ss1_sem = std(sum(veh_loom_ss1, 2) / fps, [], 3) / sqrt(size(veh_loom_ss1, 3));
veh_loom_ss2_sem = std(sum(veh_loom_ss2, 2) / fps, [], 3) / sqrt(size(veh_loom_ss2, 3));

% Plot Psilocybin Group: Loom Stimulus Set 1 vs Stimulus Set 2
figure;
hold on;
errorbar(trials, psi_loom_ss1_mean, psi_loom_ss1_sem, 'o-', 'DisplayName', 'Stimulus Set 1', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
errorbar(trials, psi_loom_ss2_mean, psi_loom_ss2_sem, 'x-', 'DisplayName', 'Stimulus Set 2', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
hold off;
title('Psilocybin Group: Loom');
xlabel('Trial Number');
ylabel('Mean Time Spent Freezing Post-Stimulus (s)');
legend('Location', 'Best');
grid on;

% Plot Vehicle Group: Loom Stimulus Set 1 vs Stimulus Set 2
figure;
hold on;
errorbar(trials, veh_loom_ss1_mean, veh_loom_ss1_sem, 'o-', 'DisplayName', 'Stimulus Set 1', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
errorbar(trials, veh_loom_ss2_mean, veh_loom_ss2_sem, 'x-', 'DisplayName', 'Stimulus Set 2', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
hold off;
title('Vehicle Group: Loom');
xlabel('Trial Number');
ylabel('Mean Time Spent Freezing Post-Stimulus (s)');
legend('Location', 'Best');
grid on;

% Set consistent axis limits
for i = 1:2
    figure(i);
    % Set consistent axis limits (modify based on data range)
    xlim([1, 20]);
    ylim([0, max([psi_loom_ss1_mean + psi_loom_ss1_sem, ...
                  psi_loom_ss2_mean + psi_loom_ss2_sem, ...
                  veh_loom_ss1_mean + veh_loom_ss1_sem, ...
                  veh_loom_ss2_mean + veh_loom_ss2_sem], [], 'all')]);
end
