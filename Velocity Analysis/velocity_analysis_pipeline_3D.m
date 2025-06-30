%% Velocity analysis from 3D data

%  sqrt(sum(diff(T_window(1:2,:), 1, 2).^2

% DLC_landmarks: 1. cable tip
%                2. back left courner (of implant)
%                3. back right courner (of implant)
%                4. nose
%                5. left ear
%                6. right ear
%                7. neck base
%                8. body anteriour
%                9. body posteriour
%                10. tail base
%                11. tail anteriour
%                12. tail posteriour
%                13. tail tip

%% Setup

% select data from which session shall be analysised
    sesh = 'extinction';
   % sesh = 'renewal';

% Setup directories   
    % folderpath to dlc data 
        %folderPath = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\data_all_mice\Extinction\data_from_DLC_iteration_3\DLC_data\';
    % folderpath to ssm fitted data 
        folderPath = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\data_all_mice\Extinction\data_from_DLC_iteration_3\SSM_fitted_data\';


    % directory to folder to save analysed data in
        save_dir = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\velocity_analysis_same_as_2D\';
    % directory to folder to save figures in
     %   fig_dir =  [ common_dir sesh '_analysis\figures\'];     


    video_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\data_all_mice\Extinction\data_from_DLC_iteration_3\video_data';


% Get a list of all .csv files in the specified folder
%fileList = dir(fullfile(folderPath, '*camera5_mouse*_extinction*'));
fileList = dir(fullfile(folderPath, ['mouse*_' sesh '*']));

% Initialise a cell arrays to store data information
T_matrices = cell(1, length(fileList)); % for each mouse (2 cells - part 1 and part 2) x and y coordinates (row) of body anterior body marker for each frame (colums)
file_names = cell(length(fileList), 1); % for each mouse (2 cells - part 1 and part 2)  filename


%% Extracting data (using DLC coordinates)

% iterate over each file, load it, and extract the T variable
for i = 1:length(fileList)
    filePath = fullfile(folderPath, fileList(i).name);
     % load the .mat file
    data = readmatrix(filePath); 
    T_matrices{i} = data(:, [23, 24])'; 
    % stores body_anterior x and y positional data in T_matrices array
    file_names{i} = fileList(i).name;
    % stores file names for sorting later
end

%% Delete 15 extra frames in mouse 1 extinxtion p2

T_matrices{2} = T_matrices{2}(:, [1:3514,3530:end]);

%% Extracting data (using SSM fitted data)

% iterate over each file, load it, and extract the T variable

for i = 1:length(fileList)
    filePath = fullfile(folderPath, fileList(i).name);

     % load T(Translation Vector) 
    load(filePath, 'Xfit'); 

    % extract x and y coordinates of anteriour body landmark
    T_matrices{i} = squeeze(Xfit(8,1:2,:));  

    % stores body_anterior x and y positional data in T_matrices array
    file_names{i} = fileList(i).name;
    % stores file names for sorting later
end

%% Sorting velocity matrix into mouse/part order

mouse_numbers = cellfun(@(x) sscanf(regexp(x, 'mouse(\d+)', 'match', 'once'), 'mouse%d'), file_names);
% extracting the mouse number from the file names

[~, sort_idx] = sort(mouse_numbers);
% sort by mouse number amd get index for sorted order

T_matrices = T_matrices(sort_idx); % sort matrix

disp('Sorted File Names:');
disp(file_names(sort_idx));

file_names = file_names(sort_idx); % sort filenames

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
num_mice = length(T_matrices); 
all_velocity_data = zeros(num_trials, num_frames, num_mice);

for mouse_idx = 1:num_mice
    x_pos = T_matrices{mouse_idx}(1,:);
    y_pos = T_matrices{mouse_idx}(2,:);

    dist = zeros(num_trials, num_frames);
    
    for trial_ix = 1:num_trials

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

% freezing_velocity_threshold = 1.1; % pixels per frame, change as needed (*0.048*15 to get cm/s value)
freezing_velocity_threshold = 0.053; % cm per frame (use this for xfit data because triangulation transforms to cm)

freezing_duration_threshold = 1*15; % convert seconds to frames (1 second × 15 FPS)


freeze_counter_matrix = all_velocity_data < freezing_velocity_threshold;

validated_freeze_matrix = zeros(size(freeze_counter_matrix));

for mouse_ix = 1:num_mice
    %iterating over each video
    for trials = 1:num_trials 
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
                validated_freeze_matrix(trials, labeled_freezing == segment, mouse_ix) = 1;
            end
        end
    end
end

% Saving Data 
%save([save_dir 'validated_freeze_matrix_xfit_' sesh], 'validated_freeze_matrix_xfit');

%clearvars -except validated_freeze_matrix common_dir save_dir fig_dir num_mice sesh 

%% Saving in same format as other behaviours (predicted by RFM)

freezing_mouse2 = validated_freeze_matrix(:,:,2);

predicted_labels_all = reshape(freezing_mouse2', [], 1);

predicted_labels = predicted_labels_all(1:10040);
save([save_dir 'xfit_053_freezing_mouse2_extinction_p1'], "predicted_labels")

predicted_labels = predicted_labels_all(10041:end);
save([save_dir 'xfit_053_freezing_mouse2_extinction_p2'], "predicted_labels")


%% STEP 8: Manual check predicted label accuracy and correct labels  

% choose a file you want to double-check the mouse behaviour of the frames predicted by the model to show the desired behaviour 
base_name = 'mouse2_extinction_p2';
predicted_labels_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\implanted_mice_analysis\velocity_analysis_same_as_2D';

    % Find the correct video and predicted labels file
    video_file_path = dir(fullfile(video_path, ['camera6_' base_name, '*.avi']));
    video_file_path = [video_path '\' video_file_path.name];
    predicted_file_path = dir(fullfile(predicted_labels_path, ['xfit_053_freezing_', base_name, '*.mat']));
    predicted_file_path = [predicted_labels_path '\' predicted_file_path.name];
    
    % fun validate predcitions function
         validate_predictions(base_name, video_file_path, predicted_file_path, save_dir);
         % should open video in which can manually scrub through frames and see if the predicted labels match actual behaviour and if necessary correct it
    
    % repeat this for a couple of videos to generate enough verified_labels data to update model (to increase its accuracy)     


%% Splitting into treatment groups, stim types and trial periods

psilocybin_freezing_data = validated_freeze_matri x(:, :, 2:2:end);
vehicle_freezing_data = validated_freeze_matrix(:, :, 1:2:end);

% 2 stimsets -> each 20 stimuli (0 = flash, 1 = loom) -> each trial
stim_set_1 = [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1,...
    0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1];

stim_set_2 = [0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0,...
    1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0 ,0 ,1 ,0];

% defining pre and post stimulus frames
pre_stim_period = 2:152;
during_stim_period = 153:201;
post_stim_period = 202:502;
early_post_stim_period = 203:352;
late_post_stim_period =  353:502;

% setting indexes for looms and flashes
stim_set_1_looms_idx = logical(stim_set_1);
stim_set_1_flashes_idx = ~stim_set_1_looms_idx;
stim_set_2_looms_idx = logical(stim_set_2);
stim_set_2_flashes_idx = ~stim_set_2_looms_idx;

% sepparating mice that get stimset 1 or stimset 2 to extraxt looms or flashes
psi_loom_ss1 = psilocybin_freezing_data(stim_set_1_looms_idx, :, 1:2:end);      % gets all the loom trials for stimset 1 mice (uneven numbers)
psi_flash_ss1 = psilocybin_freezing_data(stim_set_1_flashes_idx, :, 1:2:end);   % gets all the flash trials for stimset 1
psi_loom_ss2 = psilocybin_freezing_data(stim_set_2_looms_idx, :, 2:2:end);
psi_flash_ss2 = psilocybin_freezing_data(stim_set_2_flashes_idx, :, 2:2:end);

veh_loom_ss1 = vehicle_freezing_data(stim_set_1_looms_idx, :, 1:2:end);
veh_flash_ss1 = vehicle_freezing_data(stim_set_1_flashes_idx, :, 1:2:end);
veh_loom_ss2 = vehicle_freezing_data(stim_set_2_looms_idx, :, 2:2:end);
veh_flash_ss2 = vehicle_freezing_data(stim_set_2_flashes_idx, :, 2:2:end);

% concatenating all psi and veh freezing data
psi_loom_data = interleave(psi_loom_ss1, psi_loom_ss2);      % use interleave to preserve order 
psi_flash_data = interleave(psi_flash_ss1, psi_flash_ss2);
veh_loom_data = interleave(veh_loom_ss1, veh_loom_ss2);
veh_flash_data = interleave(veh_flash_ss1, veh_flash_ss2);

%extracting pre and post-stim freezing data
psi_post_stim_loom = psi_loom_data(:, post_stim_period, :);
psi_post_stim_flash = psi_flash_data(:, post_stim_period, :);
veh_post_stim_loom = veh_loom_data(:, post_stim_period, :);
veh_post_stim_flash = veh_flash_data(:, post_stim_period, :);

psi_pre_stim_loom = psi_loom_data(:, pre_stim_period, :);
psi_pre_stim_flash = psi_flash_data(:, pre_stim_period, :);
veh_pre_stim_loom = veh_loom_data(:, pre_stim_period, :);
veh_pre_stim_flash = veh_flash_data(:, pre_stim_period, :);

%extracting early and late post-stim freezing data

psi_early_post_stim_loom = psi_loom_data(:, early_post_stim_period, :);
psi_early_post_stim_flash = psi_flash_data(:, early_post_stim_period, :);
veh_early_post_stim_loom = veh_loom_data(:, early_post_stim_period, :);
veh_early_post_stim_flash = veh_flash_data(:, early_post_stim_period, :);

psi_late_post_stim_loom = psi_loom_data(:, late_post_stim_period, :);
psi_late_post_stim_flash = psi_flash_data(:, late_post_stim_period, :);
veh_late_post_stim_loom = veh_loom_data(:, late_post_stim_period, :);
veh_late_post_stim_flash = veh_flash_data(:, late_post_stim_period, :);


% %% Checking positions in matrixes correct for figure lableling
% 
% % checking if interleave function restored originial order of mice example mice 15 and 20
% correct_log1 = validated_freeze_matrix(stim_set_2_looms_idx, :, 15) == veh_loom_data (:,:,8);
% correct_log2 = validated_freeze_matrix(stim_set_2_flashes_idx, :, 15) == veh_flash_data (:,:,8);
% correct_log3 = validated_freeze_matrix(stim_set_2_looms_idx, :, 20) == psi_loom_data (:,:,10);
% correct_log4 = validated_freeze_matrix(stim_set_2_flashes_idx, :, 20) == psi_flash_data (:,:,10);
% 
% clearvars correct_log1 correct_log2 correct_log3 correct_log4
% 

%% plotting mean freezing after stimulus

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
fig1 = figure;
hold on;
errorbar(trials, psi_loom_trial_means, psi_loom_trial_sems, '-o', 'DisplayName', 'Loom');
errorbar(trials, psi_flash_trial_means, psi_flash_trial_sems, '-x', 'DisplayName', 'Flash');
hold off;
title('Freezing Behaviour Across Repeated Trials: Psilocybin');
legend('Location', 'Best');
grid on;

% Plot Vehicle Group
fig2 = figure;
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
    ylim([0, max([psi_loom_trial_means + psi_loom_trial_sems, ...
                  psi_flash_trial_means + psi_flash_trial_sems, ...
                  veh_loom_trial_means + veh_loom_trial_sems, ...
                  veh_flash_trial_means + veh_flash_trial_sems], [], 'all')]);
end

% save figures
saveas(fig1, [ fig_dir sesh 'freezing_over_trials_psi.fig']);
saveas(fig2, [ fig_dir sesh 'freezing_over_trials_veh.fig']);

saveas(fig1, [ fig_dir sesh 'freezing_over_trials_psi.jpg']);
saveas(fig2, [ fig_dir sesh 'freezing_over_trials_veh.jpg']);

%% plotting group variation in freezing (post-stim period only)

%Plotting median of mean freezing probability between all mice over all
%loom trials (post-stim period only)
mean_freezing_psi_loom_post = squeeze(mean(psi_post_stim_loom, 1));
mean_freezing_veh_loom_post = squeeze(mean(veh_post_stim_loom, 1));
time_seconds = [1:301] / 15; 

fig3 = figure;
plot(time_seconds, median(mean_freezing_psi_loom_post, 2)')
hold on 
plot(time_seconds, median(mean_freezing_veh_loom_post, 2)')
xlabel('time in sec');
ylabel('Median Freezing Probability')
legend('Psilocybin', 'Vehicle');
title('Freezing Per Frame, All Mice (Post-Loom Period Only)')
xlim([0,20]);
xticks(time_seconds(30:30:end))
ylim([0, 1]);
hold off


%Plotting median of mean freezing probability between all mice over all
%flash trials (post-stim period only)
mean_freezing_psi_flash_post = squeeze(mean(psi_post_stim_flash, 1));
mean_freezing_veh_flash_post = squeeze(mean(veh_post_stim_flash, 1));

fig4 = figure;
plot(time_seconds,median(mean_freezing_psi_flash_post, 2))
hold on 
plot(time_seconds,median(mean_freezing_veh_flash_post, 2))
xlabel('time in sec');
ylabel('Median Freezing Probability')
legend('Psilocybin', 'Vehicle');
title('Freezing Per Frame, All Mice (Post-Flash Period Only)')
xlim([1,20]);
xticks(time_seconds(30:30:end))
ylim([0, 1]);
hold off

save figures
saveas(fig3, [ fig_dir sesh 'freezing_over_time_loom.fig']);
saveas(fig4, [ fig_dir sesh 'freezing_over_time_flash.fig']);

saveas(fig3, [ fig_dir sesh 'freezing_over_time_loom.jpg']);
saveas(fig4, [ fig_dir sesh 'freezing_over_time_flash.jpg']);

% %% plotting group variation in freezing (over entire trials period) 
% 
% %Plotting median of mean freezing probability between all mice over all loom trials (over whole trial period)
% mean_freezing_psi = squeeze(mean(psi_loom_data, 1));
% mean_freezing_veh = squeeze(mean(veh_loom_data , 1));
% time_seconds = [1:301] / 15; 
% 
% fig3 = figure;
% 
% plot(time_seconds = [1:301] / 15; , median(mean_freezing_psi, 2)', 'LineWidth', 1.5)
% hold on 
% plot(frames, median(mean_freezing_veh, 2)', 'LineWidth', 1.5)
% xline(152, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus Start'); % Stimulus start
% xline(202, '-k', 'LineWidth', 1.5, 'DisplayName', 'Stimulus End');   % Stimulus end
% 
% xlabel('time in sec');
% ylabel('Median Freezing Probability')
% legend('Psilocybin', 'Vehicle');
% title('Freezing Per Frame (to looms), All Mice')
% xlim([1,502]);
% %xticks(time_seconds(1:30:end))
% ylim([0, 1]);
% hold off

%% plotting group variation in freezing (over entire trials period) for early and later trials

% Split data into early and later trials
early_trials_psi_loom = psi_loom_data(1:10, :, :);
early_trials_veh_loom = veh_loom_data(1:10, :, :);
late_trials_psi_loom = psi_loom_data(11:20, :, :);
late_trials_veh_loom = veh_loom_data(11:20, :, :);

%computing median of mean freezing probability between all mice over all loom trials (over whole trial period)
early_mean_freezing_psi = squeeze(mean(early_trials_psi_loom, 1));
early_mean_freezing_veh = squeeze(mean(early_trials_veh_loom, 1));
late_mean_freezing_psi = squeeze(mean(late_trials_psi_loom, 1));
late_mean_freezing_veh = squeeze(mean(late_trials_veh_loom, 1));

% Constants
framerate = 15; % frames per second
frames = 1:502; % original frame numbers
time_seconds = (frames / framerate) - (202 / framerate); % Convert frames to seconds and set Frame 202 as time 0

% Create figure
fig3 = figure;

% Plot median freezing probabilities

hold on
plot(time_seconds, median(early_mean_freezing_psi, 2)', 'Color', [0.85, 0.4, 0], 'LineWidth', 1, 'DisplayName', 'Early Psilocybin'); % Orange
hold on;
plot(time_seconds, median(early_mean_freezing_veh, 2)', 'Color', [0.4, 0.8, 1], 'LineWidth', 1, 'DisplayName', 'Early Vehicle'); % Light Blue
plot(time_seconds, median(late_mean_freezing_psi, 2)', 'Color', [1, 0, 0], 'LineWidth', 1, 'DisplayName', 'Late Psilocybin'); % Red
plot(time_seconds, median(late_mean_freezing_veh, 2)', 'Color', [0, 0, 0.8], 'LineWidth', 1, 'DisplayName', 'Late Vehicle'); % Dark Blue

% Add stimulus start and end lines
xline(0, '-k', 'LineWidth', 1, 'DisplayName', 'Stimulus End'); % Time 0 (Frame 202)
xline((152 / framerate) - (202 / framerate), '-k', 'LineWidth', 1, 'DisplayName', 'Stimulus Start'); % Time of Frame 152

% Customize plot
xlabel('Time (sec)');
ylabel('Median Freezing Probability');
legend('Early Psilocybin', 'Early Vehicle', 'Late Psilocybin', 'Late Vehicle',   'Location', 'northeast');
title('Freezing Per Frame to loom stimuli (all Mice)');
xlim([(1 / framerate) - (202 / framerate), (502 / framerate) - (202 / framerate)]); % Adjust x-axis to time range
ylim([0, 1]);
grid on; % Optional: Add gridlines

hold off;


%% Plotting group variation in freezing (over entire trials period)

% Constants
framerate = 15; % frames per second
frames = 1:502; % original frame numbers
time_seconds = (frames / framerate) - (202 / framerate); % Convert frames to seconds and set Frame 202 as time 0

% Create figure
fig3 = figure;

% Plot median freezing probabilities
plot(time_seconds, median(mean_freezing_psi, 2)', 'b', 'LineWidth', 1.5, 'DisplayName', 'Psilocybin');
hold on;
plot(time_seconds, median(mean_freezing_veh, 2)', 'r', 'LineWidth', 1.5, 'DisplayName', 'Vehicle');

% Add stimulus start and end lines
xline(0, '-k', 'LineWidth', 1, 'DisplayName', 'Stimulus End'); % Time 0 (Frame 202)
xline((152 / framerate) - (202 / framerate), '-k', 'LineWidth', 1, 'DisplayName', 'Stimulus Start'); % Time of Frame 152

% Customize plot
xlabel('Time (sec)');
ylabel('Median Freezing Probability');
legend('Psilocybin', 'Vehicle', 'Location', 'northeast');
title('Freezing Per Frame to loom stimuli (all Mice)');
xlim([(1 / framerate) - (202 / framerate), (502 / framerate) - (202 / framerate)]); % Adjust x-axis to time range
ylim([0, 1]);
grid on; % Optional: Add gridlines

hold off;


%% Plotting individual mean freezing probability between all mice over all loom trials (post-stim period only)

figure;
% Using imagesc to plot the data
imagesc([mean_freezing_psi_loom_post'; mean_freezing_veh_loom_post']); 
colorbar; % Add colour bar

% Add titles and labels
title('Freezing Per Frame Per Mouse (Post-Loom Period Only)'); % Graph title
xlabel('Frame'); % x-axis label
ylabel('Mice (Psilocybin Group and Vehicle Group)'); % y-axis label

% Label for the colour bar
c = colorbar;
c.Label.String = 'Mean Freezing Probability';

% Adjust the axis for a better fit
axis tight; % Fit the axis tightly around the data



%% Plotting mean freezing probability between all mice over all loom trials
% (all frames)
mean_freezing_psi = squeeze(mean(psi_loom_data, 1));
mean_freezing_veh = squeeze(mean(veh_loom_data, 1));

figure;
% Using imagesc to plot the data
imagesc([mean_freezing_psi'; mean_freezing_veh']); 
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

%% Comparing mean freezing (across mice and trials) in early and late post-stim periods

fps = 15; % defining frame rate

% Calculating mean freezing time (in seconds) and SEM (across all mice and all trials) 

% for early post-stim period 
early_psi_loom_trial_mean = mean(mean(sum(psi_early_post_stim_loom, 2) / fps, 3));
early_psi_flash_trial_mean = mean(mean(sum(psi_early_post_stim_flash, 2) / fps, 3));
early_veh_loom_trial_mean = mean(mean(sum(veh_early_post_stim_loom, 2) / fps, 3));
early_veh_flash_trial_mean = mean(mean(sum(veh_early_post_stim_flash, 2) / fps, 3));

early_psi_loom_trial_sem = mean(std(sum(psi_early_post_stim_loom, 2) / fps, [], 3) / sqrt(size(psi_early_post_stim_loom, 3)));
early_psi_flash_trial_sem = mean(std(sum(psi_early_post_stim_flash, 2) / fps, [], 3) / sqrt(size(psi_early_post_stim_flash, 3)));
early_veh_loom_trial_sem = mean(std(sum(veh_early_post_stim_loom, 2) / fps, [], 3) / sqrt(size(veh_early_post_stim_loom, 3)));
early_veh_flash_trial_sem = mean(std(sum(veh_early_post_stim_flash, 2) / fps, [], 3) / sqrt(size(veh_early_post_stim_flash, 3)));

% for late post-stim period 
late_psi_loom_trial_mean = mean(mean(sum(psi_late_post_stim_loom, 2) / fps, 3));
late_psi_flash_trial_mean = mean(mean(sum(psi_late_post_stim_flash, 2) / fps, 3));
late_veh_loom_trial_mean = mean(mean(sum(veh_late_post_stim_loom, 2) / fps, 3));
late_veh_flash_trial_mean = mean(mean(sum(veh_late_post_stim_flash, 2) / fps, 3));

late_psi_loom_trial_sem = mean(std(sum(psi_late_post_stim_loom, 2) / fps, [], 3) / sqrt(size(psi_late_post_stim_loom, 3)));
late_psi_flash_trial_sem = mean(std(sum(psi_late_post_stim_flash, 2) / fps, [], 3) / sqrt(size(psi_late_post_stim_flash, 3)));
late_veh_loom_trial_sem = mean(std(sum(veh_late_post_stim_loom, 2) / fps, [], 3) / sqrt(size(veh_late_post_stim_loom, 3)));
late_veh_flash_trial_sem = mean(std(sum(veh_late_post_stim_flash, 2) / fps, [], 3) / sqrt(size(veh_late_post_stim_flash, 3)));

% concatinate early and late data
psi_loom_means = [early_psi_loom_trial_mean, late_psi_loom_trial_mean];
psi_flash_means = [early_psi_flash_trial_mean, late_psi_flash_trial_mean];
veh_loom_means = [early_veh_loom_trial_mean, late_veh_loom_trial_mean];
veh_flash_means = [early_veh_flash_trial_mean, late_veh_flash_trial_mean];

psi_loom_sems = [early_psi_loom_trial_sem, late_psi_loom_trial_sem];
psi_flash_sems = [early_psi_flash_trial_sem, late_psi_flash_trial_sem];
veh_loom_sems = [early_veh_loom_trial_sem, late_veh_loom_trial_sem];
veh_flash_sems = [early_veh_flash_trial_sem, late_veh_flash_trial_sem];

%plot
figure;
hold on;
errorbar(1:2,psi_loom_means, psi_loom_sems, '-o', 'DisplayName', 'Psi Loom','Color','red');
errorbar(1:2,psi_flash_means, psi_flash_sems, '-x', 'DisplayName', 'Psi Flash', 'Color','red');
errorbar(1:2,veh_loom_means, veh_loom_sems, '-o', 'DisplayName', 'Veh Loom', 'Color','blue');
errorbar(1:2,veh_flash_means, veh_flash_sems, '-x', 'DisplayName', 'Veh Flash', 'Color','blue');

ylabel('Mean freezing time in s (across all mice and all trials)');
xticks([1,2]);
xticklabels({'frames 202-350', 'frames 351-502'})

xlim([0,3]);
ylim([0, 10]);


hold off;
title('Freezing Behaviour in early and late post-stim period');
legend('Location', 'Best');
grid on;

%% Comparing mean difference (across all mice) in early and late post-stim period for individual trials

% difference between early and late post-stim periods across trials 
diff_psi_loom_means =  mean(sum(psi_early_post_stim_loom, 2) / fps, 3) - mean(sum(psi_late_post_stim_loom, 2) / fps, 3);
diff_psi_flash_means = mean(sum(psi_early_post_stim_flash, 2) / fps, 3) - mean(sum(psi_late_post_stim_flash, 2) / fps, 3);
diff_veh_loom_means =   mean(sum(veh_early_post_stim_loom, 2) / fps, 3) - mean(sum(veh_late_post_stim_loom, 2) / fps, 3);
diff_veh_flash_means =  mean(sum(veh_early_post_stim_flash, 2) / fps, 3) - mean(sum(veh_late_post_stim_flash, 2) / fps, 3);

diff_psi_loom_trial_sems = std(sum(diff_psi_loom_means, 2) / fps, [], 3) / sqrt(size(diff_psi_loom_means, 3));
diff_psi_flash_trial_sems = std(sum(diff_psi_flash_means, 2) / fps, [], 3) / sqrt(size(diff_psi_flash_means, 3));
diff_veh_loom_trial_sems = std(sum(diff_veh_loom_means, 2) / fps, [], 3) / sqrt(size(diff_veh_loom_means, 3));
diff_veh_flash_trial_sems = std(sum(diff_veh_flash_means, 2) / fps, [], 3) / sqrt(size(diff_veh_flash_means, 3));

trials = 1:20;
% trial numbers

% Plot Psilocybin Group
figure;
hold on;
errorbar(trials, diff_psi_loom_means, diff_psi_loom_trial_sems, '-o', 'DisplayName', 'Loom');
errorbar(trials, diff_psi_flash_means,diff_psi_flash_trial_sems, '-x', 'DisplayName', 'Flash');
hold off;
title('difference in Freezing between early and late post-stim period: Psilocybin');
legend('Location', 'Best');
grid on;

% Plot Vehicle Group
figure;
hold on;
errorbar(trials, diff_veh_loom_means, diff_veh_loom_trial_sems, '-o', 'DisplayName', 'Loom');
errorbar(trials, diff_veh_flash_means , diff_veh_flash_trial_sems, '-x', 'DisplayName', 'Flash');
hold off;
title('difference in Freezing between early and late post-stim period: Vehicle');
legend('Location', 'Best');
grid on;

%% Comparing mean freezing (of all mice) in early and late post-stim periods across individual trials

fps = 15; % defining frame rate
trials = 1:20;% trial numbers

% Calculating mean freezing time (in seconds) and SEM 

% for early post-stim period 
early_psi_loom_trial_mean = mean(sum(psi_early_post_stim_loom, 2) / fps, 3);
early_psi_flash_trial_mean = mean(sum(psi_early_post_stim_flash, 2) / fps, 3);
early_veh_loom_trial_mean = mean(sum(veh_early_post_stim_loom, 2) / fps, 3);
early_veh_flash_trial_mean = mean(sum(veh_early_post_stim_flash, 2) / fps, 3);

early_psi_loom_trial_sem = std(sum(psi_early_post_stim_loom, 2) / fps, [], 3) / sqrt(size(psi_early_post_stim_loom, 3));
early_psi_flash_trial_sem = std(sum(psi_early_post_stim_flash, 2) / fps, [], 3) / sqrt(size(psi_early_post_stim_flash, 3));
early_veh_loom_trial_sem = std(sum(veh_early_post_stim_loom, 2) / fps, [], 3) / sqrt(size(veh_early_post_stim_loom, 3));
early_veh_flash_trial_sem = std(sum(veh_early_post_stim_flash, 2) / fps, [], 3) / sqrt(size(veh_early_post_stim_flash, 3));

% for late post-stim period 
late_psi_loom_trial_mean = mean(sum(psi_late_post_stim_loom, 2) / fps, 3);
late_psi_flash_trial_mean = mean(sum(psi_late_post_stim_flash, 2) / fps, 3);
late_veh_loom_trial_mean = mean(sum(veh_late_post_stim_loom, 2) / fps, 3);
late_veh_flash_trial_mean = mean(sum(veh_late_post_stim_flash, 2) / fps, 3);

late_psi_loom_trial_sem = std(sum(psi_late_post_stim_loom, 2) / fps, [], 3) / sqrt(size(psi_late_post_stim_loom, 3));
late_psi_flash_trial_sem = std(sum(psi_late_post_stim_flash, 2) / fps, [], 3) / sqrt(size(psi_late_post_stim_flash, 3));
late_veh_loom_trial_sem = std(sum(veh_late_post_stim_loom, 2) / fps, [], 3) / sqrt(size(veh_late_post_stim_loom, 3));
late_veh_flash_trial_sem = std(sum(veh_late_post_stim_flash, 2) / fps, [], 3) / sqrt(size(veh_late_post_stim_flash, 3));


%plot for looms
figure;

hold on;
errorbar(trials,early_psi_loom_trial_mean, early_psi_loom_trial_sem, '-o', 'DisplayName', 'Psi Loom early','Color','red');
errorbar(trials,late_psi_loom_trial_mean, late_psi_loom_trial_sem, '-x', 'DisplayName', 'Psi loom late', 'Color','red');
errorbar(trials,early_veh_loom_trial_mean, early_veh_loom_trial_sem, '-o', 'DisplayName', 'Veh Loom early', 'Color','blue');
errorbar(trials,late_veh_loom_trial_mean, late_veh_loom_trial_sem, '-x', 'DisplayName', 'Veh Loom late', 'Color','blue');
ylabel('Mean freezing time in s');
hold off;

title('Freezing Behaviour to Loom stimuli in early and late post-stim period');
legend('Location', 'Best');
grid on;

%plot for flash
figure;

hold on;
errorbar(trials,early_psi_flash_trial_mean, early_psi_flash_trial_sem, '-o', 'DisplayName', 'Psi Flash early','Color','red');
errorbar(trials,late_psi_flash_trial_mean, late_psi_flash_trial_sem, '-x', 'DisplayName', 'Psi Flash late', 'Color','red');
errorbar(trials,early_veh_flash_trial_mean, early_veh_flash_trial_sem, '-o', 'DisplayName', 'Veh Flash early', 'Color','blue');
errorbar(trials,late_veh_flash_trial_mean, late_veh_flash_trial_sem, '-x', 'DisplayName', 'Veh Flash late', 'Color','blue');
ylabel('Mean freezing time in s');
hold off;

title('Freezing Behaviour to flash stimuli in early and late post-stim period');
legend('Location', 'Best');
grid on;
