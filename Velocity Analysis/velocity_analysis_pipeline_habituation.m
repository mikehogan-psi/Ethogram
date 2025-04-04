% Specify the folder containing your .csv files
folderPath = 'D:\PhD 2nd Year\dlc_data_cohort2\Habituation';
% folderPath = 'D:\PhD 2nd Year\DLC Extinction Data  Cohort 1\Habituation';

% Get a list of all .csv files in the specified folder
fileList = dir(fullfile(folderPath, '*habituation*.csv'));

% identify number of mice
num_mice = length(fileList);

% Initialise a cell array to store each T matrix
T_matrices = cell(1, num_mice);
file_names = cell(num_mice, 1);

% iterate over each file, load it, and the correct columns from the .csv
% file
for i = 1:num_mice
    filePath = fullfile(folderPath, fileList(i).name);
     % load the .mat file
    data = readmatrix(filePath); 
    T_matrices{i} = data(:, [14, 15])'; 
    % for each mouse (cells) contain x and y coordinates (row) of body anterior body marker for each frame (colums)
    file_names{i} = fileList(i).name;
    % stores file names for sorting later
end

clearvars data i filePath % cleatr unnecessary variables 

%% Calculating freezing behaviour in % time spent freezing

% general freezing behaviour definition (same for all mice)

freezing_velocity_threshold = 1.1; % pixels per frame, change as needed (*0.048*15 to get cm/s value)
freezing_duration_threshold = 1*15; % convert seconds to frames (1 second × 15 FPS)

% initialise vector to store % time spent freezing for all mice
mice_freezing = zeros(1, num_mice);

% Looping through mice 
% - has to be done for each mouse individually (because all recorded for diffrent times i.e. different number of frames)
for mouse_idx = 1:num_mice
    % Step 1: Getting xy positional data and number of frames for this mouse
        
        x_pos = T_matrices{mouse_idx}(1,:);     % extracting x coordinates (for all frames)
        y_pos = T_matrices{mouse_idx}(2,:);     % extracting y coordinates (for all frames)
        num_frames = length(x_pos);             % extracting number frames (different for each mouse)
    
    % Step 2: calculating  x and y coordinate differences between frames 
        
        dx = [0, diff(x_pos)]; % difference in x coordimnates between frames
        dy = [0, diff(y_pos)]; % difference in y coordimnates between frames
                               % adds a 0 at start to maintain correct trial length
    
    % Step 3: calcuting Euclidian distance between frames    
    
        dist = sqrt(dx.^2 + dy.^2);     % calcuting Euclidian distance
          
    % Step 4: define freezing behaviour

        freeze_counter_matrix = dist < freezing_velocity_threshold;   % logical array that contains 1 dor when traveled more than treshold distance between frames (and 0 if not)
        [labeled_freezing, num_segments] = bwlabel(freeze_counter_matrix); % labeled_freezing is a vector same size as freeze_counter_matrix, but has a unique value for each 'blob' of consecutive trials
        
        segment_props = regionprops(labeled_freezing, 'Area');   % extracts length of connected regions in labeled_freezing and stores it in segment_props.Area
        segment_lengths = [segment_props.Area]; % Extract segment lengths from struct
        
        freezing_segments = segment_lengths(segment_lengths >= freezing_duration_threshold); % Filter segments that meet freezing duration threshold (15 frames or longer)
        total_freezing_frames = sum(freezing_segments); % calculates total number of frames spent in freezing
        freezing_seconds = total_freezing_frames / 15 ; % calculates total seconds spent freezing (at 15 frames per second sampling rate)

    % Step 5: calculate percentage of time spent freezing (in sec)
        total_seconds = num_frames / 15; % calculates the duration of total video in second
        percentage_freezing = freezing_seconds / total_seconds; % calulates  percentage of seconds spent freezing over entire recording

    % Step 6: save result
        mice_freezing(mouse_idx) = percentage_freezing; % stores result for each mouse 
end


%% Plotting results for vehicle and psilocybin group

% Splitting into treatment groups 
psilocybin_freezing_data = mice_freezing(2:2:end)*100;
vehicle_freezing_data = mice_freezing(1:2:end)*100;

% Calculate means and standard deviations
mean_psilocybin_habituation = mean(psilocybin_freezing_data);
std_psilocybin = std(psilocybin_freezing_data);
mean_vehicle_habituation = mean(vehicle_freezing_data);
std_vehicle = std(vehicle_freezing_data);

% Create a figure
figure;
hold on;

% Plot mean values as bars
bar(1, mean_psilocybin_habituation, 'FaceColor', 'red', 'FaceAlpha', 0.6, 'EdgeColor', 'black', 'LineWidth', 2, 'DisplayName', 'Psilocybin Mean');
bar(2, mean_vehicle_habituation, 'FaceColor', 'black', 'FaceAlpha', 0.6, 'EdgeColor', 'black', 'LineWidth', 2, 'DisplayName', 'Vehicle Mean');

% Add error bars for standard deviation
errorbar([1, 2], [mean_psilocybin_habituation, mean_vehicle_habituation], [std_psilocybin, std_vehicle], ...
    'k', 'linestyle', 'none', 'CapSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Std Dev');

% Customize plot
ylim ([0 30]);
xticks([1, 2]);
xticklabels({'Psilocybin', 'Vehicle'});

% Set axis properties
ax = gca;
ax.FontWeight = 'bold';  % Make all text in axes bold
ax.FontSize = 12;         % Increase font size for better visibility
ax.LineWidth = 2;         % Make the axes lines bold
ax.GridLineWidth = 0.5;   % Keep grid lines thin

% Labels and title
ylabel('Mean % Time Spent Freezing', 'FontWeight', 'bold');
% xlabel('Treatment Group', 'FontWeight', 'bold');
% title('Freezing Probability by Treatment Group', 'FontWeight', 'bold');

grid on;
hold off;

save_folder = 'D:\PhD 2nd Year\Figures\BNA poster figures'; 
save_path = fullfile(save_folder, 'pre_extinction_habituation.png');

% Save figure
print(gcf, save_path, '-dpng', '-r300');