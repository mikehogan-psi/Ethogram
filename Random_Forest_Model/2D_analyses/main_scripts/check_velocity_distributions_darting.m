%% Plot velocity distributions
% Folder containing the .mat files
dataDir = 'Z:\Abi\behavioral_analysis\2D_behavioural analysis\renewal_analysis\darting_analyses\RFM_darting_17\predicted_labels_data';  % <-- Change to your actual folder path
matFiles = dir(fullfile(dataDir, '*.mat'));

darting_vel = [];      % Velocities when darting (label == 1)
non_darting_vel = [];  % Velocities when NOT darting (label == 0)

% Loop through each file and extract velocities
for i = 1:length(matFiles)
    data = load(fullfile(dataDir, matFiles(i).name));
    
    if isfield(data, 'predicted_labels') && isfield(data, 'av_velocity')
        labels = data.predicted_labels_TH(:);
        velocities = data.av_velocity(:);
        
        validIdx = ~isnan(labels) & ~isnan(velocities);
        labels = labels(validIdx);
        velocities = velocities(validIdx);

        darting_vel = [darting_vel; velocities(labels == 1)];
        non_darting_vel = [non_darting_vel; velocities(labels == 0)];
    else
        warning('Missing variables in %s', matFiles(i).name);
    end
end

% Plot histograms
figure;
hold on;

% Histogram for non-darting
histogram(non_darting_vel, ... % 'Normalization', 'probability', ...
    'BinWidth', 0.5, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Histogram for darting
histogram(darting_vel, ... % 'Normalization', 'probability', ...
    'BinWidth', 0.5, 'FaceColor', [1 0.2 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

xlabel('Average Velocity');
ylabel('Probability');
title('Velocity Distributions: Renewal - threshold filters applied');
legend({'Non-Darting', 'Darting'});
grid on;

%% Plot log-transformed velocity distributions
% Folder containing the .mat files
dataDir = 'Z:\Abi\behavioral_analysis\2D_behavioural analysis\extinction_analysis\darting_analyses\RFM_darting_17\predicted_labels_data';  % <-- Change to your actual folder path
matFiles = dir(fullfile(dataDir, '*.mat'));

darting_vel = [];      % Velocities when darting (label == 1)
non_darting_vel = [];  % Velocities when NOT darting (label == 0)

% Loop through each file and extract velocities
for i = 1:length(matFiles)
    data = load(fullfile(dataDir, matFiles(i).name));

    if isfield(data, 'predicted_labels') && isfield(data, 'av_velocity')
        labels = data.predicted_labels_TH(:);
        velocities = data.av_velocity(:);

        validIdx = ~isnan(labels) & ~isnan(velocities) & velocities > 0; % Exclude zero/negative values for log
        labels = labels(validIdx);
        velocities = velocities(validIdx);

        darting_vel = [darting_vel; velocities(labels == 1)];
        non_darting_vel = [non_darting_vel; velocities(labels == 0)];
    else
        warning('Missing variables in %s', matFiles(i).name);
    end
end

% Apply log10 transformation
log_darting_vel = log10(darting_vel);
log_non_darting_vel = log10(non_darting_vel);

% Plot histograms
figure;
hold on;

% Histogram for non-darting
histogram(log_non_darting_vel, 'Normalization', 'probability', ...
    'BinWidth', 0.1, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Histogram for darting
histogram(log_darting_vel, 'Normalization', 'probability', ...
    'BinWidth', 0.1, 'FaceColor', [1 0.2 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

xlabel('log_{10}(Average Velocity)');
ylabel('Probability');
title('Log-Velocity Distributions: Extinction - thresholds applied');
legend({'Non-Darting', 'Darting'});
grid on;

%% calculate average darting time

% Set your directory
dataDir = 'Z:\Abi\behavioral_analysis\2D_behavioural analysis\extinction_analysis\darting_analyses\RFM_darting_13\predicted_labels_data';  % <-- Update this
matFiles = dir(fullfile(dataDir, '*.mat'));

all_bout_lengths = [];  % Collect bout lengths across all mice
per_mouse_avg = zeros(length(matFiles), 1);  % Store average per mouse

for i = 1:length(matFiles)
    filePath = fullfile(dataDir, matFiles(i).name);
    data = load(filePath);

    if isfield(data, 'predicted_labels_TH')
        labels = data.predicted_labels_TH(:);

        % Use bwlabel to find continuous darting bouts
        darting_runs = bwlabel(labels == 1);
        num_bouts = max(darting_runs);

        if num_bouts > 0
            bout_lengths = zeros(num_bouts, 1);
            for j = 1:num_bouts
                bout_lengths(j) = sum(darting_runs == j);
            end

            % Store all bouts for overall stats
            all_bout_lengths = [all_bout_lengths; bout_lengths];

            % Store per-mouse average
            per_mouse_avg(i) = mean(bout_lengths);

            fprintf('Mouse %02d (%s): %d bouts, avg %.2f frames\n', ...
                i, matFiles(i).name, num_bouts, per_mouse_avg(i));
        else
            fprintf('Mouse %02d (%s): No darting bouts found.\n', ...
                i, matFiles(i).name);
            per_mouse_avg(i) = NaN;
        end
    else
        warning('Missing predicted_labels in %s', matFiles(i).name);
        per_mouse_avg(i) = NaN;
    end
end

% Overall average darting bout duration
overall_avg = mean(all_bout_lengths);
fprintf('\nOverall average darting bout duration: %.2f frames\n', overall_avg);
