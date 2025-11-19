% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!!Provide session to be analysed!!!
session = 'Extinction'; % 'Extinction' or 'Renewal'

% !!! Provide trial type to analyse !!!
trial_type = 'looms'; % 'looms' or 'flashes'

% !!! Provide treatment mouse numbers for each treatment group !!!
received_psilocybin = [3; 5; 7; 8];
received_vehicle = [1; 2; 4; 6];

% !!! Provide filepath for where figures will be saved !!!
save_folder = 'D:\behaviours_results';

% !!! Decide whether to save figures or not !!!
do_save = false; % true or false
%%

behaviours = {'grooming', 'rearing', 'darting'};

% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

% Initialise 
grooming_folders = cell(length(mouse_files), 1);

num_mice = size(mouse_files, 1);
num_behaviours = numel(behaviours);

% Get behaviour folders for specified session
for mouse = 1:num_mice
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    grooming_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Grooming');
end

rearing_folders = cell(length(mouse_files), 1);

% Get behaviour folders for specified session
for mouse = 1:num_mice
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    rearing_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Rearing');
end


darting_folders = cell(length(mouse_files), 1);

% Get behaviour folders for specified session
for mouse = 1:num_mice
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    darting_folders{mouse} = fullfile(mouse_path, mouse_name,...
        session, 'Behavioural Data', 'Extracted Behaviours', 'Darting');
end

behaviour_data = cell(num_mice, num_behaviours);

for behaviour = 1:num_behaviours
    behaviour_folders = eval([behaviours{behaviour}, '_folders']);
    for mouse = 1:length(behaviour_folders)        
        current_behaviour_data_path = behaviour_folders{mouse};
        current_behaviour_data = dir(fullfile(current_behaviour_data_path, ...
            ['*_', behaviours{behaviour}, '_probabilities_' trial_type, '.mat']));
        file_to_load = fullfile(current_behaviour_data.folder, current_behaviour_data.name);
        if strcmp(trial_type, 'looms')
            load(file_to_load, 'looms_scores');
            behaviour_data{mouse, behaviour} = looms_scores;
        elseif strcmp(trial_type, 'flashes')
            load(file_to_load, 'flashes_scores');
            behaviour_data{mouse, behaviour} = flashes_scores;
        end
    end
end


%% Plotting

psi_behaviours = behaviour_data(received_psilocybin, :);
veh_behaviours = behaviour_data(received_vehicle, :);

timepoints = size(psi_behaviours{1,1}, 2); % assume same length for all
num_behaviours = size(behaviour_data, 2);

behaviour_labels = behaviours; % e.g. {'grooming','rearing','darting'}

for b = 1:num_behaviours
    % --- Psilocybin group ---
    psi_means = [];
    for m = 1:size(psi_behaviours,1)
        if isempty(psi_behaviours{m,b}), continue; end
        mouse_data = psi_behaviours{m,b}; % [trials × time]
        trial_mean = mean(mouse_data,1); % average across trials
        psi_means(end+1,:) = trial_mean; % append per mouse
    end
    psi_group_mean = mean(psi_means,1);
    psi_sem = std(psi_means,[],1) ./ sqrt(size(psi_means,1));

    % --- Vehicle group ---
    veh_means = [];
    for m = 1:size(veh_behaviours,1)
        if isempty(veh_behaviours{m,b}), continue; end
        mouse_data = veh_behaviours{m,b};
        trial_mean = mean(mouse_data,1);
        veh_means(end+1,:) = trial_mean;
    end
    veh_group_mean = mean(veh_means,1);
    veh_sem = std(veh_means,[],1) ./ sqrt(size(veh_means,1));

    % --- Plot for this behaviour ---
    figure; hold on;
    x = 1:timepoints;

    % Psilocybin shaded SEM
    fill([x fliplr(x)], [psi_group_mean+psi_sem fliplr(psi_group_mean-psi_sem)], ...
        'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

    % Vehicle shaded SEM
    fill([x fliplr(x)], [veh_group_mean+veh_sem fliplr(veh_group_mean-veh_sem)], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

    % Plot group means
    plot(psi_group_mean, 'r-', 'LineWidth', 2, 'DisplayName', 'Psilocybin');
    plot(veh_group_mean, 'b-', 'LineWidth', 2, 'DisplayName', 'Vehicle');

    % Vertical lines (not added to legend)
    xline(151, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
    xline(202, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');

    xlabel('Time (frame)');
    ylabel('Mean Probability');
    ylim([0 0.2])
    xlim([1 502])

    title_str = [session, ' ', trial_type, ' ', behaviour_labels{b}];
    title(title_str, 'Interpreter', 'none');
    legend('show');
    hold off;
    
    % Save figures 
    if do_save == true    
    save_path = fullfile(save_folder, [title_str, '.png']);
    exportgraphics(gcf, save_path, 'Resolution', 300);
    close; 
    end

end
