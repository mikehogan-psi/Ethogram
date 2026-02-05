% !!!Provide master directory with all data in!!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

% !!!Provide session to be analysed!!!
session = 'Extinction'; % 'Extinction' or 'Renewal'

% !!! Provide trial type to analyse !!!
trial_type = 'looms'; % 'looms' or 'flashes'

% !!! Provide treatment mouse numbers for each treatment group !!!
received_psilocybin = [3; 5; 7; 8];
received_vehicle = [1; 2; 4; 6; 9];

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

%% ===== Summary bar plots (mean across mice ± SEM) + individual points + stats =====

% Windows (based on your xlines)
stim_start = 151;
stim_end   = 202;

post_start = 203;
post_end   = timepoints;  % to end

% Sanity clamp (in case timepoints < 202 etc.)
stim_start = max(1, min(stim_start, timepoints));
stim_end   = max(1, min(stim_end,   timepoints));
post_start = max(1, min(post_start, timepoints));
post_end   = max(1, min(post_end,   timepoints));

stim_idx = stim_start:stim_end;
post_idx = post_start:post_end;

% Preallocate storage: per mouse value in the relevant window
% We'll store two columns: [psi_values, veh_values] in cells per behaviour
psi_window_vals = cell(num_behaviours, 1);
veh_window_vals = cell(num_behaviours, 1);

% Decide which window per behaviour: darting -> stim, others -> post
is_darting = strcmpi(behaviours, 'darting');

for b = 1:num_behaviours
    if is_darting(b)
        idx = stim_idx;
    else
        idx = post_idx;
    end

    % --- Psilocybin group per-mouse window means ---
    psi_vals = [];
    for m = 1:size(psi_behaviours,1)
        if isempty(psi_behaviours{m,b}), continue; end
        mouse_data = psi_behaviours{m,b};     % [trials x time]
        per_mouse_timecourse = mean(mouse_data, 1);  % mean across trials -> [1 x time]
        psi_vals(end+1,1) = mean(per_mouse_timecourse(idx)); % scalar per mouse
    end

    % --- Vehicle group per-mouse window means ---
    veh_vals = [];
    for m = 1:size(veh_behaviours,1)
        if isempty(veh_behaviours{m,b}), continue; end
        mouse_data = veh_behaviours{m,b};
        per_mouse_timecourse = mean(mouse_data, 1);
        veh_vals(end+1,1) = mean(per_mouse_timecourse(idx));
    end

    psi_window_vals{b} = psi_vals;
    veh_window_vals{b} = veh_vals;
end

% --- Bar plots + stats per behaviour ---
for b = 1:num_behaviours
    psi_vals = psi_window_vals{b};
    veh_vals = veh_window_vals{b};

    % Skip if missing data
    if isempty(psi_vals) || isempty(veh_vals)
        warning('Skipping %s: missing data in one group.', behaviours{b});
        continue;
    end

    psi_mean = mean(psi_vals);
    veh_mean = mean(veh_vals);

    psi_sem = std(psi_vals, 0) ./ sqrt(numel(psi_vals));
    veh_sem = std(veh_vals, 0) ./ sqrt(numel(veh_vals));

    % Stats: ranksum (robust) + ttest2 (optional)
    p_ranksum = ranksum(psi_vals, veh_vals);
    [~, p_ttest] = ttest2(psi_vals, veh_vals);

    % Build figure
    figure; hold on;

    means = [veh_mean, psi_mean];
    sems  = [veh_sem,  psi_sem];

    bh = bar(1:2, means, 'FaceAlpha', 0.8);
    bh.FaceColor = 'flat';
    bh.CData(1,:) = [0 0 1];
    bh.CData(2,:) = [1 0 0];

    errorbar(1:2, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

    % Overlay individual datapoints with jitter
    jitter_amt = 0.08;
    xVeh = 1 + (rand(size(veh_vals)) - 0.5) * 2*jitter_amt;
    xPsi = 2 + (rand(size(psi_vals)) - 0.5) * 2*jitter_amt;

    scatter(xVeh, veh_vals, 60, 'k', 'filled', 'MarkerFaceAlpha', 0.75);
    scatter(xPsi, psi_vals, 60, 'k', 'filled', 'MarkerFaceAlpha', 0.75);

    % Formatting
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Vehicle','Psilocybin'});
    ylabel('Mean Probability');
    title_str = [session,' ',trial_type,' ',behaviours{b}];

    if is_darting(b)
        win_str = sprintf('Stim window: frames %d-%d', stim_start, stim_end);
    else
        win_str = sprintf('Post-stim window: frames %d-%d', post_start, post_end);
    end
    title({title_str; win_str}, 'Interpreter','none');

    % Nice y-limits
    ymin = 0;
    ymax = max([veh_vals; psi_vals]) * 1.25;
    if ymax <= 0, ymax = 0.1; end
    ylim([ymin ymax]);

    grid on; box off;

    % Significance line + stars (ranksum)
    yStar = ymax * 0.92;
    plot([1 2], [yStar yStar], 'k-', 'LineWidth', 1.5);
    if p_ranksum < 0.001
        star = '***';
    elseif p_ranksum < 0.01
        star = '**';
    elseif p_ranksum < 0.05
        star = '*';
    else
        star = 'n.s.';
    end
    text(1.5, yStar*1.02, star, 'HorizontalAlignment','center', 'FontSize', 14);

    % ---- ADD P-VALUE TEXT ON THE GRAPH ----
    p_text = sprintf('ranksum p = %.3g\nttest2 p = %.3g', p_ranksum, p_ttest);

    % Put it top-left in axes coordinates (always in same spot)
    text(0.05, 0.95, p_text, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 11, ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.3 0.3 0.3], ...
        'Margin', 6);

    hold off;

    % Save if needed
    if do_save == true
        save_name = [title_str, '_bar.png'];
        save_path = fullfile(save_folder, save_name);
        exportgraphics(gcf, save_path, 'Resolution', 300);
        close;
    end
end
