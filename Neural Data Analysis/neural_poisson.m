%% Neural Logistic regression analyis 

% This script fits a generalized linear mixed-effects model (GLMM) in the from of a logistic function 
% to the firing patterns of individual neurons

% The goal is to identify which factors are most influencing the neurons firing pattersn
% the ivestigated factors are: 
%       - 1. timebin (i.e. time within each trial)
%       - 2. trials number 

%       - 3. behaviours: freezing, darting rearing
%       - 4. treatment groups: psilocybin or vehicle

%% Directory Setup

% Define folder path that contains neronal spiking data 
% (i.e. firing rate matrices and response groups obtained from neural_data_analysis_pipeline.m)

  neuronal_data_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\spiking_data\';

% Define folder path that contains kilosorted dats of this session
    kilosort_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\Extinction\kilosort4\'; % 


    triggers_path = 'Z:\Abi\neuronal_data\mouse_2\processed_data_extinction\concatinated_triggers\';

% Define filepath to save processed data to and name of saved data
    savefile_name_glmm = 'extinction_neuronal_cofficients_poisson_';
    savefile_name_data_tables = 'extinction_data_tables_';
    save_folder = 'Z:\Abi\neuronal_data\mouse_2\Mike processed data\';
%% Load data

% get important variables
  N_trials  = 20; % total number of trials
  fs = 30000;  % neuronal data sampling rate  
  tpre = 10;        % Time before each event to include in the window (i.e. pre-stimulus time)
  tpost = 24;       % Time after each event to include.               (i.e. post-stimulus time)

% neuronal data - converting kilosort python output files into matlab format
    clu = readNPY([kilosort_dir '\spike_clusters.npy']); % loads cluster number of each detected spike
    spk = readNPY([kilosort_dir '\spike_times.npy']);    % loads sample number of each detected spike 

    spk = double(spk)/fs;   % converts sample number (by dividing by sampling rate) into seconds
    clu_val = unique(clu);  % removes doubles (so gets list of all individual clusters i.e. potential neuronal cells - starts with 0)
    N_neurons = numel(clu_val); % get number of detected clusters i.e neurons

% load list of Loom-response-types 
  load([neuronal_data_path 'mouse2_extinction_LOOMresp_neurons'], 'resp_duringLOOM_nr', 'resp_postLOOM_nr')

% load extracted loom and flash (i.e. stimuli onset) events
  load([triggers_path 'mouse2_extinction_extracted_events'],'evt_loom')

% Load freezing data
  load("D:\PhD 2nd Year\Cohort 4 Mouse 2 DLC Data Temp\cam5\mouse2_loom_freezing.mat")  
  freeze_matrix = mouse2_freeze_loom;

% Create variable with all responsive neurons
resp_LOOM_nr = [resp_duringLOOM_nr; resp_postLOOM_nr];
resp_LOOM_nr = unique(resp_LOOM_nr);


% Specify which type of responsive neurons you want to look at and save
% name
    current_group = resp_LOOM_nr;
    append_to_save = 'resp_LOOM_nr.mat';

% Generates filepaths for saving   
    savefile_name_glmm = [savefile_name_glmm,  append_to_save];
    save_path_glmm = fullfile(save_folder, savefile_name_glmm);
    
    savefile_name_data_tables = [savefile_name_data_tables, append_to_save];
    save_path_data_tables = fullfile(save_folder, savefile_name_data_tables);
  %% Fit logistic regression of during and post-stim periods for each neuron

  % define time variables to fit glmm for duing stim and post stim periods
           period = ["during", "post"]; % -> first emement during stim, secodn post-stim 
    period_length = [50,300];           % stim period is 50 frames, post-stim is period 300 frames
          t_start = [0, 3.3];           % starts at 0 sec,  starts at 3.3 sec            
          t_end = [3.3, 23.3];          % ends at 3.3 sec, ends at 23.3 sec
          num_bins = [20, 20];          % adapt nr of bins a wanted
  
  % initialise output matrix to hold glmm results and input tables
    glmm_output = zeros(size(current_group,1), 6); % 6 columns -> hold beta values for TimeBinStd, TrialStd and Freezing for during (1-3) and post (4-6) glmm
    all_data_tables = cell(size(current_group,1), 2);

for n = 1:length(current_group) % loop through each neuron
    idx = current_group(n);
    disp(['computing GLMM of neuron: ' num2str(idx-1)]);

   for m = 1:2 % run GLMM for 1 -> during stim period, 2 -> post-stim period  

     disp([period(m) ' stimulus period' ]);
       
     % get binsize depending on how many bins you want
       num_bins_current = num_bins(m);                  % define number of bins
       bin_size = period_length(m) / num_bins_current;  % get binsize dependent on length of investigated period (if longer, greater bin size required)
       bin_size_s = bin_size/15;                  % gets binsize in seconds (15fps)
       % disp(['Binsize is ' num2str(bin_size) ' seconds'])  % display binsize

     % extracts all spikes from current cluster/neuron
       tsp = spk(clu==clu_val(idx)); 
      
     % extract the spikes that are occuring in each timebin fo each trial
       [~, ~, t_loom, fr, spikes]   = raster_NM(tsp,evt_loom,tpre,tpost,bin_size_s,false);

     % extract firing rates in time_bins that are occuring within neuronal response period
       idx_during = find(t_loom >= t_start(m) & t_loom <= t_end(m));
       fr_period = fr(:, idx_during, :);

       spike_count = reshape(fr_period', [], 1); 
       
     % Extract freezing within bin bounds
        frames_start = round(t_start(m) * 15) + 1;
        frames_end   = round(t_end(m) * 15);
        cropped_freezing = freeze_matrix(:, frames_start:frames_end);
     % Initialise binned freezing matrix
        binned_freezing = zeros(N_trials, num_bins_current);
     
     % Calculate freezing over bins
        for trial = 1:size(freeze_matrix, 1)
            for bin = 1:num_bins_current
                start_idx = round((bin-1)*bin_size+1);
                end_idx = round(bin*bin_size);
                binned_freezing(trial, bin) = mean(cropped_freezing(trial, start_idx:end_idx));
            end
        end

     % Flatten into column vector
     binned_freezing = reshape(binned_freezing, [], 1);

    % get time bin and trial indices
     time_bins = repmat((1:num_bins_current)', N_trials, 1);   % Time bin indices
     trials_binned = repelem((1:N_trials)', num_bins_current); % Trial numbers
     neuron_number = repelem(idx-1, N_trials * num_bins_current)'; % Neuron numbers

     data_binned = table(spike_count, time_bins, trials_binned, binned_freezing, ...
                'VariableNames', {'Firing', 'TimeBin', 'Trial', 'Freezing'});
     data_binned.TimeBinStd = zscore(data_binned.TimeBin);
     data_binned.TrialStd = zscore(data_binned.Trial);
     data_binned.NeuronNumber = neuron_number;
   
     glme = fitglm(data_binned, 'Firing ~ TimeBinStd + TrialStd + Freezing', ...
              'Distribution', 'poisson');
    
    % Display model summary
    % disp(['GLMM results of neuron' num2str(n) 'for' period(m) 'stimulus period'])
    % disp(glme);
    
 % Fit null model (intercept only)
    null_model = fitglm(data_binned, 'Firing ~ 1', ...
                    'Distribution', 'poisson');


    % load results for this neuron into output vector
    if m == 1
       glmm_output(n, 1:3) = glme.Coefficients.Estimate(2:4)'; % load beta values for during-stim glmm
    else
       glmm_output(n, 4:6) = glme.Coefficients.Estimate(2:4)'; % load beta values for post-stim glmm
    end
    
   % Store neuron number for clustering later (-1 to match kilosort IDs)
    glmm_output(n, 7) = idx-1;

    all_data_tables{n, m} = data_binned;

   % -------- Goodness of fit metrics --------

    % Store deviance and pseudo-R²
    model_deviance = glme.Deviance;
    null_deviance = null_model.Deviance;
    pseudo_R2 = 1 - (model_deviance / null_deviance);

    % Store deviance and pseudo-R²
    model_deviances(n, m) = model_deviance;
    pseudo_R2_vals(n, m) = pseudo_R2;

   end
   

end    



%% Find neurons with bad model fit

bad_fit_idx = pseudo_R2_vals(:,1) < 0.01 & pseudo_R2_vals(:,2) < 0.01;
good_fit_neurons = glmm_output(:, 7);
good_fit_neurons = good_fit_neurons(~bad_fit_idx);

% Exclude these from clustering
all_data_tables = all_data_tables(~bad_fit_idx, :);
glmm_output = glmm_output(~bad_fit_idx, :);

save(save_path_glmm, 'glmm_output')
save(save_path_data_tables, "all_data_tables")
%% Plot: Actual vs Predicted Freezing Across Trials

trial_val = unique(data_binned.Trial);
group_val = unique(data_binned.Group);
Nval = numel(trial_val);

% Recalculate predicted probabilities
yhat = fitted(glme);

% Preallocate arrays
ym_actual = zeros(Nval, numel(group_val));
ym_pred = zeros(Nval, numel(group_val));
sem_actual = zeros(Nval, numel(group_val));

for g = 1:numel(group_val)
    for idx = 1:Nval
        % Logical indices
        group_idx = data_binned.Group == group_val(g);
        trial_idx = data_binned.Trial == trial_val(idx);
        idx = group_idx & trial_idx;

        % Get mouse IDs in this group
        mice = unique(data_binned.Mouse(idx));
        n_mice = numel(mice);
        per_mouse_mean = zeros(n_mice, 1);

        for m = 1:n_mice
            mouse_idx = data_binned.Mouse == mice(m);
            mouse_trial_idx = idx & mouse_idx;
            per_mouse_mean(m) = mean(data_binned.Freezing(mouse_trial_idx));
        end

        % Mean & SEM across mice
        ym_actual(idx, g) = mean(per_mouse_mean);
        sem_actual(idx, g) = std(per_mouse_mean) / sqrt(n_mice);

        % Predicted freezing (mean across observations)
        ym_pred(idx, g) = mean(yhat(idx));
    end
end

% COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM ERROR BARS
figure; hold on;

errorbar(trial_val, ym_actual(:,1), sem_actual(:,1), 'ko-', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1);
plot(trial_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);

errorbar(trial_val, ym_actual(:,2), sem_actual(:,2), 'ro-', ...
    'MarkerFaceColor', 'r', 'LineWidth', 1);
plot(trial_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);

% (Rest of your plotting code stays the same)

xlabel('Trial Number', 'FontWeight', 'bold');
ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
        'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
       'Location', 'Best');
ylim([0 1]);
xlim([1 20]);
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
ax.LineWidth = 2;
ax.GridLineWidth = 0.5;
grid on;
hold off;

% print(gcf, 'freezing_across_trials_Acq_ren_GLMEM.png', '-dpng', '-r300');

% % COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM SHADED AREAS
% figure; hold on;
% 
% % Define colours
% vehicle_colour = [0 0 0];        % black
% psilo_colour   = [0.8 0 0];      % dark red
% 
% % Plot shaded SEM area: Vehicle
% x = trial_val;
% y1 = ym_actual(:,1) - sem_actual(:,1);
% y2 = ym_actual(:,1) + sem_actual(:,1);
% fill([x; flipud(x)], [y1; flipud(y2)], vehicle_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Plot shaded SEM area: Psilocybin
% y1 = ym_actual(:,2) - sem_actual(:,2);
% y2 = ym_actual(:,2) + sem_actual(:,2);
% fill([x; flipud(x)], [y1; flipud(y2)], psilo_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Now plot the means and predicted lines
% plot(trial_val, ym_actual(:,1), 'k-', 'LineWidth', 1.5);
% plot(trial_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);
% 
% plot(trial_val, ym_actual(:,2), 'r-', 'LineWidth', 1.5);
% plot(trial_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);
% 
% xlabel('Trial Number', 'FontWeight', 'bold');
% ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
% legend({'Vehicle SEM', 'Psilocybin SEM', ...
%         'Vehicle (Actual)', 'Vehicle (Predicted)', ...
%         'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
%        'Location', 'Best');
% ylim([0 1]);
% xlim([1 20]);
% ax = gca;
% ax.FontWeight = 'bold';
% ax.FontSize = 12;
% ax.LineWidth = 2;
% ax.GridLineWidth = 0.5;
% grid on;
% hold off;

% % print(gcf, 'freezing_across_trials_Acq_ren_GLMEM.png', '-dpng', '-r300');

%% Plot: Actual vs Predicted Freezing Across Time Bins
bin_val = unique(data_binned.TimeBin);
Nval = numel(bin_val);

ym_actual = zeros(Nval, numel(group_val));
ym_pred = zeros(Nval, numel(group_val));
sem_actual = zeros(Nval, numel(group_val));

for g = 1:numel(group_val)
    for idx = 1:Nval
        group_idx = data_binned.Group == group_val(g);
        time_idx = data_binned.TimeBin == bin_val(idx);
        idx = group_idx & time_idx;

        mice = unique(data_binned.Mouse(idx));
        n_mice = numel(mice);
        per_mouse_mean = zeros(n_mice, 1);

        for m = 1:n_mice
            mouse_idx = data_binned.Mouse == mice(m);
            mouse_time_idx = idx & mouse_idx;
            per_mouse_mean(m) = mean(data_binned.Freezing(mouse_time_idx));
        end

        ym_actual(idx, g) = mean(per_mouse_mean);
        sem_actual(idx, g) = std(per_mouse_mean) / sqrt(n_mice);
        ym_pred(idx, g) = mean(yhat(idx));
    end
end

% COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM ERROR BARS
figure; hold on;

errorbar(bin_val, ym_actual(:,1), sem_actual(:,1), 'ko-', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1);
plot(bin_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);

errorbar(bin_val, ym_actual(:,2), sem_actual(:,2), 'ro-', ...
    'MarkerFaceColor', 'r', 'LineWidth', 1);
plot(bin_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);




xlabel('Time Bin Number (2s)', 'FontWeight', 'bold');
ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
legend({'Vehicle (Actual)', 'Vehicle (Predicted)', ...
        'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
       'Location', 'Best');
ylim([0 1]);
xlim([1 9]);
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
ax.LineWidth = 2;
ax.GridLineWidth = 0.5;
grid on;
hold off;

% print(gcf, 'freezing_across_timebins_acq_ren_GLMEM.png', '-dpng', '-r300');

% % COMMENT/UNCOMMENT TO CHANGE ERROR BARS TO TOGGLE SEM SHADED AREAS
% figure; hold on;
% 
% % Define colours
% vehicle_colour = [0 0 0];        % black
% psilo_colour   = [0.8 0 0];      % dark red
% 
% % Plot shaded SEM area: Vehicle
% x = bin_val;
% y1 = ym_actual(:,1) - sem_actual(:,1);
% y2 = ym_actual(:,1) + sem_actual(:,1);
% fill([x; flipud(x)], [y1; flipud(y2)], vehicle_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Plot shaded SEM area: Psilocybin
% y1 = ym_actual(:,2) - sem_actual(:,2);
% y2 = ym_actual(:,2) + sem_actual(:,2);
% fill([x; flipud(x)], [y1; flipud(y2)], psilo_colour, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Now plot the means and predicted lines
% plot(bin_val, ym_actual(:,1), 'k-', 'LineWidth', 1.5);
% plot(bin_val, ym_pred(:,1), 'k--', 'LineWidth', 1.5);
% 
% plot(bin_val, ym_actual(:,2), 'r-', 'LineWidth', 1.5);
% plot(bin_val, ym_pred(:,2), 'r--', 'LineWidth', 1.5);
% 
% xlabel('Bin Number', 'FontWeight', 'bold');
% ylabel('Mean Freezing Probability', 'FontWeight', 'bold');
% legend({'Vehicle SEM', 'Psilocybin SEM', ...
%         'Vehicle (Actual)', 'Vehicle (Predicted)', ...
%         'Psilocybin (Actual)', 'Psilocybin (Predicted)'}, ...
%        'Location', 'Best');
% ylim([0 1]);
% xlim([1 8]);
% ax = gca;
% ax.FontWeight = 'bold';
% ax.FontSize = 12;
% ax.LineWidth = 2;
% ax.GridLineWidth = 0.5;
% grid on;
% hold off;

% % print(gcf, 'freezing_across_timebins_acq_ren_GLMEM.png', '-dpng', '-r300');