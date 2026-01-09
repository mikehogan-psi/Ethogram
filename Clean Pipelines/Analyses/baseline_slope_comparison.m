% BASELINE_SLOPE_COMPARISON
% Calculates baseline firing (pre stim) and then  compares this to firing
% during and post stim, plotted on log-log axes
% Then fits a linear regression model to examine relationships between
% baseline firing and firing during these periods: 
% log10(Response)=β0 + β1log10(Baseline) + β2Group + β3(log10(Baseline)*Group)

master_directory = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)";

session= 'Extinction';

mice_to_analyse = [2 3 4 5 6 7 8];

num_mice = max(mice_to_analyse);

trial_type = 1; % 1 = loom, 0 = flash

received_psi = [3 5 7 8];
received_veh = [2 4 6];

%% Load processed neural and behavioural data

mouse_files = dir(fullfile(master_directory, 'Mouse*'));

mouse_paths = cell(numel(mouse_files), 1);

for mouse = mice_to_analyse
    current_mouse_folder = mouse_files(mouse).folder;
    current_mouse_name = mouse_files(mouse).name;
    current_mouse_path = fullfile(current_mouse_folder, current_mouse_name,...
        session, 'Combined Data');
    mouse_paths{mouse} = current_mouse_path;
end

combined_data_matrices = cell(numel(mouse_files), 1);

for mouse = mice_to_analyse
    current_file_to_load = dir(fullfile(mouse_paths{mouse}, '*GLM_data*'));
    current_data = load(fullfile(current_file_to_load.folder, current_file_to_load.name));
    combined_data_matrices{mouse} = current_data;
end
%%
num_bins = 66;

pre_bins = 1:19;
during_bins = 20:27;
post_bins = 28:66;

all_baseline = [];
all_during    = [];
all_post   = [];
all_group    = [];   % 0 = veh, 1 = psi
all_mouse    = []; 

for mouse = mice_to_analyse
    current_data = combined_data_matrices{mouse}.X_table_full;
    
    % Select correct trials
    current_data = current_data(current_data.TrialIdentifier == trial_type, :);
    cellIDs = unique(current_data.CellID);

    for idx = 1:length(cellIDs)
        current_neuron = cellIDs(idx);
        data_idx = current_data.CellID == current_neuron;
        spikes = current_data.SpikeCount(data_idx);
        bins = repmat(1:66, 1, length(spikes)/num_bins)';

        % Baseline / early / late means for this cell
        pre_idx = bins >= pre_bins(1) & bins <= pre_bins(end);
        baseline = mean(spikes(pre_idx));
        during_idx = bins >= during_bins(1) & bins <= during_bins(end);
        during = mean(spikes(during_idx));
        post_idx = bins >= post_bins(1) & bins <= post_bins(end);
        post = mean(spikes(post_idx));
        
        all_baseline(end + 1, 1) = baseline;
        all_during(end + 1, 1) = during;
        all_post(end + 1, 1) = post;

        if ismember(mouse, received_psi)
            all_group(end + 1, 1) = 1;
        elseif ismember(mouse, received_veh)
            all_group(end + 1, 1) = 0;
        end
        all_mouse(end + 1, 1) = mouse;
    end

end

%%
eps_val = 1e-3;

B = all_baseline;
E = all_during;
L = all_post;
G = all_group;   % 0 = veh, 1 = psi

veh_col = [0 0.4470 0.7410];
psi_col = [1 0 0];

%% ----- BASELINE vs DURING -----
mask_BE = ~isnan(B) & ~isnan(E);

idxVeh_BE = mask_BE & (G == 0);
idxPsi_BE = mask_BE & (G == 1);

x_veh = log10(B(idxVeh_BE) + eps_val);
y_veh = log10(E(idxVeh_BE) + eps_val);

x_psi = log10(B(idxPsi_BE) + eps_val);
y_psi = log10(E(idxPsi_BE) + eps_val);

figure('Color','w'); hold on;
hVeh = scatter(x_veh, y_veh, 10, veh_col, 'filled', 'MarkerFaceAlpha', 0.3);
hPsi = scatter(x_psi, y_psi, 10, psi_col, 'filled', 'MarkerFaceAlpha', 0.3);

xlabel('log_{10}(baseline FR)');
ylabel('log_{10}(during FR)');
title(['Baseline vs during loom per cell: ' session]);
box off; ax = gca; ax.TickDir = 'out';

% fit lines (log–log) separately for each group
p_veh = polyfit(x_veh, y_veh, 1);
p_psi = polyfit(x_psi, y_psi, 1);

x_line = linspace(min([x_veh; x_psi]), max([x_veh; x_psi]), 100);
y_fit_veh = polyval(p_veh, x_line);
y_fit_psi = polyval(p_psi, x_line);

hLineVeh = plot(x_line, y_fit_veh, 'Color', veh_col, 'LineWidth', 2);
hLinePsi = plot(x_line, y_fit_psi, 'Color', psi_col, 'LineWidth', 2);

% hide lines from legend, keep only scatters
set(hLineVeh, 'HandleVisibility', 'off');
set(hLinePsi, 'HandleVisibility', 'off');
legend([hVeh hPsi], {'Vehicle','Psilocybin'}, 'Location','best');
%%
% ----- STATS: BASELINE vs DURING -----
mask_BE = ~isnan(B) & ~isnan(E);

x_BE = log10(B(mask_BE) + eps_val);
y_BE = log10(E(mask_BE) + eps_val);
g_BE = G(mask_BE);                     % 0 = veh, 1 = psi

Group_BE = categorical(g_BE, [0 1], {'Veh','Psi'});

T_BE = table(x_BE, y_BE, Group_BE, ...
    'VariableNames', {'logB','logDuring','Group'});

% Linear model with interaction: different slopes/intercepts per group
mdl_BE = fitlm(T_BE, 'logDuring ~ logB * Group');

disp('===== Baseline vs DURING (loom) =====');
disp(mdl_BE);              % coefficients, p-values, R^2
disp('ANOVA for interaction term:');
disp(anova(mdl_BE, 'summary'));


%%
% ----- BASELINE vs POST -----
mask_BL = ~isnan(B) & ~isnan(L);

idxVeh_BL = mask_BL & (G == 0);
idxPsi_BL = mask_BL & (G == 1);

x2_veh = log10(B(idxVeh_BL) + eps_val);
y2_veh = log10(L(idxVeh_BL) + eps_val);

x2_psi = log10(B(idxPsi_BL) + eps_val);
y2_psi = log10(L(idxPsi_BL) + eps_val);

figure('Color','w'); hold on;
hVeh2 = scatter(x2_veh, y2_veh, 10, veh_col, 'filled', 'MarkerFaceAlpha', 0.3);
hPsi2 = scatter(x2_psi, y2_psi, 10, psi_col, 'filled', 'MarkerFaceAlpha', 0.3);

xlabel('log_{10}(baseline FR)');
ylabel('log_{10}(post FR)');
title(['Baseline vs post stim per cell: ' session]);
box off; ax = gca; ax.TickDir = 'out';

p2_veh = polyfit(x2_veh, y2_veh, 1);
p2_psi = polyfit(x2_psi, y2_psi, 1);

x2_line = linspace(min([x2_veh; x2_psi]), max([x2_veh; x2_psi]), 100);
y2_fit_veh = polyval(p2_veh, x2_line);
y2_fit_psi = polyval(p2_psi, x2_line);

hLineVeh2 = plot(x2_line, y2_fit_veh, 'Color', veh_col, 'LineWidth', 2);
hLinePsi2 = plot(x2_line, y2_fit_psi, 'Color', psi_col, 'LineWidth', 2);

set(hLineVeh2, 'HandleVisibility', 'off');
set(hLinePsi2, 'HandleVisibility', 'off');
legend([hVeh2 hPsi2], {'Vehicle','Psilocybin'}, 'Location','best');
%%
% ----- STATS: BASELINE vs POST -----
mask_BL = ~isnan(B) & ~isnan(L);

x_BL = log10(B(mask_BL) + eps_val);
y_BL = log10(L(mask_BL) + eps_val);
g_BL = G(mask_BL);

Group_BL = categorical(g_BL, [0 1], {'Veh','Psi'});

T_BL = table(x_BL, y_BL, Group_BL, ...
    'VariableNames', {'logB','logPost','Group'});

mdl_BL = fitlm(T_BL, 'logPost ~ logB * Group');

disp('===== Baseline vs POST (late) =====');
disp(mdl_BL);
disp('ANOVA for interaction term:');
disp(anova(mdl_BL, 'summary'));
