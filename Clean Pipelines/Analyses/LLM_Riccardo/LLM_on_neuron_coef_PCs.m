session = 'Renewal';

if strcmp(session, 'Extinction')
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\Model\extinction_coefficients.mat")
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ext\cluster_assignments_extinction.mat")
    B = B_ext;
    cluster_assignments_table = cluster_assignments_extinction;
elseif strcmp(session, 'Renewal')
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\Model\renewal_coefficients.mat")
    load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ren\cluster_assignments_renewal_to_extinction.mat") 
    B = B_ren;
    cluster_assignments_table = cluster_assignments_renewal;
end

%% Calculate gain and direction of coef vectors, filter out bad cells
% Calculate coefficient vector gain per cell
B_norm = vecnorm(B, 2, 2); 
% Normalise coefficient vectors by gain to get just coefficient directions
B_dir = B ./ B_norm; 

% Filter out bad cells (nans and zeros in norm and direction)
bad_mask = ~isfinite(B_norm) | B_norm == 0 | any(~isfinite(B_dir), 2);

B(bad_mask, :) = [];
B_dir(bad_mask, :) = [];
B_norm(bad_mask) = [];
cluster_assignments_table(bad_mask, :) = [];

% Convert factors to categorical for LLM input
cluster_assignments_table.MouseID = categorical(cluster_assignments_table.MouseID);
cluster_assignments_table.ClusterID = categorical(cluster_assignments_table.ClusterID);
cluster_assignments_table.Treatment = categorical(cluster_assignments_table.Treatment);

%% Run PCA on coefficient vector directions
[coeff, score, latent, ~, explained] = pca(B_dir, 'Centered', true);

% Choose PCs that explain 90% of variance (cap at 10)
cum_exp = cumsum(explained);

% Find number of PCs that make up 90% of variance
num_pcs = find(cum_exp >= 90, 1, 'first');

% Limit max to 10
num_pcs = min(num_pcs, 10);

pcs = score(:, 1:num_pcs);

%% Run one LMM per vector direction PC and perform ANOVA on LMM results

lme_per_pc = cell(num_pcs, 1);

p_interaction = nan(num_pcs, 1);
p_treatment = nan(num_pcs, 1);

for pc = 1:num_pcs

    current_T = cluster_assignments_table;
    current_T.Component = pcs(:, pc);
    
    lme_per_pc{pc} = fitlme(current_T, 'Component ~ Treatment*ClusterID + (1|MouseID)', ...
                      'FitMethod','REML');

    anova_results = anova(lme_per_pc{pc}, 'DFMethod','Satterthwaite');
    terms = string(anova_results.Term);
    row_int = find(terms == 'Treatment:ClusterID');
    row_treat = find(terms == 'Treatment');
    p_interaction(pc) = anova_results.pValue(row_int(1));
    p_treatment(pc) = anova_results.pValue(row_treat(1));

end

% Correct p values with false-disovery rate correction
p_interaction_fdr = mafdr(p_interaction, 'BHFDR', true);
p_treatment_fdr = mafdr(p_treatment, 'BHFDR', true);

% Store the adjusted p-values in the results table
results_table_int = table((1:num_pcs)', p_interaction, p_interaction_fdr, ...
                      'VariableNames', {'PC', 'P_Value_int', 'FDR_P_Value'});

results_table_treat = table((1:num_pcs)', p_treatment, p_treatment_fdr, ...
                      'VariableNames', {'PC', 'P_Value_treat', 'FDR_P_Value'});

% Display results
% --- Per-PC ANOVA p-values (Treatment×Cluster interaction) ---
disp('============================================================')
disp('LMM on B_dir PC scores: ANOVA p-values for Treatment×ClusterID per PC')
disp('Model: Component ~ Treatment*ClusterID + (1|MouseID)')
disp('Outcome: PC scores from PCA on B_dir (direction-normalised coefficients)')
disp('============================================================')
disp(results_table_int)

% --- Per-PC ANOVA p-values (main effect of Treatment) ---
disp('============================================================')
disp('LMM on B_dir PC scores: ANOVA p-values for main effect of Treatment per PC')
disp('Model: Component ~ Treatment*ClusterID + (1|MouseID)')
disp('Outcome: PC scores from PCA on B_dir (direction-normalised coefficients)')
disp('============================================================')
disp(results_table_treat)

% Save pc identities with significant treatment:cluster interaction
signif_pcs = find(results_table_int.FDR_P_Value <= 0.05); 


%% Run LMM on vector norms (gain) and display with title

norm_T = cluster_assignments_table;
norm_T.Norm = B_norm;

lme_norm = fitlme(norm_T, 'Norm ~ Treatment*ClusterID + (1|MouseID)', 'FitMethod', 'REML');
anova_norm_results = anova(lme_norm, 'DFMethod', 'Satterthwaite');

disp('============================================================')
disp('LMM on B_norm (gain): fixed effects table')
disp('Model: Norm ~ Treatment*ClusterID + (1|MouseID)')
disp('Outcome: B_norm = L2 norm of each neuron''s coefficient vector (gain/magnitude)')
disp('============================================================')
disp(lme_norm.Coefficients)

%% Run cluster-specific contrasts per PC (within-cluster treatment effects)
pc = 5;

current_T = cluster_assignments_table;
current_T.Component = pcs(:, pc);

current_T.ClusterID = removecats(current_T.ClusterID);
clusters = categories(current_T.ClusterID);
num_clusters = numel(clusters);


lme_per_signif_pc = fitlme(current_T, 'Component ~ Treatment*ClusterID + (1|MouseID)', 'FitMethod','REML');

coef_names = string(lme_per_signif_pc.CoefficientNames);
coefs  = fixedEffects(lme_per_signif_pc);
covariance = lme_per_signif_pc.CoefficientCovariance;

main_treatment_coef = coef_names(startsWith(coef_names,"Treatment_") & ~contains(coef_names,":"));
main_treatment_coef = main_treatment_coef(1);

treatment_effect_est = nan(num_clusters, 1);
treatment_effect_se = nan(num_clusters, 1);
treatment_effect_p = nan(num_clusters, 1);

for cluster_idx = 1:num_clusters
    contrast_vector = zeros(1, numel(coef_names));
    contrast_vector(coef_names == main_treatment_coef) = 1;

    if cluster_idx ~= 1
        interaction_name = main_treatment_coef + ":ClusterID_" + string(clusters{cluster_idx});
        contrast_vector(coef_names == interaction_name) = 1;
    end

    treatment_effect_est(cluster_idx) = contrast_vector * coefs;
    treatment_effect_se(cluster_idx) = sqrt(contrast_vector * covariance * contrast_vector');
    treatment_effect_p(cluster_idx) = coefTest(lme_per_signif_pc, contrast_vector);
end

treatment_effect_p_FDR = mafdr(treatment_effect_p, 'BHFDR', true);

results_by_cluster = table(clusters, treatment_effect_est, treatment_effect_se, ...
    treatment_effect_p, treatment_effect_p_FDR, ...
    'VariableNames', {'ClusterID','TreatEffect','SE','pValue','qValue'});

disp('============================================================')
disp(['Cluster-specific treatment contrasts on B_dir PC scores (within each ClusterID) | PC = ' num2str(pc)])
disp('Model fit: Component ~ Treatment*ClusterID + (1|MouseID)')
disp('Contrast tested per cluster: (psilocybin - vehicle) within that cluster')
disp('Outcome: PC score from PCA on B_dir (direction-normalised coefficients)')
disp('============================================================')
disp(results_by_cluster)


%% Inspect loading of PCs to interpret them
coef_names = {'Grooming','Rearing','Darting','Freezing','Velocity', ...
                 'TimeBin','Trial','TrialIdentifier','StimOn','StimOn_x_Type'};
pc = 5;
loadings = coeff(:,pc);

[~, idx] = sort(abs(loadings), 'descend');
topN = 10;

disp('============================================================')
disp(['PCA loadings for B_dir | PC = ' num2str(pc) ' | Top ' num2str(topN) ' predictors by |loading|'])
disp('Loadings indicate which original GLM predictors define this PC axis in coefficient-direction space')
disp('============================================================')
disp(table(coef_names(idx)', loadings(idx), 'VariableNames',...
    {'Predictor','Loading'}))
