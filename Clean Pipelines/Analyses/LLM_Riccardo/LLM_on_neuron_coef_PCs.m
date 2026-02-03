%% Load coefficient matrix and cluster assignments table

load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\Model\extinction_coefficients.mat")
load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Clustering\GLM_output_ext\cluster_assignments_extinction.mat")

%% Assign universal variable names based on session
session = 'Extinction';

if strcmp(session, 'Extinction')
    B = B_ext;
    cluster_assignments_table = cluster_assignments_extinction;
elseif strcmp(session, 'Renewal')
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

%% Run PCA on coefficient vector directions
[coeff, score, latent, ~, explained] = pca(B_dir, 'Centered', true);

% Choose PCs that explain 90% of variance (cap at 10)
cum_exp = cumsum(explained);

% Find number of PCs that make up 90% of variance
num_pcs = find(cum_exp >= 90, 1, 'first');

% Limit max to 10
num_pcs = min(num_pcs, 10);

pcs = score(:, 1:num_pcs);

%% Run one LMM per PC and perform ANOVA on LMM results

lme_per_pc = cell(num_pcs, 1);

p_values = nan(num_pcs, 1);

for pc = 1:num_pcs

    current_T = cluster_assignments_table;
    current_T.Component = pcs(:, pc);
    
    lme_per_pc{pc} = fitlme(current_T, 'Component ~ Treatment*ClusterID + (1|MouseID)', ...
                      'FitMethod','REML');

    anova_results = anova(lme_per_pc{pc}, 'DFMethod','Satterthwaite');
    row = find(contains(string(anova_results.Term),'Treatment:ClusterID'));
    p_values(pc) = anova_results.pValue(row(1));

end

% Correct p values with false-disovery rate correction
p_values_fdr = mafdr(p_values, 'BHFDR', true);

% Store the adjusted p-values in the results table
results_table = table((1:num_pcs)', p_values, p_values_fdr, ...
                      'VariableNames', {'PC', 'P_Value', 'FDR_P_Value'});

disp(results_table)