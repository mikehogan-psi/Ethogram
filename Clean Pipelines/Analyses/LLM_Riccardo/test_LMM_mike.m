function[] = test_LMM_mike()

%% Simulation parameters
Ncell = 2000;      % number of cells
Nmice = 8;        % number of animals
Ncluster = 10;      % number of clusters
clu_bias = 0.02*randn(1,Ncluster);  
clu_treatment_bias = 0.1*randn(1,Ncluster);        

%% Assign animals
ind_mice = sort(unidrnd(Nmice, Ncell, 1), 'ascend');
U = zeros(Ncell, Nmice);
for n = 1:Ncell
    U(n, ind_mice(n)) = 1;
end

%% Assign clusters
ind_clu = unidrnd(Ncluster, Ncell, 1);
C = zeros(Ncell, Ncluster);
for n = 1:Ncell
    C(n, ind_clu(n)) = 1;
end

%% Random animal intercepts
animal = 0.1*randn(1, Nmice);

%% Treatment (assigned at animal level)
treatment = double(ind_mice > Nmice/2); % 0 = Control, 1 = Treated

%% Residual noise
noise = 0.05*randn(Ncell,1);

%% Generate outcome
y = C*clu_bias' + (C*clu_treatment_bias').*treatment + U*animal' + noise;

%% Create table for fitlme
T = table(y(:), categorical(treatment, [0 1], {'Control','Treated'}), ...
          categorical(ind_clu), categorical(ind_mice), ...
          'VariableNames', {'y','treatment','cluster','animal'});

%% Fit linear mixed-effects model
lme = fitlme(T, 'y ~ treatment*cluster + (1|animal)', 'FitMethod','REML');

%% Extract estimated fixed effect for treatment
beta_hat = fixedEffects(lme);
ref_idx = find(strcmp(lme.CoefficientNames, 'treatment_Treated'));
ref_effect = beta_hat(ref_idx);  % effect in reference cluster

interaction_idx = contains(lme.CoefficientNames, 'treatment_Treated:cluster_');
interaction_effects = beta_hat(interaction_idx);

% Absolute treatment effect per cluster:
beta_hat_absolute = zeros(Ncluster,1);
beta_hat_absolute(1) = ref_effect;                 % reference cluster
beta_hat_absolute(2:end) = ref_effect + interaction_effects;


%% Display results
figure; hold on;
plot(clu_treatment_bias,'.:');
plot(beta_hat_absolute,'.:');
line([0 Ncluster+1],[0 0],'Color','k')
legend('original','estimated');
xlabel('#cluster')
ylabel('fixed effects: treatment*cluster')


%% Show all stats
Tcoef = lme.Coefficients;
disp(Tcoef)

%% Display p-values
figure; hold on;
pval = Tcoef.pValue([ref_idx find(interaction_idx)])';
pval(pval<0.00001) = 0.000001;
plot(log10(pval),'.:');
line([0 Ncluster+1],-2*[1 1],'Color','k')
xlabel('#cluster')
ylabel('log10(p-values): treatment*cluster');