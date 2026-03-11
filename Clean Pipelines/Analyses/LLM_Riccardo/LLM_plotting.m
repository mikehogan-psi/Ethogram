figure('Color','w'); 
tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

% Scree
nexttile; 
plot(1:numel(explained), explained, '-o', 'LineWidth', 1.5);
xlabel('PC'); ylabel('% variance explained');
title('PCA on direction vectors: Scree');
xlim([1 10])
grid on;

% Cumulative
nexttile;
plot(1:numel(explained), cumsum(explained), '-o', 'LineWidth', 1.5); hold on;
yline(90,'--','90%','LineWidth',1.2);
xline(num_pcs,'--',sprintf('Chosen = %d', num_pcs),'LineWidth',1.2);
xlabel('PC'); ylabel('Cumulative % explained');
title('PCA on direction vectors: Cumulative');
xlim([1 10])
ylim([0 100]);
grid on;

%%
figure('Color','w'); 
tiledlayout(2,1, 'Padding','compact', 'TileSpacing','compact');

% Interaction
nexttile;
pcx = results_table_int.PC;
plot(pcx, results_table_int.P_Value_int, '-o', 'LineWidth', 1.5); hold on;
plot(pcx, results_table_int.FDR_P_Value, '-s', 'LineWidth', 1.5);
yline(0.05,'--','0.05');
xlabel('PC'); ylabel('p-value');
title('LMM ANOVA: Treatment×ClusterID per PC');
legend({'raw p','FDR q'}, 'Location','best');
set(gca,'YScale','log'); grid on;
ylim([0.004 1])

% Treatment main effect
nexttile;
pcx = results_table_treat.PC;
plot(pcx, results_table_treat.P_Value_treat, '-o', 'LineWidth', 1.5); hold on;
plot(pcx, results_table_treat.FDR_P_Value, '-s', 'LineWidth', 1.5);
yline(0.05,'--','0.05');
xlabel('PC'); ylabel('p-value');
title('LMM ANOVA: main Treatment effect per PC');
legend({'raw p','FDR q'}, 'Location','best');
set(gca,'YScale','log'); grid on;

sig = results_table_int.PC(results_table_int.FDR_P_Value <= 0.05);
disp("Significant PCs (interaction, FDR<=0.05): " + strjoin(string(sig), ", "));

%%
pc = 4;
T = cluster_assignments_table;
T.Component = pcs(:,pc);

% Collapse to one value per Mouse×Cluster (mean across cells)
G = findgroups(T.MouseID, T.ClusterID);
mouseMean = splitapply(@mean, T.Component, G);
mouseTreat = splitapply(@(x) x(1), T.Treatment, G);
mouseClu   = splitapply(@(x) x(1), T.ClusterID, G);

Tm = table(mouseMean, mouseTreat, mouseClu, ...
    'VariableNames', {'Component','Treatment','ClusterID'});

Tm.ClusterID = removecats(Tm.ClusterID);
Tm.Treatment = removecats(Tm.Treatment);

clusters = categories(Tm.ClusterID);
treats   = categories(Tm.Treatment);

nClu = numel(clusters);
nCols = ceil(sqrt(nClu));
nRows = ceil(nClu/nCols);

figure('Color','w');
for c = 1:nClu
    subplot(nRows,nCols,c); hold on;
    Tc = Tm(Tm.ClusterID==clusters{c},:);

    % Jittered scatter by treatment
    for t = 1:numel(treats)
        yt = Tc.Component(Tc.Treatment==treats{t});
        if isempty(yt), continue; end

        x0 = t;
        xj = x0 + 0.08*randn(size(yt));   % jitter
        plot(xj, yt, 'o', 'LineWidth', 1.2);

        % mean ± SEM
        mu = mean(yt,'omitnan');
        sem = std(yt,'omitnan') / sqrt(sum(isfinite(yt)));
        errorbar(x0, mu, sem, 'k', 'LineWidth', 1.5, 'CapSize', 10);
    end

    xlim([0.5 numel(treats)+0.5]);
    xticks(1:numel(treats));
    xticklabels(string(treats));
    ylabel(sprintf('PC%d score (mouse mean)', pc));
    title("Cluster " + string(clusters{c}));
    grid on;
end
sgtitle(sprintf('Between-subject: PC%d by Treatment within Cluster (mouse means)', pc));

%%
pc = 4;
R = results_by_cluster;

est = R.TreatEffect;
se  = R.SE;
ci95 = 1.96 * se;

y = 1:height(R);

figure('Color','w'); hold on;
errorbar(est, y, ci95, 'horizontal', 'o', 'LineWidth', 1.8, 'CapSize', 10);
xline(0,'--','LineWidth',1.2);

yticks(y);
yticklabels("Cluster " + string(R.ClusterID));
ylim([0.5 height(R)+0.5]);   % tight y

% Tight x-lims around the CI range
xmin = min(est - ci95); 
xmax = max(est + ci95);
pad = 0.15 * (xmax - xmin + eps);
xlim([xmin - pad, xmax + pad]);

xlabel(sprintf('Treatment effect (psilocybin - vehicle) on PC%d score', pc));
title(sprintf('Cluster-specific treatment contrasts (PC%d)', pc));
grid on;

% Mark FDR-significant clusters
sig = R.qValue <= 0.05;
plot(est(sig), y(sig), 's', 'MarkerSize', 9, 'LineWidth', 1.8);

legend({'Estimate \pm 95% CI','0','FDR sig'}, 'Location','best');

%%
Tn = cluster_assignments_table;
Tn.Norm = B_norm;

% Reduce to one value per Mouse×Treatment×Cluster (mean across cells)
G = findgroups(Tn.MouseID, Tn.Treatment, Tn.ClusterID);
meanNorm = splitapply(@mean, Tn.Norm, G);
mouseG   = splitapply(@(x) x(1), Tn.MouseID, G);
treatG   = splitapply(@(x) x(1), Tn.Treatment, G);
cluG     = splitapply(@(x) x(1), Tn.ClusterID, G);

Tmouse = table(mouseG, treatG, cluG, meanNorm, ...
    'VariableNames', {'MouseID','Treatment','ClusterID','MeanNorm'});

figure('Color','w');
boxchart(categorical(Tmouse.ClusterID), Tmouse.MeanNorm, 'GroupByColor', categorical(Tmouse.Treatment));
xlabel('Cluster'); ylabel('Mean B\_norm (per mouse)');
title('Gain (B\_norm): per-mouse means by Cluster and Treatment');
grid on;
legend(categories(categorical(Tmouse.Treatment)), 'Location','best');

%%
% --- Bar chart of PCA loadings for a chosen PC ---
pc = 4;                          % choose PC
loadings = coeff(:, pc);         % coeff from pca()
names = string(coef_names(:));   % predictor names (same order as columns in B_dir)

% Sort by absolute loading (so biggest contributors are at top)
[~, idx] = sort(abs(loadings), 'descend');
load_s = loadings(idx);
names_s = names(idx);

% Optional: show only top N
topN = min(10, numel(load_s));
load_s = load_s(1:topN);
names_s = names_s(1:topN);

figure('Color','w'); hold on;
barh(load_s, 'FaceColor', 'flat');   % horizontal bars
yline(0);                            % (optional) not super meaningful for barh; can remove
xline(0,'--','LineWidth',1.2);

yticks(1:topN);
yticklabels(names_s);
set(gca, 'YDir','reverse');          % biggest at top
xlabel(sprintf('Loading on PC%d', pc));
title(sprintf('Top %d PCA loadings (by |loading|) for PC%d', topN, pc));
grid on;
