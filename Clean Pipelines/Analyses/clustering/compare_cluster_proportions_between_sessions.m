% Load your assignment tables
load("D:\PhD 3rd Year\poisson_GLM_data_21_01_26\cluster_assignments_extinction.mat")
load("D:\PhD 3rd Year\renewal_glm_data\cluster_assignments_retention_to_extinction.mat")

cluster_assignments.Session = repmat("extinction", height(cluster_assignments), 1);
cluster_assignments_retention.Session = repmat("retention", height(cluster_assignments_retention), 1);

T = vertcat(cluster_assignments, cluster_assignments_retention);

T.Treatment = categorical(T.Treatment, [0 1], {'Veh','Psi'});
T.Session   = categorical(T.Session, {'extinction','retention'});
T.ClusterID = categorical(T.ClusterID);

% Handy group label for plotting
T.Group = categorical(strcat(string(T.Session), " - ", string(T.Treatment)));

%%
% counts per Group x Cluster
GC = groupcounts(T, {'Group','ClusterID'});
GC.Properties.VariableNames{'GroupCount'} = 'N';

% total per Group
GT = groupcounts(T, {'Group'});
GT.Properties.VariableNames{'GroupCount'} = 'N_total';

% join totals + compute proportions
GC = outerjoin(GC, GT, 'Keys', {'Group'}, 'MergeKeys', true);
GC.Prop = GC.N ./ GC.N_total;

% ---- CRUCIAL LINE: drop extra columns before unstack ----
GC2 = GC(:, {'Group','ClusterID','Prop'});

% wide: one row per Group, one column per cluster
P = unstack(GC2, 'Prop', 'ClusterID');

% fill missing props (clusters absent in some groups)
is_num = varfun(@isnumeric, P, 'OutputFormat','uniform');
P(:, is_num) = fillmissing(P(:, is_num), 'constant', 0);

xlab = P.Group;                   % now should be 4 unique groups
Y = table2array(P(:, is_num));

figure;
bar(xlab, Y, 'stacked');
xlabel('Session - Treatment');
ylabel('Proportion of cells');
title('Cluster composition (proportions)');

leg = string(P.Properties.VariableNames(is_num));
leg = "C" + erase(leg,"x");
legend(leg, 'Location','eastoutside');
%%
% Count cells per MouseID×Session×Treatment×ClusterID
G = groupcounts(T, {'MouseID','Session','Treatment','ClusterID'});
G.Properties.VariableNames{'GroupCount'} = 'N';

% Total cells per MouseID×Session×Treatment
Gtot = groupcounts(T, {'MouseID','Session','Treatment'});
Gtot.Properties.VariableNames{'GroupCount'} = 'N_total';

% Join totals back in
G = outerjoin(G, Gtot, 'Keys', {'MouseID','Session','Treatment'}, 'MergeKeys', true);

% Proportion per cluster within each mouse/session/treatment
G.Prop = G.N ./ G.N_total;

% ---- Plot: per cluster, dots per mouse + mean±SEM ----
clusters = categories(T.ClusterID);
x_groups = categories(categorical(strcat(string(G.Session), " - ", string(G.Treatment))));

figure;
for k = 1:numel(clusters)
    subplot(1, numel(clusters), k); hold on;

    this = G(G.ClusterID == clusters{k}, :);
    this.Group = categorical(strcat(string(this.Session), " - ", string(this.Treatment)));

    % jittered scatter
    for gi = 1:numel(x_groups)
        y = this.Prop(this.Group == x_groups{gi});
        x = gi + 0.08*randn(size(y));
        scatter(x, y, 25, 'filled');
        
        m = mean(y, 'omitnan');
        s = std(y, 'omitnan') / sqrt(sum(isfinite(y)));
        errorbar(gi, m, s, 'k', 'LineWidth', 1.5);
    end

    xlim([0.5 numel(x_groups)+0.5]);
    ylim([0 1]);
    set(gca, 'XTick', 1:numel(x_groups), 'XTickLabel', x_groups);
    xtickangle(30);
    ylabel('Proportion of cells');
    title("Cluster " + clusters{k});
end
sgtitle('Per-mouse cluster proportions (dots = mice)');
%%
clusters = categories(T.ClusterID);

figure;
for k = 1:numel(clusters)
    subplot(1, numel(clusters), k); hold on;

    this = G(G.ClusterID == clusters{k}, :);

    % Make wide table: one row per MouseID×Treatment, columns = extinction/retention
    W = unstack(this(:, {'MouseID','Treatment','Session','Prop'}), 'Prop', 'Session');

    % Some mice might be missing a session -> drop those
    if ~all(ismember({'extinction','retention'}, W.Properties.VariableNames))
        warning('Missing session columns for cluster %s', clusters{k});
        continue
    end

    ok = isfinite(W.extinction) & isfinite(W.retention);
    W = W(ok, :);

    d = W.retention - W.extinction;  % Δ per mouse

    % plot per mouse, split by treatment
    treatments = categories(W.Treatment);
    for ti = 1:numel(treatments)
        idx = W.Treatment == treatments{ti};
        x = ti + 0.08*randn(sum(idx),1);
        scatter(x, d(idx), 25, 'filled');
        m = mean(d(idx), 'omitnan');
        s = std(d(idx), 'omitnan')/sqrt(sum(idx));
        errorbar(ti, m, s, 'k', 'LineWidth', 1.5);
    end

    yline(0,'--k');
    xlim([0.5 2.5]);
    set(gca,'XTick',1:2,'XTickLabel',treatments);
    ylabel('\Delta proportion (ret - ext)');
    title("Cluster " + clusters{k});
end
sgtitle('Per-mouse change: retention − extinction');

%%
