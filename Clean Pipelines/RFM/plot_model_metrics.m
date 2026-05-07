load("W:\Mike\Neuropixels_Fear_Conditioning\Models\RFMs\Darting\Models\Darting_model_3.mat")
do_save = true;
behaviour = 'Darting'; % set manually (or parse from filename)

% ---- helper ----
mean_sem = @(x) deal(mean(x,'omitnan'), std(x,'omitnan')/sqrt(sum(isfinite(x))));

% ---- fold summaries ----
[mPRA, sePRA] = mean_sem(prauc_scores);
[mROC, seROC] = mean_sem(rocauc_scores);
[mBSS, seBSS] = mean_sem(brier_skill_scores);
[mPrev,sePrev]= mean_sem(prev_scores);

% ---- pooled PR curve + baseline ----
p = p1_all_test(:);
y = y_all_test(:);
p = min(max(p, 1e-6), 1-1e-6);

[rec, prec, ~, prAUC_pooled] = perfcurve(y, p, 1, 'xCrit','reca','yCrit','prec');
prev_pooled = mean(y==1);

% ---- pooled calibration ----
nBins = 10;
edges = linspace(0,1,nBins+1);
binMeanP = nan(nBins,1);
binFracY = nan(nBins,1);
binN     = zeros(nBins,1);

for bb = 1:nBins
    idx = (p >= edges(bb) & p < edges(bb+1));
    if bb == nBins
        idx = (p >= edges(bb) & p <= edges(bb+1));
    end
    binN(bb) = sum(idx);
    if binN(bb) > 0
        binMeanP(bb) = mean(p(idx));
        binFracY(bb) = mean(y(idx));
    end
end
ece = nansum((binN/sum(binN)) .* abs(binFracY - binMeanP));

% ---- styling ----
set(groot, ...
    'defaultAxesLineWidth', 1.5, ...
    'defaultAxesFontWeight','bold', ...
    'defaultAxesFontSize', 14, ...
    'defaultTextFontWeight','bold', ...
    'defaultTextFontSize', 14, ...
    'defaultLineLineWidth', 2);

% proPalette = [
%     0.090 0.380 0.670  % navy-blue
%     0.000 0.580 0.620  % teal
%     0.400 0.760 0.650  % seafoam
% ];

proPalette = [
    0.20 0.20 0.20  % dark grey
    0.45 0.45 0.45  % mid grey
    0.70 0.70 0.70  % light grey
];



%% ---------- FIGURE 1: Summary bars ----------
fig1 = figure('Color','w','Position',[100 100 450 380]); hold on;

vals  = [mPRA, mROC, mBSS];
errs  = [sePRA, seROC, seBSS];
names = {'PR-AUC','ROC-AUC','Brier skill'};

bh = bar(vals, 'FaceColor','flat');
bh.CData = proPalette(1:numel(vals), :);
bh.EdgeColor = [0 0 0];
bh.LineWidth = 1.8;

er = errorbar(1:numel(vals), vals, errs, 'k.', 'LineWidth', 1.4);
er.CapSize = 10;

xticks(1:numel(vals));
xticklabels(names);
xtickangle(30);
ylabel('Score');
ylim([0 max(vals+errs)*1.15]);
grid off; box off;
ax = gca;
ax.TickDir = 'out';        % little ticks
ax.TickLength = [0.02 0.02];

% title(sprintf('%s: CV summary (mean \\pm SEM)\nPrev=%.2f%%', behaviour, 100*mPrev));

%% ---------- FIGURE 2: PR curve ----------
fig2 = figure('Color','w','Position',[600 100 450 380]); hold on;

plot(rec, prec, 'LineWidth', 2.5, 'Color', 'r');
hBL = yline(prev_pooled, '--', sprintf('Baseline prevalence = %.3f', prev_pooled), ...
    'LabelHorizontalAlignment','left', FontSize=14);
hBL.LineWidth = 2;
try
    hBL.Label.FontWeight = 'bold';
    hBL.Label.FontSize = 13;
catch
    hBL.FontWeight = 'bold';
end

xlabel('Recall');
ylabel('Precision');
% title(sprintf('%s: PR curve (pooled), AUC=%.3f', behaviour, prAUC_pooled));
grid on; box off;
xlim([0 1]); ylim([0 1]);
ax = gca;
ax.TickDir = 'out';        % little ticks
ax.TickLength = [0.02 0.02];

%% ---------- FIGURE 3: Calibration ----------
fig3 = figure('Color','w','Position',[1100 100 450 380]); hold on;

plot([0 1],[0 1],'k--','LineWidth',1.5);
plot(binMeanP, binFracY, '-o', 'LineWidth', 2);

xlabel('Mean predicted probability (bin)');
ylabel('Observed fraction positive (bin)');
% title(sprintf('%s: Calibration (pooled), ECE=%.3f', behaviour, ece));
grid on; box off;
xlim([0 1]); ylim([0 1]);

dy = 0.03;      % default vertical offset
dx_first = 0.03; % push first label right
dy_last  = -0.04; % push last label downward
ax = gca;
ax.TickDir = 'out';        % little ticks
ax.TickLength = [0.02 0.02];

for bb = 1:numel(binMeanP)
    if ~isnan(binMeanP(bb)) && binN(bb) > 0

        x = binMeanP(bb);
        y = binFracY(bb) + dy;

        % Special-case adjustments
        if bb == 1
            x = x + dx_first;     % move right
            y = y+0.02;
        elseif bb == numel(binMeanP)
            y = binFracY(bb) + dy_last;  % move down
        end

        % Clamp to axes limits so it stays visible
        x = min(max(x, 0.02), 0.98);
        y = min(max(y, 0.02), 0.98);

        text(x, y, sprintf('n=%d', binN(bb)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontWeight','bold', ...
            'FontSize', 12);
    end
end


%%
% ---------- SAVE FIGURES (high-quality PNG) ----------
if do_save
    save_folder = "D:\PhD 3rd Year\Thesis Figures\Model performance metrics";

    prefix = lower(behaviour);
    dpi = 600;

    set([fig1 fig2 fig3], 'Color', 'w');

    exportgraphics(fig1, fullfile(save_folder, prefix + "_perf_summary.png"), 'Resolution', dpi);
    exportgraphics(fig2, fullfile(save_folder, prefix + "_pr_curve.png"),      'Resolution', dpi);
    exportgraphics(fig3, fullfile(save_folder, prefix + "_calibration.png"),   'Resolution', dpi);

    disp("Saved figures to: " + save_folder);
end
