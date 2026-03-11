% --- settings ---
maxPlotFrames = 4000;
alphaCloud    = 0.1;
ptSizeCloud   = 8;
ptSizeMean    = 80;
zThresh       = 8;      % robust z cutoff (6–10 typical)

[Np,~,Nt] = size(X);

% pick frames
if Nt > maxPlotFrames
    idxPlot = randperm(Nt, maxPlotFrames);
else
    idxPlot = 1:Nt;
end

Xplot = X(:,:,idxPlot);     % Np x 3 x Nplot

% ---- robust centre/scale from ALL points in Xplot ----
Xp_all = reshape(Xplot, [], 3);
Xp_all = Xp_all(all(isfinite(Xp_all),2),:);

mu  = median(Xp_all, 1);
mad = median(abs(Xp_all - mu), 1);
mad(mad==0) = 1e-6;          % avoid divide-by-zero

% ---- create filtered version with outliers set to NaN (for RANSAC + plotting) ----
Xplot_f = Xplot;
for m = 1:Np
    Xm = squeeze(Xplot_f(m,:,:))';                         % Nplot x 3
    good = all(isfinite(Xm),2);
    z_m = abs((Xm - mu) ./ (1.4826*mad));                  % robust z
    bad = good & any(z_m >= zThresh, 2);                   % only mark finite rows
    Xm(bad,:) = NaN;
    Xplot_f(m,:,:) = permute(Xm, [2 1]);                   % back to 1x3xN
end

% Robust centroid used by your pipeline
template_ransac = Estimate_mean_RANSAC(Xplot_f, false);     % Np x 3

% ---- colours ----
if Np == 11
    markerColors = [
        0.5000 0.3000 0.0000
        0.9020 0.6235 0.0000
        0.3373 0.7059 0.9137
        0.0000 0.6196 0.4510
        0.9412 0.8941 0.2588
        0.0000 0.4471 0.6980
        0.8353 0.3686 0.0000
        0.8000 0.4745 0.6549
        0.4000 0.4000 0.4000
        0.3922 0.0000 0.8000
        0.6000 0.2000 0.1000
        0.7631 0.1342 0.5321
    ];
else
    markerColors = hsv2rgb([linspace(0,1,Np)' 0.85*ones(Np,1) 0.95*ones(Np,1)]);
end

% ---- plot ----
fig = figure('Color','w','Renderer','opengl');
ax = axes('Parent',fig); hold(ax,'on');
set(fig,'Units','pixels','Position',[100 100 1300 800]);   % bigger window
set(ax,'Units','normalized','Position',[0.06 0.08 0.72 0.86]); % leave room right for legend + margins

% clouds (use FILTERED Xplot_f so outliers don't wreck the view)
for m = 1:Np
    Xm = squeeze(Xplot_f(m,:,:))';                         % Nplot x 3
    Xm = Xm(all(isfinite(Xm),2),:);
    if isempty(Xm), continue; end

    scatter3(ax, Xm(:,1), Xm(:,2), Xm(:,3), ptSizeCloud, ...
        'filled', 'MarkerFaceColor', markerColors(m,:), ...
        'MarkerEdgeColor','none', 'MarkerFaceAlpha', alphaCloud);
end

% robust mean pose points (RANSAC centroid) - nicer aesthetics
for m = 1:Np
    ctr = template_ransac(m,:);
    if any(~isfinite(ctr)), continue; end

    % 1) halo (under)
    scatter3(ax, ctr(1), ctr(2), ctr(3), ptSizeMean*1.35, ...
        'o', 'filled', ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', [1 1 1], ...
        'LineWidth', 1.0);

    % 2) coloured centroid (over)
    scatter3(ax, ctr(1), ctr(2), ctr(3), ptSizeMean, ...
        'o', 'filled', ...
        'MarkerFaceColor', markerColors(m,:), ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineWidth', 0.8);
end

% axis limits from FILTERED points
Xp = reshape(Xplot_f, [], 3);
Xp = Xp(all(isfinite(Xp),2),:);
lo = prctile(Xp, 1, 1);
hi = prctile(Xp, 99, 1);
pad = 0.10 * max(hi - lo);

xlim(ax, [lo(1)-pad hi(1)+pad]);
ylim(ax, [lo(2)-pad hi(2)+pad]);
zlim(ax, [lo(3)-pad hi(3)+pad]);

axis(ax,'equal'); axis(ax,'vis3d');
camproj(ax,'orthographic');
camup(ax,[0 0 1]);
camtarget(ax, [mean(xlim(ax)) mean(ylim(ax)) mean(zlim(ax))]);
camzoom(ax, 1.2);
view(ax, -60, 25);

% --- re-frame camera so the whole cloud stays in view ---
axis(ax,'equal'); axis(ax,'vis3d');
camproj(ax,'orthographic');

% centre camera on the middle of your axis limits
camtarget(ax, [mean(xlim(ax)) mean(ylim(ax)) mean(zlim(ax))]);

% zoom OUT a bit (values <1 zoom out)
camzoom(ax, 0.85);     % try 0.8–0.9

% optional: widen limits a little more (extra padding)
xl = xlim(ax); yl = ylim(ax); zl = zlim(ax);
pad2 = 0.15;  % 15% extra
xlim(ax, xl + pad2*[-1 1]*range(xl));
ylim(ax, yl + pad2*[-1 1]*range(yl));
zlim(ax, zl + pad2*[-1 1]*range(zl));

labels = { ...
    'Cable Tip', ...
    'Box Left Back Corner', ...
    'Box Right Back Corner', ...
    'Nose', ...
    'Left Ear', ...
    'Right Ear', ...
    'Neck Base', ...
    'Body Anterior', ...
    'Body Posterior', ...
    'Tail Base' ...
    'Tail Anterior'...
};

% Create invisible dummy handles for a clean legend
hLeg = gobjects(Np,1);
for m = 1:Np
    hLeg(m) = scatter3(ax, NaN, NaN, NaN, ptSizeMean, ...
        'filled', 'MarkerFaceColor', markerColors(m,:), ...
        'MarkerEdgeColor','k', 'LineWidth', 0.6);
end

% Make legend markers bigger
for m = 1:Np
    hLeg(m).SizeData = 200;      % try 150–400 (area in points^2)
    hLeg(m).LineWidth = 1.0;     % thicker edge in legend
end

set(fig,'Units','pixels','Position',[100 100 1300 800]);
set(ax,'Units','normalized','Position',[0.06 0.08 0.72 0.86]);

lgd = legend(ax, hLeg, labels, 'Location','northeast', 'Box','off');
lgd.FontSize = 20;
lgd.FontWeight = 'bold';
lgd.ItemTokenSize = [18 12];
grid(ax,'on'); box(ax,'on');
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
set(hLeg(m), 'SizeData', 120);  % legend marker size (area in points^2)

% --- prevent clipping of 3D axes/labels ---
set(ax,'Units','normalized');

% Give MATLAB more padding around axes content
ax.LooseInset = max(ax.LooseInset, [0.08 0.08 0.08 0.08]);   % [left bottom right top]

% Make sure the axes outer box fits inside the figure nicely
ax.ZLabel.String = 'Z';
ax.ZLabel.Units = 'normalized';
ax.ZLabel.Position = [-0.08 0.5 0];   % move left/inward; tweak -0.08 to -0.12 if needed
set(ax,'Units','normalized');
ax.Position = [0.18 0.10 0.52 0.82];   % more left margin
ax.FontSize   = 14;        % tick label size
ax.FontWeight = 'bold';    % tick labels bold
ax.XLabel.FontSize = 16; ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 16; ax.YLabel.FontWeight = 'bold';
ax.ZLabel.FontSize = 16; ax.ZLabel.FontWeight = 'bold';
ax.LineWidth = 1.8;   % thickness of the box/axes lines