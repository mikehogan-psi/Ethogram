function plotPseudoR2Hist(pseudoR2, binEdges, nCells)
% plotPseudoR2Hist  Plot % of cells vs pseudoR2 with adjustable bin edges.
% 
% Usage:
%   plotPseudoR2Hist(pseudoR2, binEdges)            % nCells defaults to 1774
%   plotPseudoR2Hist(pseudoR2, binEdges, nCells)    % set your own total
%
% Inputs:
%   pseudoR2 : vector of values (NaNs ignored)
%   binEdges : vector of bin edges (e.g., 0:0.05:1 or linspace(0,1,21))
%   nCells   : scalar, total cells to divide by (default 1774)
%
% The y-axis is % of cells = (count / nCells) * 100.

    if nargin < 3 || isempty(nCells), nCells = 1774; end
    if nargin < 2 || isempty(binEdges)
        % Default: 20 bins spanning the finite data range
        x = pseudoR2(:);
        x = x(isfinite(x));
        if isempty(x), error('No finite pseudoR2 values.'); end
        binEdges = linspace(min(x), max(x), 21);
    end

    % Keep finite values only
    x = pseudoR2(:);
    x = x(isfinite(x));

    % Histogram counts using specified edges
    counts = histcounts(x, binEdges);

    % Convert to % of cells (divide by nCells, not by #valid)
    pct = (counts / nCells) * 100;

    % Bar plot at bin centers, with bar widths matching bin widths
    centers = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    widths  = diff(binEdges);

    bar(centers, pct, 'BarWidth', 1, 'EdgeColor', 'k', 'FaceColor', [0.7 0.7 0.7]);
    set(gca, 'XLim', [binEdges(1), binEdges(end)]);
    xlabel('pseudoR2');
    ylabel('% of cells');
    box on;
    xlim([0.05 1])
end
