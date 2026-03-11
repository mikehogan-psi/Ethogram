function plot3d_video_test(Xfit, make_lines)
% X: Np x 3 x Ntt

[Np, D, Ntt] = size(Xfit);
assert(D == 3, 'X must be Np x 3 x Ntt');

% Define edges (no index > Np)
edges = [
    1 2
    1 3
    2 3
    4 5
    4 6
    5 6
    7 2
    7 3
    7 8
    8 9
    9 10
    10 11
];
edges = edges(all(edges <= Np, 2), :);

% ---- Auto axis limits from ALL data ----
Xall = reshape(Xfit, [], 3);                 % (Np*Ntt) x 3
good = all(isfinite(Xall), 2);
Xall = Xall(good, :);

if isempty(Xall)
    error('All points are NaN/Inf. Nothing to plot.');
end

mins = min(Xall, [], 1);
maxs = max(Xall, [], 1);
ctr  = (mins + maxs) / 2;
span = (maxs - mins);
pad  = 0.10 * max(span);                  % 10% padding based on max span
if pad == 0, pad = 1; end                 % handle degenerate case

xL = [ctr(1)-span(1)/2-pad, ctr(1)+span(1)/2+pad];
yL = [ctr(2)-span(2)/2-pad, ctr(2)+span(2)/2+pad];
zL = [ctr(3)-span(3)/2-pad, ctr(3)+span(3)/2+pad];

figure('Color','w');

for t = 1:Ntt
    cla;

    Xt = Xfit(:,:,t); % Np x 3

    % If this frame is completely invalid, skip it
    if ~any(isfinite(Xt(:)))
        title(sprintf('t=%d (all NaN/Inf)', t));
        drawnow;
        continue;
    end

    % Plot points
    plot3(Xt(:,1), Xt(:,2), Xt(:,3), '.', 'MarkerSize', 20);
    hold on;

    % Plot lines
    if make_lines
        for e = 1:size(edges,1)
            i = edges(e,1); j = edges(e,2);
            if all(isfinite(Xt([i j],:)), 'all')
                line([Xt(i,1) Xt(j,1)], [Xt(i,2) Xt(j,2)], [Xt(i,3) Xt(j,3)], ...
                    'Color','k','LineWidth',2);
            end
        end
    end

    xlim(xL); ylim(yL); zlim(zL);
    axis vis3d;
    grid on;
    view([-70 20]);
    title(sprintf('t=%d', t));
    grid("on")

    drawnow;
    pause(0.04);
end
end