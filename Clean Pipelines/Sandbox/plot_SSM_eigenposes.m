% SSM_save_path = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Models\SSMs\SSM_3D_implant_mouse1_head_fixed.mat";
SSM_save_path = "W:\Mike\Neuropixels_Fear_Conditioning\Models\SSMs\SSM_mouse1_to_9.mat";

S = load(SSM_save_path);
template = S.template;     % Np x 3
eignV3D  = S.eignV3D;      % Np x 3 x Ndim
lambda   = S.lambda;       % Ndim x 1

% ---- settings ----
k = 2;                  % amplitude in SD units
fps = 30;               % animation framerate
seconds_per_pose = 6;   % show each eigenpose for this long
nCycles = 2;            % sine cycles per pose (1 = back-and-forth once)
viewAzEl = [-60 25];    % camera view

[Np,~,Ndim] = size(eignV3D);

% Fixed axis limits (so each pose is framed consistently)
mins = min(template,[],1);
maxs = max(template,[],1);
pad  = 0.25 * max(maxs - mins);
xL = [mins(1)-pad maxs(1)+pad];
yL = [mins(2)-pad maxs(2)+pad];
zL = [mins(3)-pad maxs(3)+pad];

nFramesPose = max(2, round(seconds_per_pose * fps));
pattern = sin(linspace(0, 2*pi*nCycles, nFramesPose));

for d = 1:Ndim
    % New figure for each eigenpose
    fig = figure('Color','w');
    ax = axes('Parent',fig); hold(ax,'on');
    axis(ax,'equal'); axis(ax,'vis3d');
    grid(ax,'on'); box(ax,'on');
    camproj(ax,'orthographic');
    view(ax, viewAzEl(1), viewAzEl(2));
    xlim(ax,xL); ylim(ax,yL); zlim(ax,zL);
    xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');

    % Single handle: only ONE pose shown at a time
    h = plot3(ax, template(:,1), template(:,2), template(:,3), ...
        'o', 'MarkerSize', 7, 'LineWidth', 1.2);

    for t = 1:nFramesPose
        a = k * pattern(t);
        Xfig = template + a*sqrt(lambda(d)) * eignV3D(:,:,d);

        set(h, 'XData', Xfig(:,1), 'YData', Xfig(:,2), 'ZData', Xfig(:,3));
        title(ax, sprintf('Eigenpose %d / %d   (%.2f SD)', d, Ndim, a));
        drawnow;
        pause(1/fps);
    end

    % Close before moving to the next eigenpose
    close(fig);
end