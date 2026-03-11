function scrub3d_and_export(Xfit, edges, make_lines)
% scrub3d_and_export(X, edges, make_lines)
% X:     Np x 3 x Ntt (e.g., 11 x 3 x 10040)
% edges: Ne x 2 list of connections (indices into points)
% make_lines: true/false

if nargin < 3, make_lines = true; end

[Np, D, Ntt] = size(Xfit);
assert(D==3, 'X must be Np x 3 x Ntt');

% edges = [
%     1 2
%     1 3
%     2 3
%     4 5
%     4 6
%     5 6
%     7 2
%     7 3
%     7 8
%     8 9
%     9 10
%     10 11
% ];

% Safety: drop edges out of range
edges = edges(all(edges <= Np, 2), :);

% Auto-axis limits from all finite points
Xall = reshape(Xfit, [], 3);
good = all(isfinite(Xall), 2);
Xall = Xall(good, :);
if isempty(Xall), error('All X are NaN/Inf.'); end

mins = min(Xall, [], 1);
maxs = max(Xall, [], 1);
ctr  = (mins + maxs)/2;
span = maxs - mins;
pad  = 0.10 * max(span); if pad==0, pad=1; end

xL = [ctr(1)-span(1)/2-pad, ctr(1)+span(1)/2+pad];
yL = [ctr(2)-span(2)/2-pad, ctr(2)+span(2)/2+pad];
zL = [ctr(3)-span(3)/2-pad, ctr(3)+span(3)/2+pad];

% Figure
fig = figure('Color','w','Name','3D scrubber','NumberTitle','off');
ax = axes('Parent',fig);
grid(ax,'off'); axis(ax,'vis3d');
xlim(ax,xL); ylim(ax,yL); zlim(ax,zL);
view(ax,[-70 20]);

% Plot objects (create once, update each frame)
hold(ax,'on');
% One distinct color per marker:
% Colorblind-friendly palette (Okabe–Ito)
okabeIto = [
    0.5000 0.3000 0.0000  % black
    0.9020 0.6235 0.0000  % orange
    0.3373 0.7059 0.9137  % sky blue
    0.0000 0.6196 0.4510  % bluish green
    0.9412 0.8941 0.2588  % yellow
    0.0000 0.4471 0.6980  % blue
    0.8353 0.3686 0.0000  % vermillion
    0.8000 0.4745 0.6549  % reddish purple
    0.4000 0.4000 0.4000  % dark grey  (NEW)
    0.3922 0.0000 0.8000  % deep purple (NEW)
];


% Assign to markers (repeats if Np > 8)
markerColors = okabeIto(mod(0:Np-1, size(okabeIto,1)) + 1, :);

hPts = scatter3(ax, nan(Np,1), nan(Np,1), nan(Np,1), ...
    50, markerColors, 'filled');   % 80 = marker size (tweak)

% Pre-create line handles (optional)
hLn = gobjects(size(edges,1),1);
if make_lines
    for e = 1:size(edges,1)
        hLn(e) = line(ax, [nan nan], [nan nan], [nan nan], 'Color','k','LineWidth',2);
    end
end

% UI controls
uicontrol(fig,'Style','text','String','Frame','Units','normalized',...
    'Position',[0.02 0.01 0.05 0.04],'BackgroundColor','w');

hEdit = uicontrol(fig,'Style','edit','String','1','Units','normalized',...
    'Position',[0.08 0.01 0.08 0.05],'Callback',@jumpFrame);

hSlider = uicontrol(fig,'Style','slider','Min',1,'Max',Ntt,'Value',1,...
    'SliderStep',[1/(Ntt-1) 50/(Ntt-1)], 'Units','normalized',...
    'Position',[0.18 0.015 0.55 0.04],'Callback',@slideFrame);

hSave = uicontrol(fig,'Style','pushbutton','String','Save frame (HQ)',...
    'Units','normalized','Position',[0.75 0.01 0.22 0.05],...
    'Callback',@saveFrame);

% View buttons
uicontrol(fig,'Style','pushbutton','String','Front', ...
    'Units','normalized','Position',[0.18 0.06 0.10 0.05], ...
    'Callback',@(~,~) setView('Front'));

uicontrol(fig,'Style','pushbutton','String','Side', ...
    'Units','normalized','Position',[0.29 0.06 0.10 0.05], ...
    'Callback',@(~,~) setView('Side'));

uicontrol(fig,'Style','pushbutton','String','Top', ...
    'Units','normalized','Position',[0.40 0.06 0.10 0.05], ...
    'Callback',@(~,~) setView('Top'));

uicontrol(fig,'Style','pushbutton','String','Iso', ...
    'Units','normalized','Position',[0.51 0.06 0.10 0.05], ...
    'Callback',@(~,~) setView('Iso'));

% optional: show index in title
t = 1;
updateFrame(t);

% ---- callbacks ----
    function slideFrame(~,~)
        t = round(get(hSlider,'Value'));
        set(hSlider,'Value',t);
        set(hEdit,'String',num2str(t));
        updateFrame(t);
    end

    function jumpFrame(~,~)
        t = round(str2double(get(hEdit,'String')));
        if isnan(t), t = 1; end
        t = max(1, min(Ntt, t));
        set(hEdit,'String',num2str(t));
        set(hSlider,'Value',t);
        updateFrame(t);
    end

    function updateFrame(t)
       Xt = Xfit(:,:,t);
    set(hPts, 'XData', Xt(:,1), 'YData', Xt(:,2), 'ZData', Xt(:,3)); % edit back to add removed point after

        if make_lines
            for e = 1:size(edges,1)
                i = edges(e,1); j = edges(e,2);
                if all(isfinite(Xt([i j],:)), 'all')
                    set(hLn(e), 'XData',[Xt(i,1) Xt(j,1)], ...
                                'YData',[Xt(i,2) Xt(j,2)], ...
                                'ZData',[Xt(i,3) Xt(j,3)], ...
                                'Visible','on');
                else
                    set(hLn(e),'Visible','off');
                end
            end
        end

        title(ax, sprintf('Frame %d / %d', t, Ntt));
        drawnow;
    end

    function saveFrame(~,~)
        t = round(get(hSlider,'Value'));

        % Choose folder/name
        [file, path] = uiputfile({'*.png';'*.tif';'*.pdf';'*.svg'}, ...
            'Save current frame as...', sprintf('mouse3D_frame_%05d.png', t));
        if isequal(file,0), return; end
        out = fullfile(path,file);

        % Make export consistent
        set(fig,'InvertHardcopy','off');     % keep white background
        set(fig,'Renderer','opengl');        % good for 3D

        % Export: PNG/TIF at high DPI; PDF/SVG vector (note: 3D can be quirky in vector)
        [~,~,ext] = fileparts(out);
        ext = lower(ext);

        switch ext
            case {'.png','.tif'}
                exportgraphics(ax, out, 'Resolution', 600); % 600 dpi for figures
            case {'.pdf','.svg'}
                exportgraphics(ax, out, 'ContentType','vector'); % may rasterize some 3D elements depending on renderer
            otherwise
                exportgraphics(ax, out, 'Resolution', 600);
        end

        fprintf('Saved: %s\n', out);
    end

function setView(name)
    switch lower(name)
        case 'front'
            view(ax, -60, 25);
        case 'side'
            view(ax,  30, 25);
        case 'top'
            view(ax, 120, 25); 
        case 'iso'
            view(ax, 210, 25); 
        otherwise
            warning('Unknown view "%s"', name);
    end
    drawnow;
end

end

