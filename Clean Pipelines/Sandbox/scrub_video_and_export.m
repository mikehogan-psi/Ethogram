function scrub_video_and_export(vidPath)
% scrub_video_and_export(vidPath)
% Scrub through a video and export selected frames in high quality.
%
% Example:
%   scrub_video_and_export("D:\data\vid.avi")

if nargin < 1 || isempty(vidPath)
    [f,p] = uigetfile({'*.avi;*.mp4;*.mov;*.m4v','Video files'; '*.*','All files'}, ...
                      'Select a video');
    if isequal(f,0), return; end
    vidPath = fullfile(p,f);
end

vr = VideoReader(vidPath);

% Some codecs make NumFrames unreliable; compute from duration.
fps = vr.FrameRate;
nFrames = max(1, floor(vr.Duration * fps));  % good enough for scrubbing
dur = vr.Duration;

% Read first frame to init display
vr.CurrentTime = 0;
frame = readFrame(vr);

fig = figure('Color','w','Name','Video scrubber','NumberTitle','off');
ax  = axes('Parent',fig);
hIm = imshow(frame, 'Parent', ax);
axis(ax,'image'); axis(ax,'off');
title(ax, sprintf('Frame 1 / ~%d   (t = %.3f s)', nFrames, 0));

% --- UI controls ---
uicontrol(fig,'Style','text','String','Frame','Units','normalized', ...
    'Position',[0.01 0.01 0.06 0.04],'BackgroundColor','w');

hEditFrame = uicontrol(fig,'Style','edit','String','1','Units','normalized', ...
    'Position',[0.07 0.01 0.08 0.05],'Callback',@jumpToFrame);

uicontrol(fig,'Style','text','String','Time (s)','Units','normalized', ...
    'Position',[0.16 0.01 0.07 0.04],'BackgroundColor','w');

hEditTime = uicontrol(fig,'Style','edit','String','0','Units','normalized', ...
    'Position',[0.23 0.01 0.08 0.05],'Callback',@jumpToTime);

hSlider = uicontrol(fig,'Style','slider','Min',1,'Max',nFrames,'Value',1, ...
    'SliderStep',[1/max(1,nFrames-1) 50/max(1,nFrames-1)], ...
    'Units','normalized','Position',[0.33 0.015 0.42 0.04], ...
    'Callback',@slideToFrame);

hSave = uicontrol(fig,'Style','pushbutton','String','Save frame (HQ)', ...
    'Units','normalized','Position',[0.76 0.01 0.23 0.05], ...
    'Callback',@saveFrame);

% Cache last good frame/time so the UI feels snappy
currentFrameIdx = 1;
currentTimeSec  = 0;
currentImg      = frame;

% -------- nested helpers ----------
    function setFrameByIndex(idx)
        % Clamp
        idx = round(idx);
        idx = max(1, min(nFrames, idx));

        % Convert frame index -> time.
        % Using (idx-1)/fps aligns frame 1 at t=0.
        t = (idx-1) / fps;
        t = max(0, min(dur - 1e-6, t)); % avoid exactly Duration

        % Seek and read. readFrame reads the NEXT frame from CurrentTime.
        % Set CurrentTime then readFrame to approximate frame at time t.
        vr.CurrentTime = t;
        try
            fr = readFrame(vr);
        catch
            % Fallback: if readFrame errors at end, step slightly back
            vr.CurrentTime = max(0, t - (1/fps));
            fr = readFrame(vr);
        end

        currentFrameIdx = idx;
        currentTimeSec  = t;
        currentImg      = fr;

        % Update display + UI
        set(hIm, 'CData', fr);
        set(hEditFrame, 'String', num2str(idx));
        set(hEditTime,  'String', sprintf('%.3f', t));
        set(hSlider,    'Value', idx);

        title(ax, sprintf('Frame %d / ~%d   (t = %.3f s)', idx, nFrames, t));
        drawnow;
    end

    function slideToFrame(~,~)
        setFrameByIndex(get(hSlider,'Value'));
    end

    function jumpToFrame(~,~)
        idx = str2double(get(hEditFrame,'String'));
        if isnan(idx), idx = currentFrameIdx; end
        setFrameByIndex(idx);
    end

    function jumpToTime(~,~)
        t = str2double(get(hEditTime,'String'));
        if isnan(t), t = currentTimeSec; end
        t = max(0, min(dur - 1e-6, t));
        idx = round(t * fps) + 1;
        setFrameByIndex(idx);
    end

    function saveFrame(~,~)
        idx = currentFrameIdx;

        % Ask where to save
        [file, path] = uiputfile({'*.png';'*.tif';'*.jpg'}, ...
            'Save current frame as...', sprintf('frame_%05d.png', idx));
        if isequal(file,0), return; end
        out = fullfile(path,file);

        % Save at ORIGINAL pixel resolution (best for images)
        % This avoids figure scaling / DPI weirdness.
        [~,~,ext] = fileparts(out);
        ext = lower(ext);
        switch ext
            case '.png'
                imwrite(currentImg, out, 'png');
            case '.tif'
                imwrite(currentImg, out, 'tif');
            case '.jpg'
                imwrite(currentImg, out, 'jpg', 'Quality', 98);
            otherwise
                imwrite(currentImg, out, 'png');
        end

        fprintf('Saved frame %d to: %s\n', idx, out);
    end

% Init
setFrameByIndex(1);

end