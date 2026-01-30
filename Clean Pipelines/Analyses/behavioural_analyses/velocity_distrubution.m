v = all_velocity_data(:);

% Avoid log(0)
v_log = log10(v);

%% Plot density
figure;
histogram(v_log, 200, 'Normalization', 'pdf');   % 200 bins, adjust as needed
xlabel('log_{10}(velocity)');
ylabel('Probability density');
title('Velocity distribution (all mice)');

grid on;

%%
v = all_velocity_data(:);
%% Flatten and log-transform
v = v(~isnan(v) & v > 0);          % remove NaNs and zeros
v_log = log10(v);

% Restrict to central range to avoid crazy tails
x_grid = linspace(prctile(v_log,1), prctile(v_log,99), 512);
[f,xi] = ksdensity(v_log, x_grid);

%% Find peaks (local maxima) in the density
[pk, locs_pk] = findpeaks(f, xi);

if numel(pk) < 2
    warning('Only one clear peak found in velocity distribution; threshold will be heuristic.');
    % Fallback: e.g. use a low percentile
    log_thresh = prctile(v_log, 10);
else
    % Take the two highest peaks
    [~, order] = sort(pk, 'descend');
    peak_x = locs_pk(order(1:2));
    
    % Sort them so left = freezing, right = movement
    peak_x = sort(peak_x);
    x1 = peak_x(1);
    x2 = peak_x(2);

    % Find minimum of f between the two peaks
    mask = xi >= x1 & xi <= x2;
    [~, idx_min_local] = min(f(mask));
    xi_local = xi(mask);
    log_thresh = xi_local(idx_min_local);
end

velocity_threshold = 10.^log_thresh;

fprintf('Log threshold: %.3f\n', log_thresh);
fprintf('Velocity threshold: %.5f cm/frame\n', velocity_threshold);

%% Plot to check
figure; hold on;
plot(xi, f, 'LineWidth', 2);
yL = ylim;
plot([log_thresh log_thresh], yL, 'r--', 'LineWidth', 1.5);
xlabel('log_{10}(velocity)');
ylabel('Density');
title('Velocity KDE with chosen threshold');
grid on;
