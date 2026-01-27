sr_raw = 30000;
sr_lfp = 1000;
d_sample_factor = sr_raw/sr_lfp;

raw_path = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Concatenated Data\mouse8_renewal_concatenated_neural_data.dat";
kilosort_path = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Concatenated Data\kilosort4";
evt_path = "Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 8\Renewal\Neural Data\Triggers\mouse8_renewal_checkerboard_triggers.mat";
chan_map = readNPY(fullfile(kilosort_path,'channel_map.npy'));   % [nCh x 1]
chan_pos = readNPY(fullfile(kilosort_path,'channel_positions.npy')); % [nCh x 2]

num_channels = numel(chan_map);

% Simple lowpass for LFP (<150 Hz)
lfp_lowpass_filter = designfilt('lowpassiir', ...
               'FilterOrder', 4, ...
               'HalfPowerFrequency', 150, ...
               'SampleRate', sr_raw);

t_pre = 0.2;
t_post = 0.5;

n_pre_raw  = round(t_pre  * sr_raw);   % samples @ 30 kHz
n_post_raw = round(t_post * sr_raw);

window_raw = -n_pre_raw:n_post_raw;
num_timepoints_raw = numel(window_raw);

t_pre_lfp  = t_pre;
t_post_lfp = t_post;
n_pre_lfp  = round(t_pre_lfp  * sr_lfp);
n_post_lfp = round(t_post_lfp * sr_lfp);
num_timepoints_lfp = n_pre_lfp + n_post_lfp + 1;

%%

load(evt_path)
evt = evt(2:end);

file_id = fopen(raw_path, 'r');

lfp_sum_across_events = zeros(num_channels, num_timepoints_lfp);
number_of_events_used = 0;

n_events = numel(evt);
fprintf('Processing %d events...\n', n_events);

for event_idx = 1:numel(evt)
    current_evt_sample = evt(event_idx);
    
    % Calculate start and end sample numbers
    start_sample = current_evt_sample - n_pre_raw;
   
    if start_sample < 1, continue
    end

    end_sample = current_evt_sample + n_post_raw;
    
    % Number of raw samples in the window
    num_samples_window  = end_sample - start_sample + 1;
    
    byte_offset = (start_sample - 1) * num_channels * 2;

    status = fseek(file_id, byte_offset, "bof");

    if status ~= 0
        warning('fseek failed for event %d, skipping', event_idx);
            continue;
    end

    raw_data_chunk = fread(file_id, [num_channels, num_samples_window], "int16")';
    raw_data_chunk = double(raw_data_chunk);   % % [num_samples_window x num_channels]
    
    lfp_data_chunk = filtfilt(lfp_lowpass_filter, raw_data_chunk);
    
    lfp_data_chunk_ds = lfp_data_chunk(1:d_sample_factor:end, :).';  % now [channel x time]

    current_chunk_length = size(lfp_data_chunk_ds, 2);

    if current_chunk_length > num_timepoints_lfp
        % Trim extra samples
        lfp_data_chunk_ds = lfp_data_chunk_ds(:, 1:num_timepoints_lfp);
    
    elseif current_chunk_length < num_timepoints_lfp
        % Pad with NaNs
        padded_chunk = nan(num_channels, num_timepoints_lfp);
        padded_chunk(:, 1:current_chunk_length) = lfp_data_chunk_ds;
        lfp_data_chunk_ds = padded_chunk;
    end


    lfp_sum_across_events = lfp_sum_across_events + lfp_data_chunk_ds;
    number_of_events_used = number_of_events_used + 1;

    if mod(event_idx, 10) == 0 || event_idx == n_events
        fprintf('  Event %d / %d (%.1f%%)\n', ...
            event_idx, n_events, 100*event_idx/n_events);
    end


end

fclose(file_id);

lfp_avg_per_channel = lfp_sum_across_events / number_of_events_used;


%% Build time axis for LFP (in seconds)
t_lfp = (-n_pre_lfp : n_post_lfp) / sr_lfp;   % 1 x num_timepoints_lfp

% Sanity check dimensions
assert(size(lfp_avg_per_channel, 2) == numel(t_lfp), ...
    'Time axis and LFP data have mismatched lengths.');

%% Sort channels by depth (for all subsequent plots)
channel_depths = chan_pos(:, 2);              % y-coordinate in µm
[sorted_depths, sort_idx] = sort(channel_depths, 'ascend');

%% 1) Raw checkerboard-locked LFP heatmap (no baseline / CAR)
lfp_sorted_raw = lfp_avg_per_channel(sort_idx, :);   % [channels x time]

figure;
imagesc(t_lfp, sorted_depths, lfp_sorted_raw);
set(gca, 'YDir', 'normal');  % so depth increases downward in a sensible way
colormap('parula');
colorbar;
xlabel('Time (s)');
ylabel('Depth (\mum)');
title('Checkerboard-locked LFP (mean across events)');
hold on;
xline(0, 'k--', 'LineWidth', 1.5);  % vertical line at event time
hold off;

%% 2) Example single-channel trace
example_ch = 200;  % try a few: 50, 150, 250, 350, etc.
figure;
plot(t_lfp, lfp_avg_per_channel(example_ch, :));
xline(0, '--k');
xlabel('Time (s)');
ylabel('LFP (a.u.)');
title(sprintf('Channel %d: checkerboard-locked LFP', example_ch));

%% 3) Baseline subtraction per channel (-200 to 0 ms)
baseline_idx = t_lfp >= -0.2 & t_lfp < 0;          % -200 to 0 ms
baseline = mean(lfp_avg_per_channel(:, baseline_idx), 2, 'omitnan');  % [nCh x 1]

lfp_bs = lfp_avg_per_channel - baseline;           % [nCh x nTime]

lfp_bs_sorted = lfp_bs(sort_idx, :);

figure;
imagesc(t_lfp, sorted_depths, lfp_bs_sorted);
set(gca, 'YDir', 'normal');
colormap parula;
colorbar;
xlabel('Time (s)');
ylabel('Depth (\mum)');
title('Baseline-subtracted checkerboard-locked LFP');
xline(0, '--k');

max_abs = max(abs(lfp_bs_sorted(:)), [], 'omitnan');
caxis([-max_abs max_abs]);

%% 4) Spectrum of a single channel (optional debug)
spec_ch = 200;     % pick any channel
x_spec  = lfp_avg_per_channel(spec_ch, :);
Fs      = sr_lfp;  % 1000 Hz

N  = numel(x_spec);
f  = (0:N-1) * (Fs/N);
X  = abs(fft(x_spec));

figure;
plot(f(1:floor(N/2)), X(1:floor(N/2)));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(sprintf('Average LFP spectrum, channel %d', spec_ch));

%% 5) Common-average reference (CAR) on baseline-subtracted LFP
car_trace = mean(lfp_bs, 1, 'omitnan');      % [1 x nTime], average across channels
lfp_car   = lfp_bs - car_trace;              % [nCh x nTime]

lfp_car_sorted = lfp_car(sort_idx, :);

figure;
imagesc(t_lfp, sorted_depths, lfp_car_sorted);
set(gca, 'YDir', 'normal');
colormap parula;
colorbar;
xlabel('Time (s)');
ylabel('Depth (\mum)');
title('CAR + baseline-subtracted checkerboard-locked LFP');
xline(0, '--k');

max_abs = max(abs(lfp_car_sorted(:)), [], 'omitnan');
caxis([-max_abs max_abs]);

%% 6) CSD-ish (second spatial derivative along depth)
nCh = size(lfp_car_sorted, 1);
csd = nan(size(lfp_car_sorted));   % same size, NaN at edges

% V(z-1) - 2*V(z) + V(z+1)
csd(2:end-1, :) = ...
      lfp_car_sorted(1:end-2, :) ...
    - 2*lfp_car_sorted(2:end-1, :) ...
    + lfp_car_sorted(3:end, :);

% (Optional) convert to physical units with depth spacing:
% dz = median(diff(sorted_depths));
% csd = csd / (dz^2);

figure;
imagesc(t_lfp, sorted_depths, csd);
set(gca, 'YDir', 'normal');
colormap parula;
colorbar;
xlabel('Time (s)');
ylabel('Depth (\mum)');
title('CSD-ish (second spatial derivative) map');
xline(0, '--k');

max_abs = max(abs(csd(:)), [], 'omitnan');
caxis([-max_abs max_abs]);

%% 7) Response (0–200 ms) minus baseline (-200–0 ms) per channel (CAR LFP)
baseline_idx = t_lfp >= -0.2 & t_lfp < 0;
resp_idx     = t_lfp >= 0    & t_lfp <= 0.2;

baseline_mean = mean(lfp_car(:, baseline_idx), 2, 'omitnan');  % [nCh x 1]
resp_mean     = mean(lfp_car(:, resp_idx), 2, 'omitnan');      % [nCh x 1]
resp_diff     = resp_mean - baseline_mean;                     % [nCh x 1]

figure;
plot(resp_diff(sort_idx), sorted_depths, '-o');
set(gca, 'YDir', 'reverse');   % deepest at bottom if you prefer
xlabel('Response (post - pre) [a.u., CAR]');
ylabel('Depth (\mum)');
title('Checkerboard response per channel (CAR LFP)');

%%
baseline_idx = t_lfp >= -0.2 & t_lfp < 0;
resp_idx     = t_lfp >= 0    & t_lfp <= 0.2;

baseline_mean = mean(lfp_car(:, baseline_idx), 2, 'omitnan');
resp_mean     = mean(lfp_car(:, resp_idx), 2, 'omitnan');
resp_diff     = resp_mean - baseline_mean;

plot(resp_diff(sort_idx), sorted_depths, '-o');

%%
