function [spike_counts, bin_edges, bin_centres] = raster_binned(tsp, evt, t_start, t_end, num_bins)
% RASTER_BINNED: Bin spike trains around events into fixed trial × bin counts
%
% Inputs:
%   tsp      : spike timestamps (in seconds)
%   evt      : event times (in seconds)
%   t_start  : window start relative to event (e.g., -10)
%   t_end    : window end relative to event (e.g., 20)
%   num_bins : number of bins for the window
%
% Outputs:
%   spike_counts : [nTrials × num_bins] spike counts
%   bin_edges    : edges of time bins relative to event
%   bin_centres  : centres of time bins relative to event

% Define bin edges
bin_edges   = linspace(t_start, t_end, num_bins+1);
bin_centres = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

nTrials = numel(evt);
spike_counts = zeros(nTrials, num_bins);

for trial = 1:nTrials
    % Get spike times relative to event
    rel_spikes = tsp - evt(trial);

    % Keep only spikes inside window
    rel_spikes = rel_spikes(rel_spikes >= t_start & rel_spikes <= t_end);

    % Histogram spike counts into bins
    spike_counts(trial, :) = histcounts(rel_spikes, bin_edges);
end

end
