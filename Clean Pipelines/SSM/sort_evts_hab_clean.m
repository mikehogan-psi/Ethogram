function evt_hab = sort_evts_hab_clean(evt, n_events)
% SORT_EVTS_HAB - Generate pseudo-events for continuous habituation session
%
% INPUTS:
%   evt      - vector of event timestamps in SECONDS
%   n_events - number of evenly spaced events across habituation
%
% OUTPUT:
%   evt_hab  - vector of pseudo-event timestamps in SECONDS
% 
% This function is called in:
%   neural_data_preprocessing_pipeline.m

% Total duration of habituation
total_time = evt(end) - evt(1);

% Generate n_events evenly spaced times across the session
evt_hab = linspace(evt(1), evt(end), n_events)';

end