function[evt_loom,evt_flash] = sort_evts_clean(evt, n_trials, onset, which_epoch, which_stim_set)
% SORT_EVTS - Extract loom and flash event times from trigger signal
%
% INPUTS:
%   evt           - vector of event timestamps (e.g. in seconds)
%   n_trials      - number of stimulus trials in this part (usually 20)
%   onset         - frame index (within a trial) corresponding to stimulus onset
%   which_epoch   - part index (2 = p1, 3 = p2)
%   which_stim_set- which stimulus set this mouse received (1 or 2)
%
% OUTPUTS:
%   evt_loom      - vector of loom stimulus onset times (1 per trial)
%   evt_flash     - vector of flash stimulus onset times (1 per trial)
% 
% This function is called in:
%   neural_data_preprocessing_pipeline.m



% Reshape the event vector into a [n_trials × trial_length] matrix
% Each row = one trial, each column = one sample
trial_length = numel(evt)/n_trials;
evt_mat = zeros(n_trials,trial_length);
for trial = 1:n_trials
    evt_mat(trial,:) = evt((trial-1)*trial_length+1:trial*trial_length);
end

if which_stim_set == 1
    if which_epoch == 2
       stim = [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1];
    elseif which_epoch == 3
       stim = [0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1];
    end

    elseif which_stim_set == 2
        if which_epoch == 2
        stim = [0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0];
        elseif which_epoch == 3
        stim = [1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0 ,0 ,1 ,0];
        end
end

evt_loom = evt_mat(stim==1, round(onset));
evt_flash = evt_mat(stim==0, round(onset));

end

