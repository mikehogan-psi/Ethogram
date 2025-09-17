
function [] = extract_neural_data_triggers(hab_path, p1_path, p2_path,...
    check_path, hab_path_triggers, p1_path_triggers, p2_path_triggers,...
    check_path_triggers, save_path_triggers, save_name_triggers)
%EXTRACT_NEURAL_DATA_TRIGGERS Extract trigger event timings from Neuropixels recordings.
%   reads trigger event timings from .npy files produced during Neuropixels recordings and saves them 
%   as .mat files. The function aligns triggers across sessions and 
%   handles habituation, extinction parts 1 and 2, and optionally a 
%   checkerboard recording if present (for renewal).
%
%   The function:
%     • Loads trigger sample numbers and states from the *_triggers folders.
%     • Aligns events across consecutive recordings using cumulative offsets.
%     • Selects only rising-edge triggers (states == 1).
%     • Saves trigger events as .mat files in SAVE_PATH_TRIGGERS, with names
%       prefixed by SAVE_NAME_TRIGGERS.
%     • Skips extraction if a corresponding .mat file already exists.
%
%   Inputs:
%       HAB_PATH            - Path to habituation recording
%       P1_PATH             - Path to extinction part 1 recording
%       P2_PATH             - Path to extinction part 2 recording
%       CHECK_PATH          - Path to checkerboard recording (empty if none)
%       HAB_PATH_TRIGGERS   - Path to habituation triggers folder
%       P1_PATH_TRIGGERS    - Path to part 1 triggers folder
%       P2_PATH_TRIGGERS    - Path to part 2 triggers folder
%       CHECK_PATH_TRIGGERS - Path to checkerboard triggers folder (empty if none)
%       SAVE_PATH_TRIGGERS  - Directory where .mat files will be saved
%       SAVE_NAME_TRIGGERS  - Prefix for saved .mat file names
%
%   Outputs:
%       None (trigger times are saved to disk as .mat files).
%
%   To call this function with all inputs preformatted, use
%   neural_data_preprocessing_pipeline.mat

% Check whether checkerboard recording is present
if isempty(check_path) && isempty(check_path_triggers)

    cont0 = 0;
    
    % Extract trigger timings for habituation
    evt = readNPY(fullfile(hab_path_triggers, 'sample_numbers.npy')); % change rest to this format
    states = readNPY([hab_path_triggers 'states.npy']);
    cont = readNPY([hab_path 'sample_numbers.npy']);
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);
    if exist(fullfile(save_path_triggers, [save_name_triggers, 'habituation_triggers.mat']), 'file')
        warning([save_name_triggers, 'habituation_triggers already exists, skipping extraction'])
    else
        save(fullfile(save_path_triggers, [save_name_triggers, 'habituation_triggers.mat']),'evt');
        disp([save_name_triggers, 'habituation_triggers has been saved!'])
    end
    
    % Extract trigger timings for p1
    evt = readNPY([p1_path_triggers 'sample_numbers.npy']); 
    states = readNPY([p1_path_triggers 'states.npy']);
    cont = readNPY([p1_path 'sample_numbers.npy']);
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);
    if exist(fullfile(save_path_triggers, [save_name_triggers, 'p1_triggers.mat']), 'file')
        warning([save_name_triggers, 'p1_triggers already exists, skipping extraction'])
    else
        save(fullfile(save_path_triggers, [save_name_triggers, 'p1_triggers.mat']),'evt');
        disp([save_name_triggers, 'p1_triggers has been saved!'])
    end

    % Extract trigger timings for p2
    evt = readNPY([p2_path_triggers 'sample_numbers.npy']); 
    states = readNPY([p2_path_triggers 'states.npy']);
    cont = readNPY([p2_path 'sample_numbers.npy']);
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);
    if exist(fullfile(save_path_triggers, [save_name_triggers, 'p2_triggers.mat']), 'file')
        warning([save_name_triggers, 'p2_triggers already exists, skipping extraction'])
    else
        save(fullfile(save_path_triggers, [save_name_triggers, 'p2_triggers.mat']),'evt');
        disp([save_name_triggers, 'p2_triggers has been saved!'])
    end

else
    cont0 = 0;
    
    % Extract trigger timings for habituation
    evt = readNPY([hab_path_triggers 'sample_numbers.npy']); 
    states = readNPY([hab_path_triggers 'states.npy']);
    cont = readNPY([hab_path 'sample_numbers.npy']);
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);
    if exist(fullfile(save_path_triggers, [save_name_triggers, 'habituation_triggers.mat']), 'file')
        warning([save_name_triggers, 'habituation_triggers already exists, skipping extraction'])
    else
        save(fullfile(save_path_triggers, [save_name_triggers, 'habituation_triggers.mat']),'evt');
        disp([save_name_triggers, 'habituation_triggers has been saved!'])
    end
    
    % Extract trigger timings for p1
    evt = readNPY([p1_path_triggers 'sample_numbers.npy']); 
    states = readNPY([p1_path_triggers 'states.npy']);
    cont = readNPY([p1_path 'sample_numbers.npy']);
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);
    if exist(fullfile(save_path_triggers, [save_name_triggers, 'p1_triggers.mat']), 'file')
        warning([save_name_triggers, 'p1_triggers already exists, skipping extraction'])
    else
        save(fullfile(save_path_triggers, [save_name_triggers, 'p1_triggers.mat']),'evt');
        disp([save_name_triggers, 'p1_triggers has been saved!'])
    end
    
    % Extract trigger timings for p2
    evt = readNPY([p2_path_triggers 'sample_numbers.npy']); 
    states = readNPY([p2_path_triggers 'states.npy']);
    cont = readNPY([p2_path 'sample_numbers.npy']);
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);
    if exist(fullfile(save_path_triggers, [save_name_triggers, 'p2_triggers.mat']), 'file')
        warning([save_name_triggers, 'p2_triggers already exists, skipping extraction'])
    else
        save(fullfile(save_path_triggers, [save_name_triggers, 'p2_triggers.mat']),'evt');
        disp([save_name_triggers, 'p2_triggers has been saved!'])
    end

    % Extract trigger timings for checkerboard
    evt = readNPY([check_path_triggers 'sample_numbers.npy']); 
    states = readNPY([check_path_triggers 'states.npy']);
    cont = readNPY([check_path 'sample_numbers.npy']);
    evt = evt-cont(1)+cont0;
    evt = evt(states==1);
    cont0 = cont0+numel(cont);
    if exist(fullfile(save_path_triggers, [save_name_triggers, 'checkerboard_triggers.mat']), 'file')
        warning([save_name_triggers, 'checkerboard_triggers already exists, skipping extraction'])
    else
        save(fullfile(save_path_triggers, [save_name_triggers, 'checkerboard_triggers.mat']),'evt');
        disp([save_name_triggers, 'checkerboard_triggers has been saved!'])
    end

end