% NEURAL_DATA_PREPROCESSING_PIPELINE
% Concatenates raw .dat neural data files and saves in appropriate file
% using function (concatenate_neural_data_clean)
% location on shared drive
% Then extracts trigger timings in raw sampling rate terms (30KHz) using 
% extract_neural_data_triggers
% Then extracts event timings (redundant code now, do not use)

%% Setup 

% !!! Specify which session is to be analysed !!!
% Acquisition, Extinction, or Renewal
session = 'Renewal';

% !!! Provide directory with all data in !!!
master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

%!!!Specify stim set for each mouse (mouse number only)!!!
received_stim_set_1 = [1 2 3 6 7];
received_stim_set_2 = [4 5 8 9];
%% Get directories for .dat files and trigger.npy files and concantenate/extract
% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

% Get master neural data folders for specified session
for mouse = 2:length(mouse_files) % Always skip mouse 1 (data wrong format)
    mouse_name = mouse_files(mouse).name;
    mouse_path = mouse_files(mouse).folder;
    current_raw_neural_data_folder = fullfile(mouse_path, mouse_name,...
        session, 'Neural Data', 'Raw Data');
    
    % Construct list of filepaths for raw neural data folders for current mouse
    raw_folder_dir = dir(fullfile(current_raw_neural_data_folder, 'mouse*'));
    raw_folder_list = cell(length(raw_folder_dir), 1);
      
    for part = 1:length(raw_folder_dir)
        raw_folder_list{part} = fullfile(raw_folder_dir(part).folder,...
            raw_folder_dir(part).name);
    end
    
    % Sort parts into correct order prior to concatenation
    if strcmp(session, 'Extinction')
        idx = zeros(3,1);
        key_words = {'habituation', 'p1', 'p2'};
        for part = 1:numel(key_words)
            idx(contains(raw_folder_list, key_words{part})) = part;
        end
    elseif strcmp(session, 'Renewal')
        idx = zeros(4,1);
        key_words = {'habituation', 'p1', 'p2', 'checkerboard'};
        for part = 1:numel(key_words)
            idx(contains(raw_folder_list, key_words{part})) = part;
        end
    end

    [~, sort_idx] = sort(idx);
    raw_folder_list = raw_folder_list(sort_idx);


    % Assign filepaths for .dat files in each folder to a cell array
    dat_files = cell(length(raw_folder_list), 1);

    for part = 1:length(raw_folder_list)
        % Handle inconsistent experiment and recording numbers
        exp_dir = dir(fullfile(raw_folder_list{part}, 'Record Node 101', 'experiment*'));
        rec_dir = dir(fullfile(exp_dir(1).folder, exp_dir(1).name, 'recording*'));
        current_dat = fullfile(rec_dir(1).folder, rec_dir(1).name,...
            'continuous', 'Neuropix-PXI-100.ProbeA\');
    dat_files{part} = current_dat;
    end
    
    % Extract file paths and assign names for concatenation function
    hab_path = dat_files{1};
    p1_path = dat_files{2};
    p2_path = dat_files{3};

    % Check whether this is the renewal session (has an extra recording)
    if strcmp(session, 'Renewal')
        check_path = dat_files{4};
    else
        check_path = [];
    end
    
   % Assign filepaths for TTL files in each folder to a cell array
    TTL_files = cell(length(raw_folder_list), 1);

    for part = 1:length(raw_folder_list)
        % Handle inconsistent experiment and recording numbers (can be
        % variable)
        exp_dir = dir(fullfile(raw_folder_list{part}, 'Record Node 101', 'experiment*'));
        rec_dir = dir(fullfile(exp_dir(1).folder, exp_dir(1).name, 'recording*'));
        current_TTL = fullfile(rec_dir(1).folder, rec_dir(1).name,...
            'events', 'Neuropix-PXI-100.ProbeA', 'TTL\');
    TTL_files{part} = current_TTL;
    end

    hab_path_triggers = TTL_files{1};
    p1_path_triggers = TTL_files{2};
    p2_path_triggers = TTL_files{3};
    
    if strcmp(session, 'Renewal')
        check_path_triggers = TTL_files{4};
    else
        check_path_triggers = [];
    end    
   
    % Extract mouse/session basename for file save names
    base_name = regexp(dat_files{1}, ['mouse\d+_' session], 'match', 'ignorecase');
    
    % Dynamically construct save paths for input into concatenation and 
    % trigger extraction functions
    save_path = fullfile(mouse_path, mouse_name, session, 'Neural Data', 'Concatenated Data');
    save_name = [base_name{1} '_concatenated_neural_data.dat'];
        
    save_path_triggers = fullfile(mouse_path, mouse_name, session, 'Neural Data', 'Triggers');
    save_name_triggers = [base_name{1} '_'];
    
    % Extract triggers (in original neural sample rate terms)
    disp('Extracting trigger timings...')
    if exist(save_path_triggers, 'dir')
        if strcmp(session, 'Extinction')
            extract_neural_data_triggers(hab_path, p1_path, p2_path,...
        [], hab_path_triggers, p1_path_triggers, p2_path_triggers,...
        [], save_path_triggers, save_name_triggers)
        elseif strcmp(session, 'Renewal')
            extract_neural_data_triggers(hab_path, p1_path, p2_path,...
        check_path, hab_path_triggers, p1_path_triggers, p2_path_triggers,...
        check_path_triggers, save_path_triggers, save_name_triggers)
        end
    else
        warning('%s does not exist, skipping trigger extraction', save_path_triggers)
    end
    
    % Check whether data for this mouse/session has already been
    % concantenated, skip if so
    if exist(fullfile(save_path, save_name), 'file')
        warning('%s already exists, skipping concatenation', save_name)
        continue
    end
    
    % Check whether save folder for concatenated data exists, skip if not
    if ~exist(save_path, 'dir')
        warning('%s does not exist, skipping concatenation', save_path)
        continue
    end
    
    hab_path = [hab_path 'continuous.dat'];
    p1_path = [p1_path 'continuous.dat'];
    p2_path = [p2_path 'continuous.dat'];
    if strcmp(session, 'Renewal')
        check_path = [check_path 'continuous.dat'];
    end

    % Run concatenation function, checking for extra renewal recording
    disp(['Concatenating neural data for ', base_name{1}, '...'])
    if strcmp(session, 'Extinction')
        concatenate_neural_data_clean(hab_path, p1_path, p2_path, [], save_path, save_name)
        disp([save_name, ' has been saved!'])
    elseif strcmp(session, 'Renewal')
        concatenate_neural_data_clean(hab_path, p1_path, p2_path, check_path, save_path, save_name)
        disp([save_name, ' has been saved!'])
    end   

end

%% Extract event times for aligning data from trigger timing files

% % Constants
% fs = 30000;  % Neuronal data sampling rate (30kHz) 
% n_trials = 20; % Number of trials in each session part
% onset = 10; % Second of onset for stimuli
% pseudo_events = 20; % Number of events to extract from habituation
% 
% % Get master neural data folders for specified session
% for mouse = 2:length(mouse_files) % Always skip mouse 1 (data wrong format)
%     mouse_name = mouse_files(mouse).name;
%     mouse_path = mouse_files(mouse).folder;
%     current_trigger_folder = fullfile(mouse_path, mouse_name,...
%         session, 'Neural Data', 'Triggers');
% 
%     % Construct list of filepaths for trigger folders for current mouse
%     trigger_folder_dir = dir(fullfile(current_trigger_folder, 'mouse*'));
%     trigger_file_list = cell(length(trigger_folder_dir), 1);
% 
%     for part = 1:length(trigger_folder_dir)
%         trigger_file_list{part} = fullfile(trigger_folder_dir(part).folder,...
%             trigger_folder_dir(part).name);
%     end
% 
%     % Sort parts into correct order prior to extraction
%     if strcmp(session, 'Extinction')
%         idx = zeros(3,1);
%         key_words = {'habituation', 'p1', 'p2'};
%         for part = 1:numel(key_words)
%             idx(contains(trigger_file_list, key_words{part})) = part;
%         end
%     elseif strcmp(session, 'Renewal')
%         idx = zeros(4,1);
%         key_words = {'habituation', 'p1', 'p2', 'checkerboard'};
%         for part = 1:numel(key_words)
%             idx(contains(trigger_file_list, key_words{part})) = part;
%         end
%     end
% 
%     [~, sort_idx] = sort(idx);
%     trigger_file_list = trigger_file_list(sort_idx);
% 
%     % Check which stim set each mouse received
%     if ismember(mouse, received_stim_set_1)
%         which_stim_set = 1;
%     elseif ismember(mouse, received_stim_set_2)
%         which_stim_set = 2;
%     end
% 
%     evt_flash = [];
%     evt_loom = [];
%     evt_hab = [];
%     evt_checkerboard = [];
% 
%     save_name = [lower(mouse_name), '_', lower(session), '_extracted_events.mat'];
%     save_file = fullfile(current_trigger_folder, save_name);
%     if exist(save_file, 'file')
%         warning('%s already exists, skipping extraction', save_file)
%     else
%         disp('Extracting event timings...')
%         % Iterate over all session parts
%         for part = 1:length(trigger_file_list) 
% 
%             load(trigger_file_list{part}, 'evt'); 
%             % Convert sampling number into seconds (divide by sampling rate)
%             evt = double(evt)/fs; 
% 
%             % Remove first event (artefact)
%             evt = evt(2:end);     
% 
%             % Extract arbitrary, evenly-spaced events from habituation (no stimuli)
%             if part == 1
%                 evt_hab = sort_evts_hab_clean(evt, pseudo_events); 
% 
%             % Extract loom and flash events (i.e. triggers of stimulus onset)
%             elseif part == 2 || part == 3    
%                 [evt_loom1,evt_flash1] = sort_evts_clean(evt, n_trials, onset, part, ...
%                     which_stim_set);
% 
%                 evt_flash = cat(1,evt_flash,evt_flash1);
%                 evt_loom =  cat(1,evt_loom ,evt_loom1);
% 
%                 % Extract every trigger for checkerboard stim
%             elseif part == 4
%                 evt_checkerboard = evt;    
%             end 
% 
%         end   
% 
%         % Check session and save extracted events
%         if strcmp(session, 'Extinction')
%             save(save_file, 'evt_hab', 'evt_loom', 'evt_flash');
%             disp([save_name, ' has been saved!'])
%         elseif strcmp(session, 'Renewal')
%             save(save_file, 'evt_hab', 'evt_loom', 'evt_flash', 'evt_checkerboard');
%             disp([save_name, ' has been saved!'])
%         end     
% 
%     end
% 
% end
% 
% 
