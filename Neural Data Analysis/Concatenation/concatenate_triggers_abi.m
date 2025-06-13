%% Concatinating triggers
% This script allows the postprocessing and interpretation of datafiles
% containing event/trigger information
% 
%  1. read the NPY files conta generated 




%% setup filepaths

% define which mouse
   mouse = 'mouse2';

% define which session
    sesh = 'extinction';
  % sesh = 'renewal';

% Define folder path to specific session that you want the data to be concatenated from
common_dir = 'Z:\Abi\neuronal_data\mouse_2\Neural Data\Extinction'; % !!! CHANGE THIS !!!

% get folder paths to habituation, part1 and part2 (of Extinction or renewal session) 
    concistent_path = '\Record Node 101\experiment1\recording1\'; % this should be the same for every recording folder

    hab_path      = dir(fullfile(common_dir, '*habituation*'));
    hab_path      = [common_dir '\' hab_path(1).name concistent_path ];
    hab_path_TTL  = [hab_path 'events\Neuropix-PXI-100.ProbeA\TTL\'];
    hab_path_cont = [hab_path 'continuous\Neuropix-PXI-100.ProbeA\'];
    
    p1_path      = dir(fullfile(common_dir, '*p1*'));
    p1_path      = [common_dir '\' p1_path(1).name concistent_path ];
    p1_path_TTL  = [p1_path 'events\Neuropix-PXI-100.ProbeA\TTL\'];
    p1_path_cont = [p1_path 'continuous\Neuropix-PXI-100.ProbeA\'];
    
    p2_path      = dir(fullfile(common_dir, '*p2*'));
    p2_path      = [common_dir '\' p2_path(1).name concistent_path ];
    p2_path_TTL  = [p2_path 'events\Neuropix-PXI-100.ProbeA\TTL\'];
    p2_path_cont = [p2_path 'continuous\Neuropix-PXI-100.ProbeA\'];

% create output folder where data will be saved
filepath_out = [common_dir '\processed_triggers\'];
    mkdir(filepath_out);
    addpath(filepath_out);

% clear unnecessary variables
    clearvars common_dir concistent_path hab_path p1_path p2_path

%% Step 1: converting python files into MATLAB variables

% habituation data 
evt_hab    = readNPY([hab_path_TTL  'timestamps.npy']); 
states_hab = readNPY([hab_path_TTL  'states.npy']);
cont_hab   = readNPY([hab_path_cont 'timestamps.npy']);

% session part 1 data
evt_p1    = readNPY([p1_path_TTL  'timestamps.npy']); 
states_p1 = readNPY([p1_path_TTL  'states.npy']);
cont_p1   = readNPY([p1_path_cont 'timestamps.npy']);

% session part 2 data
evt_p2    = readNPY([p2_path_TTL  'timestamps.npy']); 
states_p2 = readNPY([p2_path_TTL  'states.npy']);
cont_p2   = readNPY([p2_path_cont 'timestamps.npy']);


%% concatinating triggers timeline (not alined to cont/sampling data timestamps)

% put all data into joined variabes (for loop)
    evt{1} = evt_hab; 
    evt{2} = evt_p1;
    evt{3} = evt_p2;
    states{1} = states_hab; 
    states{2} = states_p1;
    states{3} = states_p2;
    cont{1} = cont_hab; 
    cont{2} = cont_p1;
    cont{3} = cont_p2;

% initialize variables
t_start = 0;
evt_concat = [];

%conactinate evt data
for s = 1:3

evt{s} = evt{s}(states{s}==1);
evt{s} = evt{s}(2:end); 
evt{s} = evt{s} - evt{s}(1) + t_start;
t_start = evt{s}(end) + 30;
evt_concat = [evt_concat; evt{s}]; 

end

%% concatinating triggers from extrinction only (excluding habituation)

% put all data into joined variabes (for loop)
    evt{1} = evt_hab; 
    evt{2} = evt_p1;
    evt{3} = evt_p2;
    states{1} = states_hab; 
    states{2} = states_p1;
    states{3} = states_p2;
    cont{1} = cont_hab; 
    cont{2} = cont_p1;
    cont{3} = cont_p2;

% initialize variables
t_start = 0;
evt_extinction = [];

%conactinate evt data
for s = 2:3

evt{s} = evt{s}(states{s}==1);
evt{s} = evt{s}(2:end); 
evt{s} = evt{s} - evt{s}(1) + t_start;
t_start = evt{s}(end) + 30;
evt_extinction  = [evt_extinction ; evt{s}]; 

end


%% concatinating triggers timeline (alined to cont/sampling data timestamps)

% put all data into joined variabes (for loop)
    evt{1} = evt_hab; 
    evt{2} = evt_p1;
    evt{3} = evt_p2;
    states{1} = states_hab; 
    states{2} = states_p1;
    states{3} = states_p2;
    cont{1} = cont_hab; 
    cont{2} = cont_p1;
    cont{3} = cont_p2;

evt_concat_cont = [];

for s = 1:3

evt{s} = evt{s}(states{s}==1);
evt{s} = evt{s}(2:end); 
evt_concat_cont = [evt_concat_cont; evt{s}]; 

end


%% Step 2: Aligning trigger timestamps relative to recording start 

% put all data into joined variabes (for loop)
    evt{1} = evt_hab; 
    evt{2} = evt_p1;
    evt{3} = evt_p2;
    states{1} = states_hab; 
    states{2} = states_p1;
    states{3} = states_p2;
    cont{1} = cont_hab; 
    cont{2} = cont_p1;
    cont{3} = cont_p2;

evt_concat_rec_start = [];
cont0 = 0; % set start recording time to 0

for s = 1:3

evt{s} = evt{s}(states{s}==1);
evt{s} = evt{s}(2:end); 

evt{s} = evt{s}-cont{s}(1) +cont0;
cont0 = cont{s}(end);

evt_concat_rec_start = [evt_concat_rec_start; evt{s}]; 

end

%% Step 3: saving data
% 
save([filepath_out mouse '_' sesh 'concatinated_timestamps'],'evt_concat', 'evt_concat_cont', 'evt_concat_rec_start');
save([filepath_out mouse '_' sesh 'raw_data'],'evt_hab', 'evt_p1','evt_p2', 'cont_hab', 'cont_p1','cont_p2' , 'states_hab', 'states_p1', 'states_p2');


%% Plotting trigger timestamps

% Create a new figure
figure;
hold on;

% Set y-axis limits for the vertical lines
y_limits = [0 1];  % You can change this to fit your data context

% Loop through each event and plot a vertical line
for i = 1:length(evt_concat)
    x = evt_concat(i);
    line([x x], y_limits, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
end

% Label axes
xlabel('Time (s)');
ylabel('Trigger Marker');
title('Trigger Event Times');

% Improve appearance
ylim(y_limits);
xlim([min(evt_concat)-1, max(evt_concat)+1]);  % Add a bit of padding
grid on;

%%
% Convert evt_concat (in seconds) to minutes
evt_minutes = evt_concat / 60;

% Create a new figure
figure;
hold on;

% Split into first 4500 events (habituation) and the rest
n_events = length(evt_concat);
cutoff = min(4500, n_events);  % In case there are fewer than 4500 events

% Plot first 4500 triggers in red
xline(evt_minutes(1:cutoff), 'r');

% Plot remaining triggers (if any) in blue
if cutoff < n_events
    xline(evt_minutes(cutoff+1:end), 'b');
end

% Format plot
xlabel('Time (minutes)');
ylabel('Trigger Marker');
title('Trigger Event Times (Red: First 4500, Blue: After)');

% Set x-axis ticks every minute
xlim([0, ceil(max(evt_minutes)) + 1]);
xticks(0:1:ceil(max(evt_minutes)));
grid on;

% Optional aesthetic y-limits
ylim([0 1]);
