% convert settings.xml to kilosort channel map

[filename, path] = uigetfile(".xml", "Select the settings file");
data = readstruct(strcat(path, '\', filename));


name = char(data.SIGNALCHAIN(1).PROCESSOR(1).EDITOR.NP_PROBE.electrodeConfigurationPresetAttribute); % produces char rather than string
x_pos_str = data.SIGNALCHAIN(1).PROCESSOR(1).EDITOR.NP_PROBE.ELECTRODE_XPOS;
x_pos_chann = fieldnames(x_pos_str);
pattern = '\d+';
in= regexp(x_pos_chann,pattern,'match');
chanMap0ind = cellfun(@str2num,[in{:}])'; % chanMap0ind should start at 0
chanMap = chanMap0ind+1; %chanMap should start at 1
xcoords = cell2mat(struct2cell(x_pos_str)); % extracts correctly
connected = true([384,1]);
kcoords = double(ones([384,1])); % same as NP2
y_pos_str = data.SIGNALCHAIN(1).PROCESSOR(1).EDITOR.NP_PROBE.ELECTRODE_YPOS;
ycoords = cell2mat(struct2cell(y_pos_str)); % extracts correctly


save('ChanMap.mat','chanMap','xcoords','ycoords','chanMap0ind','connected','kcoords','name')

% plotting

pos_col = 0:15:9585;
n_shanks = 4;
y_pos = repmat(pos_col',[1,n_shanks*2]);
x_pos = [8 40 258 290 508 540 758 790];
x_pos = repmat(x_pos,[640,1]);

y_pos = reshape(y_pos,1,[]);
x_pos = reshape(x_pos,1,[]);

% plot
figure
scatter(x_pos,y_pos,'.k')
hold on
scatter(xcoords,ycoords,"ob")
title(name)