% Paths
kilosort_path = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 4\Extinction\Neural Data\Concatenated Data\kilosort4';

chanMap      = readNPY(fullfile(kilosort_path,'channel_map.npy'));        % length = nCh
chanPos      = readNPY(fullfile(kilosort_path,'channel_positions.npy'));  % [nCh x 2], columns: x,y (µm)
spikeClu     = readNPY(fullfile(kilosort_path,'spike_clusters.npy'));     % [nSpikes x 1]
spikeTemp    = readNPY(fullfile(kilosort_path,'spike_templates.npy'));    % [nSpikes x 1]
templates    = readNPY(fullfile(kilosort_path,'templates.npy'));          % [nTemplates x nTime x nChan]
amplitudes   = readNPY(fullfile(kilosort_path,'amplitudes.npy'));         % [nSpikes x 1]

%%
% Peak-to-peak per channel per template
ptp = squeeze(max(templates,[],2) - min(templates,[],2));  % [nTemplates x nChan]

[ptpMax, bestChanIdx] = max(ptp,[],2);                     % best channel *index within recorded channels*
bestChan = chanMap(bestChanIdx);                           % map to original channel IDs if needed

% Template depth (y) using peak channel:
templY = chanPos(bestChanIdx,2);  % µm

% Assign each cluster a template (via its modal spike_template), then inherit depth:
cluIds = unique(spikeClu);
clu2templ = zeros(numel(cluIds),1);
cluDepth  = zeros(numel(cluIds),1);

for i = 1:numel(cluIds)
    c = cluIds(i);
    idx = (spikeClu == c);
    t = mode(spikeTemp(idx));     % dominant template for this cluster
    clu2templ(i) = t;
    cluDepth(i)  = templY(t+1);   % +1 if templates are 0-indexed in NPY
end

%%
T = table(cluIds, cluDepth, 'VariableNames', {'unit_id','depth_um'});
writetable(T, fullfile(kilosort_path,'unit_depths_um.csv'));
save(fullfile(kilosort_path,'unit_depths_um.mat'),'T');
