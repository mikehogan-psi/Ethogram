% use this to save velocity for a video after labelling a video with
% freezing using the RFM labelling and training GUI
% first caluclate velocity using the freezing extractor pipeline and select
% mouse you want/have labelled
all_velocity = all_velocity';

all_velocity = all_velocity(:);

p1 = all_velocity(1:10040);

extracted_features = p1;

save_path  = 'W:\Mike\Neuropixels_Fear_Conditioning\Models\RFMs\Freezing\Labelled Data';
save_name = 'mouse8_extinction_p1_freezing_labels';


save(fullfile(save_path, save_name), 'labels', 'extracted_features');
%%

freezing_idx = logical(labels);
non_idx = ~freezing_idx;

freezing_vel = extracted_features(freezing_idx);
non_vel = extracted_features(non_idx);

mean_non_vel = mean(non_vel);
mean_freeze_vel = mean(freezing_vel);

%%

% Choose common bin edges
edges = -4:0.1:4;

log_freeze = log10(freezing_vel);
log_non = log10(non_vel);

figure
h1 = histogram(log_freeze, edges, 'Normalization', 'probability');
hold on
h2 = histogram(log_non, edges, 'Normalization', 'probability');

% Set distinct colours and transparency
h1.FaceColor = [0 0.4470 0.7410]; % blue
h1.EdgeColor = 'none';
h1.FaceAlpha = 0.6;

h2.FaceColor = [0.8500 0.3250 0.0980]; % orange
h2.EdgeColor = 'none';
h2.FaceAlpha = 0.6;

hold off
xlabel('Value')
ylabel('Probability')
legend('Freeze Velocities','Non-Freeze Velocities')
