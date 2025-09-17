%% Extract individual frame from  video and use it to measure width of arena
% Load the video
obj = VideoReader("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort_4_06_05_25 (SC PAG Implanted Animals)\Extinction\Mouse 2\Video Data\camera5_mouse2_extinction_p1_2025-06-06-132222-0000.avi"); % Replace with your video file path

% Specify the frame number you want to save
frame_number = 100; % Replace with your desired frame number

% Read the frame
frame = read(obj, frame_number);

% Save the frame as a PNG
imwrite(frame, 'frame_100.png'); % Replace 'frame_100.png' with your desired filename

disp('Frame saved as PNG successfully!');

% Load an image
frame = imread('frame_100.png'); % Replace with your image file
imshow(frame);
title('Click two points to measure the width');

% Get user input for two points
[x, y] = ginput(2);

% Calculate pixel distance
arena_width_pixels = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);
disp(['Arena width in pixels: ', num2str(arena_width_pixels)]);
%% Plotting log of velocity data to find appropriate freezing velocity

% Transform trianguled data (cm) back into px
  % all_velocity_data = all_velocity_data / 0.048;



log_v_data = log10(all_velocity_data);
figure;
histogram(log_v_data, 60); 
xlabel('Log Velocity (pix/frame)');
xlim([-5 1]);
ylim([0 2200]);
ylabel('Frequency (Frames)');
title('Velocity Distribution - X data');
%print(gcf, 'velocity_distribution.png', '-dpng', '-r300');

%% Plotting (linear) velocity data
figure;
histogram(all_velocity_data); 
xlabel('Velocity (pix/frame)');
ylabel('Frequency (Frames)');
title('Velocity Distribution - X data');
print(gcf, 'velocity_distribution.png', '-dpng', '-r300');

%%
% Flatten and remove NaNs (if applicable)
vX = all_velocity_data(:);
vX = vX(~isnan(vX));

vXfit = all_velocity_data(:);
vXfit = vXfit(~isnan(vXfit));

mean_X = mean(vX);
median_X = median(vX);

mean_Xfit = mean(vXfit);
median_Xfit = median(vXfit);

fprintf('Original X Data:  Mean = %.4f cm/frame, Median = %.4f cm/frame\n', mean_X, median_X);
fprintf('Fitted Xfit Data: Mean = %.4f cm/frame, Median = %.4f cm/frame\n', mean_Xfit, median_Xfit);


%% Plotting two guassian curves over histogram of velocity data

% Assuming log_velocity_data is your log-transformed velocity data:
log_v_data = log10(all_velocity_data(:)); % Flatten to column vector
log_v_data = log_v_data(~isinf(log_v_data) & ~isnan(log_v_data)); % Remove invalid values

% Fit a GMM with 2 components (bimodal assumption)
GMModel = fitgmdist(log_v_data, 2);

% Extract the means and variances of the Gaussian components
mu1 = GMModel.mu(1); % Mean of first Gaussian
mu2 = GMModel.mu(2); % Mean of second Gaussian
sigma1 = sqrt(GMModel.Sigma(1)); % Std deviation of first Gaussian
sigma2 = sqrt(GMModel.Sigma(2)); % Std deviation of second Gaussian
pi1 = GMModel.ComponentProportion(1); % Weight of first Gaussian
pi2 = GMModel.ComponentProportion(2); % Weight of second Gaussian

% Define the Gaussian probability density functions
x = linspace(min(log_v_data), max(log_v_data), 1000);
g1 = pi1 * normpdf(x, mu1, sigma1); % First Gaussian
g2 = pi2 * normpdf(x, mu2, sigma2); % Second Gaussian

% Plot the histogram and GMM
figure;
histogram(log_v_data, 'Normalization', 'pdf'); % Normalized histogram
hold on;
plot(x, g1, 'r-', 'LineWidth', 2); % First Gaussian
plot(x, g2, 'b-', 'LineWidth', 2); % Second Gaussian
% plot(x, g1 + g2, 'k-', 'LineWidth', 2); % Combined GMM
legend('', ' Freezing Gaussian', 'Movement Guassian', 'GMM Fit');
xlabel('Log Velocity');
ylabel('Probability Density');
title('Gaussian Mixture Model for Log Velocity Data');
hold off;

% Find the intersection of the two Gaussians
syms v;
pdf1 = pi1 * (1 / (sqrt(2 * pi) * sigma1)) * exp(-0.5 * ((v - mu1) / sigma1)^2);
pdf2 = pi2 * (1 / (sqrt(2 * pi) * sigma2)) * exp(-0.5 * ((v - mu2) / sigma2)^2);

% Solve for the intersection point
intersection = vpasolve(pdf1 == pdf2, v, [min(x) max(x)]);
intersection_real_val = double(exp(intersection)); % Convert back to original scale


%% 

% Choose the number of standard deviations
n_std = 0.25; % Adjust this value as needed

% Identify the "freezing" Gaussian (lower mean)
if mu1 < mu2
    freezing_mean = mu1;
    freezing_std = sigma1;
else
    freezing_mean = mu2;
    freezing_std = sigma2;
end

% Calculate the adjusted threshold by taking away 1 stdv from gaussians'
% intersection point
adjusted_threshold_log = double(intersection) - n_std * freezing_std; % In log scale
adjusted_threshold = exp(adjusted_threshold_log); % Convert back to original scale
adjusted_threshold_cms = adjusted_threshold * 0.048 * 15;
disp(['Adjusted threshold velocity for freezing: ', num2str(adjusted_threshold), ' pix/frame ('...
    num2str(adjusted_threshold_cms), ' cm/s)']);

%%

freezing_mean + freezing_std
%%
% Plot the histogram and GMM
figure('Position', [100, 100, 1000, 600]); % Larger figure for poster

% Histogram with transparency (no legend entry)
h = histogram(log_v_data, 'Normalization', 'pdf', 'FaceColor', [0.7, 0.7, 0.7], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8); 
hold on;

% Plot the Gaussians
p1 = plot(x, g1, 'r-', 'LineWidth', 3); % Freezing Gaussian
p2 = plot(x, g2, 'b-', 'LineWidth', 3); % Movement Gaussian


% Labels and title
xlabel('Log Velocity (Pixels/Frame)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 18, 'FontWeight', 'bold');
% title('Gaussian Mixture Model for Log Velocity Data', 'FontSize', 20, 'FontWeight', 'bold');

% Improved legend (excluding the histogram)
legend([p1, p2], 'Freezing Gaussian', 'Movement Gaussian', ...
    'FontSize', 16, 'Location', 'northeast');

% Grid and visual tweaks
grid on;
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'off');
hold off;

% % Save figure
% print(gcf, 'velocity_distribution.png', '-dpng', '-r300');

%%

%all_velocity_data_real_units = all_velocity_data*0.048*15;
all_velocity_data_real_units = all_velocity_data;

% Assuming log_velocity_data is your log-transformed velocity data:
log_v_data_real_units = log10(all_velocity_data_real_units(:)); % Flatten to column vector
log_v_data_real_units = log_v_data_real_units(~isinf(log_v_data_real_units) & ~isnan(log_v_data_real_units)); % Remove invalid values

% Fit a GMM with 2 components (bimodal assumption)
GMModel = fitgmdist(log_v_data_real_units, 2);

% Extract the means and variances of the Gaussian components
mu1 = GMModel.mu(1); % Mean of first Gaussian
mu2 = GMModel.mu(2); % Mean of second Gaussian
sigma1 = sqrt(GMModel.Sigma(1)); % Std deviation of first Gaussian
sigma2 = sqrt(GMModel.Sigma(2)); % Std deviation of second Gaussian
pi1 = GMModel.ComponentProportion(1); % Weight of first Gaussian
pi2 = GMModel.ComponentProportion(2); % Weight of second Gaussian

% Define the Gaussian probability density functions
x = linspace(min(log_v_data_real_units), max(log_v_data_real_units), 1000);
g1 = pi1 * normpdf(x, mu1, sigma1); % First Gaussian
g2 = pi2 * normpdf(x, mu2, sigma2); % Second Gaussian

% Plot the histogram and GMM
figure('Position', [100, 100, 1000, 600]); % Larger figure for poster

% Histogram with transparency (no legend entry)
h = histogram(log_v_data_real_units, 'Normalization', 'pdf', 'FaceColor', [0.7, 0.7, 0.7], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8); 
hold on;

% Plot the Gaussians
p1 = plot(x, g1, 'r-', 'LineWidth', 3); % Freezing Gaussian
p2 = plot(x, g2, 'b-', 'LineWidth', 3); % Movement Gaussian


% Labels and title
xlabel('Log Velocity (cm/s)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 18, 'FontWeight', 'bold');
% title('Gaussian Mixture Model for Log Velocity Data', 'FontSize', 20, 'FontWeight', 'bold');

% Improved legend (excluding the histogram)
legend([p1, p2], '"Freezing" Gaussian', '"Movement" Gaussian', ...
    'FontSize', 16, 'Location', 'northeast');

% Grid and visual tweaks
grid on;
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'off');
hold off;


% Save figure
print(gcf, 'velocity_distribution.png', '-dpng', '-r300');