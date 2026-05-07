load("Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)\Mouse 10\Extinction\Behavioural Data\SSM Fitted Data\mouse10_extinction_p1_SSM_fitted.mat") % load SSM file

frame_range = 1000:2000; % adjust to segment of interest

Tseg = T(:, frame_range);

dx = [0, diff(Tseg(1,:))];
dy = [0, diff(Tseg(2,:))];
vel = sqrt(dx.^2 + dy.^2);

bseg = b(:, frame_range);

missing_seg = missing(:, frame_range);
missing_sum = sum(missing_seg, 1); % number of missing points per frame

figure;

% --- Velocity ---
subplot(3,1,1)
plot(frame_range, vel, 'k')
ylabel('Velocity')
title('Translation (T) velocity')

% --- Shape (b coefficients) ---
subplot(3,1,2)
plot(frame_range, bseg')
ylabel('b coefficients')
legend('b1','b2','b3') % adjust if more

% --- Missingness ---
subplot(3,1,3)
plot(frame_range, missing_sum, 'r')
ylabel('Missing points')
xlabel('Frame')
title('Tracking quality')