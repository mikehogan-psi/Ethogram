function[] = SSM_est_tail_only(triangulated_data_path, SSM_model_name_path)

disp('Starting tail-only SSM estimation...');

% options
options_SSM.good_landm = 10:13;
options_SSM.ind_tail = 1:4;
options_SSM.K = 3;
options_SSM.TH_Eigen = 0.8;

% load data
disp('Loading triangulated .mat files...');
filepath = triangulated_data_path;
fname = dir(filepath);
fname = fname(3:end);  % skip '.' and '..'
X = []; W = [];
for n = 1:numel(fname)
    disp(['  Loading file ' num2str(n) '/' num2str(numel(fname)) ': ' fname(n).name]);
    temp = load(fullfile(filepath, fname(n).name));
    X = cat(3, X, temp.X);
    W = cat(2, W, temp.W);
end

disp('Selecting tail landmarks...');
X = X(options_SSM.good_landm, :, :);
W = W(options_SSM.good_landm, :);

% filter good frames
disp('Filtering good frames...');
min_visible_landmarks = 2;
Nframe = size(X, 3);
num_visible = zeros(1, Nframe);
for i = 1:Nframe
    num_visible(i) = sum(all(~isnan(X(:,:,i)), 2));
end
good_frames = find(num_visible >= min_visible_landmarks);
Nsamp = min(4500, numel(good_frames));
indrand = randsample(good_frames, Nsamp);
X = X(:, :, indrand);
W = W(:, indrand);
disp(['  Retained ' num2str(size(X,3)) ' valid frames for model building.']);

% estimate mean and align
disp('Estimating template and aligning...');
template = Estimate_mean_RANSAC(X, false);
X = Alignment(X, template);
template = nanmean(X, 3);
X = Alignment(X, template);

% compute eigenvectors
disp('Running PPCA to compute shape model...');
X_KNN = Near_NaN_Euclidian(X, options_SSM.K, false);
[~, Ndim, Cov_SSM, lambda, eignVectors, sigma2] = probPCA(X_KNN, options_SSM.TH_Eigen, true);
lambda = lambda(1:Ndim);
eignVectors = eignVectors(:,1:Ndim);
for n = 1:Ndim
    eignV3D(:,:,n) = reshape(eignVectors(:,n), size(eignVectors,1)/3, 3);
end

% save
disp(['Saving model to: ' SSM_model_name_path]);
save(SSM_model_name_path, 'lambda', 'template', 'eignV3D', 'sigma2', 'options_SSM', 'Cov_SSM', 'X');

% visualise poses
disp('Plotting aligned tail poses...');
figure; hold on;
plot3(squeeze(X(:,1,:)), squeeze(X(:,2,:)), squeeze(X(:,3,:)), '.k', 'MarkerSize', 5);
plot3(template(:,1), template(:,2), template(:,3), '+r', 'MarkerSize', 30, 'LineWidth', 3);
xlim([-10 10]); ylim([-10 10]); zlim([-10 10]);

% visualise eigenvectors
disp('Animating deformation along eigenvectors...');
figure; 
h = subplot(1,1,1); hold on;
pattern = 6 * sin(2 * pi * (1:100) / 20);
for d = 1:Ndim 
    for t = 1:numel(pattern)
        subplot(h); hold on;
        Xfig = template + pattern(t) * sqrt(lambda(d)) * eignV3D(:,:,d);
        plot3(Xfig(:,1), Xfig(:,2), Xfig(:,3), '.k', 'MarkerSize', 18);
        xlim([-10 10]); ylim([-10 10]); zlim([-10 10]); title(['Mode ' num2str(d)]);
        view([60 40]);
        drawnow; pause(0.05); cla;
    end
end

disp('SSM estimation complete.');

end
