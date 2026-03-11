function[] = SSM_est_clean(master_directory, SSM_save_path)

%options
options_SSM.good_landm = [1:11];
options_SSM.ind_head = [1:7];
options_SSM.ind_body = [8:11];
options_SSM.K = 3;
options_SSM.TH_Eigen = 0.8;


% --- Progress bar setup ---
hwb = waitbar(0,'Starting SSM estimation...','Name','SSM progress');
cleanupObj = onCleanup(@() safeCloseWaitbar(hwb)); % ensures it closes on error

% Load data

all_files = dir(fullfile(master_directory, '**', '*triangulated.mat'));

X = [];
W = [];

for n = 1:numel(all_files)
    temp = load(fullfile(all_files(n).folder, all_files(n).name));
    X = cat(3, X, temp.X);
    W = cat(2, W, temp.W);
    
    % Waitbar loop
    if mod(n, max(1, round(numel(all_files)/200))) == 0 || n == numel(all_files)
        waitbar(0.02 + 0.18*(n/numel(all_files)), hwb, ...
            sprintf('Loading files: %d / %d', n, numel(all_files)));
    end
end


%remove bad head points
waitbar(0.22, hwb, 'Selecting landmarks...');
good_landm = options_SSM.good_landm;
X = X(good_landm,:,:);
W = W(good_landm,:);


% %take a random sample
% Nsamp = 4500;
% Nframe = size(X,3);
% indrand = randperm(Nframe);
% indrand = indrand(1:Nsamp);
% X = X(:,:,indrand);
% W = W(:,indrand);

% take sample (but ONLY frames that have NO missing datapoints!!)
min_visible_landmarks = 8;
Nsamp = 4500;
% Count number of visible landmarks per frame
waitbar(0.25, hwb, 'Counting visible landmarks per frame...');
Nframe = size(X, 3);
num_visible = zeros(1, Nframe);

for i = 1:Nframe
    % Landmark is visible if it has all 3 coordinates non-NaN
    visible_landmarks = sum(all(~isnan(X(:,:,i)), 2));
    num_visible(i) = visible_landmarks;
    % Waitbar loop
    if mod(i, max(1, round(Nframe/200))) == 0 || i == Nframe
        waitbar(0.25 + 0.20*(i/Nframe), hwb, ...
            sprintf('Scanning frames: %d / %d', i, Nframe));
    end
end

% Find frames with at least 8 visible landmarks
waitbar(0.45, hwb, 'Sampling good frames...');
good_frames = find(num_visible >= min_visible_landmarks);
% Randomly select up to Nsamp of them
Nsamp = min(Nsamp, numel(good_frames));  % make sure you don't oversample
indrand = randsample(good_frames, Nsamp);
% Subset X and W
X = X(:,:,indrand);
W = W(:,indrand);


%fit head
waitbar(0.52, hwb, 'Fitting head template...');
template_head = Estimate_mean_RANSAC(X(options_SSM.ind_head,:,:),false);
X(options_SSM.ind_head,:,:) = fit_template_v1(W(options_SSM.ind_head,:),X(options_SSM.ind_head,:,:),template_head);

%estimate mean
waitbar(0.60, hwb, 'Estimating mean pose (RANSAC) and aligning...');
template = Estimate_mean_RANSAC(X,false); % estimates mean_pose_3D (here termed template)
X = Alignment(X, template);

%recalculate mean
waitbar(0.68, hwb, 'Recalculating mean pose and realigning...');
template = nanmean(X,3);
X = Alignment(X, template);

%estimate eigenvectors
waitbar(0.75, hwb, 'Computing KNN for probPCA...');
X_KNN = Near_NaN_Euclidian(X, options_SSM.K, false);

waitbar(0.82, hwb, 'Running probPCA (this may take a while)...');
[~, Ndim, Cov_SSM, lambda, eignVectors, sigma2] = probPCA(X_KNN,options_SSM.TH_Eigen,true);
lambda = lambda(1:Ndim);
eignVectors = eignVectors(:,1:Ndim);

waitbar(0.92, hwb, 'Reshaping eigenvectors...');
for n = 1:Ndim
    eignV3D(:,:,n) = reshape(eignVectors(:,n),size(eignVectors,1)/3,3);
end

% Save
waitbar(0.97, hwb, 'Saving SSM...');
save(SSM_save_path,'lambda','template','eignV3D','sigma2', 'options_SSM', "Cov_SSM", 'X');
disp(['SSM saved as ' SSM_save_path])

waitbar(1.0, hwb, 'Done!');
pause(0.2);
safeCloseWaitbar(hwb);

%figure all poses
figure; hold on;
plot3(squeeze(X(options_SSM.ind_head,1,:)),squeeze(X(options_SSM.ind_head,2,:)),squeeze(X(options_SSM.ind_head,3,:)),'.b','MarkerSize',5)
plot3(squeeze(X(options_SSM.ind_body,1,:)),squeeze(X(options_SSM.ind_body,2,:)),squeeze(X(options_SSM.ind_body,3,:)),'k.','MarkerSize',5)
plot3(template(:,1),template(:,2),template(:,3),'+r','MarkerSize',30,'LineWidth',3);
xlim([-10 10]); ylim([-10 10]); zlim([-10 10]);

%show eigenvectors
figure; 
h = subplot(1,1,1); hold on;
pattern = 6*sin(2*pi*[1:100]/20);
for d = 1: Ndim 
    for t = 1:numel(pattern)
        subplot(h); hold on;
        Xfig = template+pattern(t)*sqrt(lambda(d))*eignV3D(:,:,d);
        plot3(Xfig(options_SSM.ind_head,1),Xfig(options_SSM.ind_head,2),Xfig(options_SSM.ind_head,3),'.b','MarkerSize',18);
        plot3(Xfig(options_SSM.ind_body,1),Xfig(options_SSM.ind_body,2),Xfig(options_SSM.ind_body,3),'.k','MarkerSize',18);
        xlim([-10 10]); ylim([-10 10]); zlim([-10 10]); title(d);
        view([60 40])
        %view([60 90])
        drawnow; pause(0.05); cla;
    end
end


function safeCloseWaitbar(hwb)
    if ~isempty(hwb) && isvalid(hwb)
        close(hwb);
    end
end

end