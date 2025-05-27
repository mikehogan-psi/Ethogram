function[] = SSM_est_mike_abi_implanted(triangulated_data_path, SSM_model_name_path) 

% this script had been adapted to the DLC model used to track implanted mice
% DLC_landmarks: 1. cable tip
%                2. back left courner (of implant)
%                3. back right courner (of implant)
%                4. nose
%                5. left ear
%                6. right ear
%                7. neck base
%                8. body anteriour
%                9. body posteriour
%                10. tail base
%                11. tail anteriour
%                12. tail posteriour
%                13. tail tip

%options
options_SSM.good_landm = [1:11];
options_SSM.ind_head = [1:7];      % including implant markers (1:3)
options_SSM.ind_body = [8:11];
options_SSM.K = 3;                 % number of eigenvectors to be included in SSM
options_SSM.TH_Eigen = 0.8;        % threshold for how much varianve eigenvectors need to explain

%load data
filepath = triangulated_data_path;
fname = dir(filepath);
fname = fname(3:end);
X = []; W = [];
for n = 1:numel(fname)
    temp = load([filepath '\' fname(n).name]);
    X = cat(3,X,temp.X);
    W = cat(2,W,temp.W);
end

%remove bad head points
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
Nframe = size(X, 3);
num_visible = zeros(1, Nframe);
for i = 1:Nframe
    % Landmark is visible if it has all 3 coordinates non-NaN
    visible_landmarks = sum(all(~isnan(X(:,:,i)), 2));
    num_visible(i) = visible_landmarks;
end
% Find frames with at least 8 visible landmarks
good_frames = find(num_visible >= min_visible_landmarks);
% Randomly select up to Nsamp of them
Nsamp = min(Nsamp, numel(good_frames));  % make sure you don't oversample
indrand = randsample(good_frames, Nsamp);
% Subset X and W
X = X(:,:,indrand);
W = W(:,indrand);


% %fit head 
% template_head = Estimate_mean_RANSAC(X(options_SSM.ind_head,:,:),false);
% X(options_SSM.ind_head,:,:) = fit_template_v1(W(options_SSM.ind_head,:),X(options_SSM.ind_head,:,:),template_head);

%estimate mean
template = Estimate_mean_RANSAC(X,false); % estimates mean_pose_3D (here termed template)
X = Alignment(X, template);

%recalculate mean
template = nanmean(X,3);
X = Alignment(X, template);

%estimate eigenvectors
X_KNN = Near_NaN_Euclidian(X, options_SSM.K, false);
[~, Ndim, Cov_SSM, lambda, eignVectors, sigma2] = probPCA(X_KNN,options_SSM.TH_Eigen,true);
lambda = lambda(1:Ndim);
eignVectors = eignVectors(:,1:Ndim);
for n = 1:Ndim
    eignV3D(:,:,n) = reshape(eignVectors(:,n),size(eignVectors,1)/3,3);
end

%save
save(SSM_model_name_path,'lambda','template','eignV3D','sigma2', 'options_SSM', "Cov_SSM", 'X');

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

