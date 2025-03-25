function[] = Estimate_SSM_final_abi(which_part, all_DLC_files, SSM_model_name)
%INPUT: which_part is 'body' or 'tail' depending whether you re going to
%estimate Statistical Shape Model for body or tail. They will be saved separately. 

%initial parameters
TH = 0.9; %likelihood threshold for outliers
Npose = 4500; %number of poses to estimate SSM, 4000-5000 is appropriate
K = 3; %K nearest neighbour

% intialise cell array to contain all filenames
filename = cell(1,length(all_DLC_files));

% loop through all filenames and load them into cell array
for i = 1:length(all_DLC_files)
    
    filename{i} = all_DLC_files(i).name;

end

%MIKE: Input body part
if strcmp(which_part,'body')
   ind = [1:7 11 12];
elseif strcmp(which_part,'tail')
   ind = [7:10]; 
else
   disp('Valid INPUTs: "body" or "tail"');
   return
end

%load coordinates
[D_train,lik] = load_coordinates(filename);
[Nbp,Nframe] = size(lik);

%remove some body parts (body or tail)
D_train = D_train(ind,:,:);
lik = lik(ind,:); 
Nbp = numel(ind);

%To generate the model identify poses with no outliers and select a random
%selction of Npose from them;
ind_good = find(min(lik)>TH);
Ngood = numel(ind_good);
ind_rand = randperm(Ngood);
Npose = min(Npose,Ngood);
ind_rand = ind_rand(1:Npose);
D_train = D_train(:,:,ind_rand);

%Estimate Mean Pose
mean_pose_2D = Estimate_mean_RANSAC(D_train, false);
mean_pose_2D(:,1) = mean_pose_2D(:,1)-mean(mean_pose_2D(:,1));
mean_pose_2D(:,2) = mean_pose_2D(:,2)-mean(mean_pose_2D(:,2));
%%%
Data_aligned=Alignment(D_train,mean_pose_2D);
%%%
Data_KNN = Near_NaN_Euclidian(Data_aligned, K, false);
%%%
[mean_pose_ppca, ~, Cov_pPCA, eignValues, eignVectors] = pPCA(Data_KNN,true);
mean_pose_2D_ppca = reshape(mean_pose_ppca,[2,Nbp,1])';

%save
save(SSM_model_name);

%final figure
figure; hold on;
for n = 1:Npose
    plot(Data_aligned(:,1,n),Data_aligned(:,2,n),'.','MarkerSize',8,'Color',0.666*[1 1 1]);
end
plot(mean_pose_2D_ppca(:,1),mean_pose_2D_ppca(:,2),'r.','MarkerSize',18);
xlim([-100 100]); ylim([-100 100]);

%animated video
T = 100;
figure; 
h = subplot(1,1,1); hold on;
for n = 1:3
    for m = 1:T
        pose = mean_pose_ppca + 3*sqrt(eignValues(n))*eignVectors(:,n)*sin(2*pi*m/(0.5*T));   
        pose_2D = reshape(pose,[2,Nbp,1])';
        plot(pose_2D(:,1),pose_2D(:,2),'r.','MarkerSize',18);
        xlim([-120 120]); ylim([-120 120]);
        title([which_part ' eigpose ' num2str(n)])
        drawnow; pause(0.05);
        cla;
    end
end