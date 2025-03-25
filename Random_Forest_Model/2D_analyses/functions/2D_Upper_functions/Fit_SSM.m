function[Xfit,b,t] = Fit_SSM(filename, which_part, save_path, SSM_model_name)

% save variables from script
 filename = {filename}; 

%which_part is 'body' or 'tail' depending whether you re going to
%estimate Statistical Shape Model for body or tail. They will be saved separately. 

%init parameters
Ndim = 2; THeig = 0.9;
THlik = 0.8; %outliers

%load SSM
%MIKE: 
if strcmp(which_part,'body')
   ind = [1:7 11 12];
   load(SSM_model_name,'mean_pose_2D_ppca','eignValues','eignVectors');
   min_num = 6;
elseif strcmp(which_part,'tail')
   ind = [7:10]; 
   load(SSM_model_name,'mean_pose_2D_ppca','eignValues','eignVectors');
   min_num = 3;
else
   disp('Valid INPUTs: "body" or "tail"');
   return
end

%determine the number of eigenvalues
Neig = min(find(cumsum(eignValues)/sum(eignValues)>THeig));

%load coordinates
[X,lik] = load_coordinates(filename);
[Nbp,Nframe] = size(lik);

%remove some body parts (body or tail)
X = X(ind,:,:);
lik = lik(ind,:); 
Nbp = numel(ind);

%remove outliers
for n = 1:Nframe
    ind_out = find(lik(:,n)<THlik);
    if numel(ind_out)
        X(ind_out,:,n) = NaN;
    end
end

%rearrange model 
lambda = eignValues(1:Neig);
var_res = mean(eignValues(Neig+1:end));
mean_pose = mean_pose_2D_ppca;
eigen2=reshape(eignVectors(:,:),Ndim,Nbp,Nbp*Ndim);
for i=1:Nbp*Ndim
    P(:,:,i)=(eigen2(:,:,i)');
end
P = P(:,:,1:Neig);

%fit the 2D data
Xfit = zeros(Nbp,2,Nframe); 
b = zeros(Neig,Nframe);
C = zeros(1,Nframe);
A = zeros(1,Nframe);
T = zeros(2,Nframe);
missing = true(1,Nframe);
for n = 1:Nframe
    [Xfit(:,:,n),b(:,n),A(n),T(:,n),C(n),missing(n)] = fit_data(X(:,:,n),lambda,mean_pose,P,var_res,min_num);
    disp(sprintf('frame %s of %s',num2str(n),num2str(Nframe)));
end

%save data
save([save_path '\' filename{1}(1:end-4) '_' which_part '_fit'],'b','Xfit','C','X','missing','A','T'); 

function[Xfit,b,A,T,C,missing] = fit_data(X,lambda,mean_pose,P,var_res,min_num)
alpha_reg = 0.1;% 0.001;
Nshape = numel(lambda);
Nbp = size(X,1);
stop_search = false;
options = optimoptions('fminunc','Display','none');
%init shape parameters
b0 = zeros(Nshape,1);
%fit the SSM
ind_num = find(~isnan(X(:,1)));
missing = false;
if numel(ind_num)>=min_num
    b = fminunc(@(b) fit_SSM(b, X(ind_num,:), mean_pose(ind_num,:), lambda, P(ind_num,:,:), alpha_reg, var_res), b0, options);
    [C, ~, R, T] = fit_SSM(b, X(ind_num,:), mean_pose(ind_num,:), lambda, P(ind_num,:,:), alpha_reg, var_res);
    Xfit = mean_pose;
    for n = 1:Nshape
        Xfit = Xfit + b(n)*P(:,:,n);
    end
    Xfit = Xfit*R + repmat(T,Nbp,1);
    A = rotm2eul([R(1,:) 0; R(2,:) 0; 0 0 1]);
    A = A(1);
else
    missing = true;
    Xfit = NaN*mean_pose;
    T = [NaN; NaN];
    A = NaN;
    b = NaN*ones(Nshape,1);
    C = 1000;
end

%show results
make_fig = false;
if make_fig & ~missing & ind_num<6
    fig1 = figure; hold on
    plot(X(:,1),X(:,2),'k.','MarkerSize',16); 
    plot(Xfit(:,1),Xfit(:,2),'bo','MarkerSize',8,'LineWidth',2); 
    title(['mean fit: C = ' num2str(C)]);
    mX = nanmean(X);
    xlim([mX(1)-150 mX(1)+150]);
    ylim([mX(2)-150 mX(2)+150]);
    ginput(); close all;
end

function[C, Xfit, R, T] = fit_SSM(b, X, mu, lambda, P, alpha_reg, var_res)
Nshape = numel(lambda);
%
Xfit0 = mu;
for n = 1:Nshape
    Xfit0 = Xfit0+b(n)*P(:,:,n);
end
%
[~, Xfit, tr] = procrustes(X,Xfit0,'Reflection',false, 'Scaling',false);
%
dist = sum((X(:)-Xfit(:)).^2)/var_res;
%
reg = alpha_reg*sum((b.^2)./lambda);
%
C = dist+reg;
%
R = tr.T;
T = tr.c(1,:);
%%
