function[] = calculate_head_body_mike_abi(options)

options.Ltrial = 502;
options.wind = [0.2 0.6 0.2];
options.good_landm = [1:8];
options.ind_head = [1:4];
options.ind_body = [5:8];
options.alpha_reg = 0.1;
options.TH_missing = 0.5;
options.Lwind_missing = [1:3];

%
load('mouse3_extinction_p2_3D_triangulated.mat')
save_path = 'C:\Users\Abi Hogan\Documents\Psychedelics_Internship\behavior_analysis\pre_aquisistion_3D_analysis\results_head_body\';
%remove bad landmarks
X = X(options.good_landm,:,1:1506);
W = W(options.good_landm,1:1506);
%load SSM
load('SSM_3D_head_fitted_good_frames.mat','lambda','template','eignV3D','sigma2');
%fit SSM
disp(sprintf('Fit SSM'));
[b,X,R,T] = fit_SSM_Neuropix1_v2(X,W,template,lambda,eignV3D,options.alpha_reg,sigma2,options);
%load('so_far');
%separate into trials
disp(sprintf('Dividing data into trials'));
[Xt,Rt,Tt,bt] = divide_into_trials(X,R,T,b,options.Ltrial);
%extrapolate missing values
disp(sprintf('Extrapolating missing values'));
[Tt,Rt,bt] = extrapolate_missing_data(Tt,Rt,bt,options);
%regenerate X
Xt = generate_X(bt,Rt,Tt,template,eignV3D);
%extract good trials
good_trials = extract_good_trials(Xt);
%smooth time series
disp(sprintf('Smooth time series'));
[Rtf,Ttf,btf] = smooth_data(Rt,Tt,bt,good_trials,options.wind);
%regenerate X smoothed
Xtf = generate_X(btf,Rtf,Ttf,template,eignV3D);
%save
disp(sprintf('Saving results'));
save([save_path 'results_head_body.mat'],'Xtf','Ttf','Rtf','btf','template','eignV3D','lambda','good_trials');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%



function[X1,R1,T1,b1] = divide_into_trials(X,R,T,b,Ltrial)
Ntrial = floor(size(X,3)/Ltrial);
Np = size(X,1);
X1 = zeros(Np,3,Ltrial,Ntrial);
T1 = zeros(Ltrial,3,Ntrial);
R1 = zeros(3,3,Ltrial,Ntrial);
if size(b)>0
    Nb = size(b,2);
    b1 = zeros(Ltrial,Nb,Ntrial);
else 
    b1 = [];
end
for t = 1:Ntrial
    for p = 1:Np
        for d = 1:3
            X1(p,d,:,t) = X(p,d,Ltrial*(t-1)+1:Ltrial*t);
        end
    end
    for n = 1:3
        for m = 1:3
            R1(n,m,:,t) = R(n,m,Ltrial*(t-1)+1:Ltrial*t);
        end
    end
    for d = 1:3
        T1(:,d,t) = T(Ltrial*(t-1)+1:Ltrial*t,d);
    end
    if size(b)>0
        for d = 1:size(b,2)
            b1(:,d,t) = b(Ltrial*(t-1)+1:Ltrial*t,d);
        end
    end
end

function[T1,R1,b1] = extrapolate_missing_data(T1,R1,b1,options)
TH = options.TH_missing;
Lwind = options.Lwind_missing;
[Ltrial,Nd,Ntrial] = size(T1);
for t = 1:Ntrial
    for s = 1:3
        T1(:,s,t) = extrapolate_missing(squeeze(T1(:,s,t)),Lwind,TH);
    end
    if size(b1)
        for s = 1:size(b1,2)
            b1(:,s,t) = extrapolate_missing(squeeze(b1(:,s,t)),Lwind,TH);
        end
    end
    for r = 1:3
        for s = 1:3
            R1(r,s,:,t) = extrapolate_missing(squeeze(R1(r,s,:,t)),Lwind,TH);
        end
    end
    for l = 1:Ltrial
        [Ur,Sr,Vr] = svd(squeeze(R1(:,:,l,t)));
        R1(:,:,l,t) = Ur*Vr';
    end

end


