function[b,Xfit,R,T] = fit_SSM_Neuropix1_v2(X,W,mean_pose,lambda,P,alpha_reg,var_res,options)
Nshape = numel(lambda);
[Np,Nd,Nt] = size(X);
b0 = zeros(1,Nshape);
Xfit = NaN*zeros(Np,Nd,Nt);
R = NaN*zeros(3,3,Nt);
T = NaN*zeros(Nt,3);
b = NaN*zeros(Nt,Nshape);
C = NaN*zeros(1,Nt);
%Here generate all combinations
ind_body = options.ind_body;
ind_head = options.ind_head;
Nbody = numel(ind_body);
Nhead = numel(ind_head);
[head_num,body_num] = ndgrid(1:Nhead,1:Nbody);
comb = [head_num(:) body_num(:)];
comb = comb(sum(comb,2)>2,:);
comb = comb(sum(comb,2)<=(Nbody+Nhead-2),:);
Ncomb = size(comb,1);
%init outputs
Xfit1 = zeros(Np,Nd,Ncomb);
R1 = zeros(3,3,Ncomb);
T1 = zeros(Ncomb,3);
b1 = zeros(Ncomb,Nshape);
%start search
for t = 1:Nt
    [Wsort,ind_sort] = sort(W(:,t),'descend');
    ind_sort = ind_sort(Wsort>0);
    ind_sort_head = intersect(ind_sort,ind_head);
    ind_sort_body = intersect(ind_sort,ind_body);
    num_best_head = numel(ind_sort_head); 
    num_best_body = numel(ind_sort_body);
    if ((num_best_head+num_best_body)>2)
        gof = zeros(1,Ncomb);
        for c = 1:Ncomb
            if ((num_best_head>=comb(c,1))&(num_best_body>=comb(c,2)))
                ind1 = [ind_sort_head(1:comb(c,1)); ind_sort_body(1:comb(c,2))];
                btemp = fminunc(@(btemp) fit_SSM(btemp, X(ind1,:,t), mean_pose(ind1,:), lambda, P(ind1,:,:), alpha_reg, var_res), b0, optimoptions('fminunc','Display','none'));
                if all(abs(btemp)<10*std(lambda))
                    b1(c,:) = btemp;
                    [~, ~, R1(:,:,c), T1(c,:)] = fit_SSM(btemp, X(ind1,:,t), mean_pose(ind1,:), lambda, P(ind1,:,:), alpha_reg, var_res);
                    Xfit1(:,:,c) = mean_pose;
                    for m = 1:Nshape
                        Xfit1(:,:,c) = Xfit1(:,:,c)+b1(c,m)*P(:,:,m);
                    end
                    Xfit1(:,:,c) = Xfit1(:,:,c)*R1(:,:,c)+repmat(T1(c,:),Np,1);
                    temp1 = (squeeze(X(ind_sort,:,t)) - squeeze(Xfit1(ind_sort,:,c)));
                    dist = sum(temp1.^2,2);
                    gof(c) = 1/sum(dist.*W(ind_sort,t));
                end
            end
        end
        if max(gof)>0
            ind_best = find(gof == max(gof),1);
            Xfit(:,:,t) = Xfit1(:,:,ind_best);
            R(:,:,t) = R1(:,:,ind_best);
            T(t,:) = T1(ind_best,:);
            b(t,:) = b1(ind_best,:);
        end
    end
    if mod(t,100) == 0
        disp(sprintf('Pose %s of %s',num2str(t),num2str(Nt)));
    end
    graph = false;
    if graph
        figure; hold on;
        plot3(X(:,1,t),X(:,2,t),X(:,3,t),'.','MarkerSize',20);
        plot3(Xfit(:,1,t),Xfit(:,2,t),Xfit(:,3,t),'.','MarkerSize',20);
        xlim([-25 25]); ylim([-25 25]); zlim([-5 25]);
        ginput(); close all;
    end
end

function[C, Xfit, R, t] = fit_SSM(b, X, mu, lambda, P, alpha_reg, var_res)
Nshape = numel(lambda);
%
Xfit0 = mu;
for n = 1:Nshape
    Xfit0 = Xfit0+b(n)*squeeze(P(:,:,n));
end
%
[~, Xfit, tr] = procrustes(X,Xfit0,'Reflection',false, 'Scaling',false);
%
dist = sum((X(:)-Xfit(:)).^2)/var_res;
%
reg = alpha_reg*sum((b.^2)./lambda');
%
C = dist+reg;
%
R = tr.T;
t = tr.c(1,:);
