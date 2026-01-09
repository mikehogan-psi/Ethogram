function[Xfit,R,T,W] = fit_template_v1(W,X,mean_template)
[Np,Nd,Nt] = size(X);
Xfit = NaN*zeros(Np,Nd,Nt);
R = NaN*zeros(Nd,Nd,Nt);
T = NaN*zeros(Nt,Nd);
if numel(W) == 0
    W = ones(Np,Nt);
end
%combinations
Ncomb = 2^Np;
comb_str = dec2bin(0:Ncomb-1);
comb = false(size(comb_str,1),size(comb_str,2));
comb(comb_str=='1') = true;
comb = comb(find(sum(comb,2)>=3),:);
comb = comb(find(sum(comb,2)<=4),:);
Ncomb = size(comb,1);
Xfit1 = zeros(Np,3,Ncomb);
R1 = zeros(3,3,Ncomb);
T1 = zeros(Ncomb,3);
%start search
for t = 1:Nt
    [~,ind_sort] = sort(W(:,t),'descend');
    Nok = sum(W(:,t)>0);
    ind_sort = ind_sort(1:Nok);
    if Nok>=3
        dist = zeros(1,Ncomb);
        for n = 1:Ncomb
            ind1 = find(comb(n,:));
            if sum(isnan(X(ind1,1,t))) == 0
                [~,~,temp] = procrustes(squeeze(X(ind1,:,t)),mean_template(ind1,:),'Scaling',false,'Reflection',false);
                R1(:,:,n) = temp.T; T1(n,:) = temp.c(1,:);
                Xfit1(:,:,n) = mean_template*R1(:,:,n) + repmat(T1(n,:),Np,1);
                temp1 = (squeeze(X(ind_sort,:,t)) - Xfit1(ind_sort,:,n));
                dist_temp1 = sqrt(sum(temp1.^2,2));
                dist(n) = 1/sum(dist_temp1.*W(ind_sort,t));
            end
        end
        indmax = find(dist == nanmax(dist),1);
        if dist(indmax)>0
            R(:,:,t) = R1(:,:,indmax); 
            T(t,:) = T1(indmax,:);
            Xfit(:,:,t) = Xfit1(:,:,indmax); 
            W(:,t) = mean(W(:,t));
        end
    end 
    if mod(t,100) == 0
       disp(sprintf('Fit template to %s frames out of %s', num2str(t),num2str(Nt)));
    end
end
