function[] = analyse_results_311025()
load('C:\Users\G71044MH\OneDrive - The University of Manchester\Documents\GitHub\Ethogram\Neural Data Analysis\PoissonNN\results\results.mat')
%
ind = find(pseudoR2>0.1);
pseudoR2 = pseudoR2(ind);
W = weights(ind,:);
[Ncell,Nw] = size(W);
%norm
for n = 1:Ncell
    W(n,:) = W(n,:)/norm(W(n,:));
end
%pca on weight
[coeff,Wpc,lambda] = pca(W);
Npc = min(find(cumsum(lambda)/sum(lambda) > 0.8));
Wpc = Wpc(:,1:Npc);
%figure 1
figure; hold on;
errorbar(mean(Wpc),std(Wpc)/sqrt(Ncell),'.')
xlim([0 Npc+1]);
xlabel('#PC');
ylabel('PC weights');
%figure 2
figure;
hist(pseudoR2,[-0.5:0.02:1])
hist(pseudoR2,[-1:0.02:1])
%figure 3
figure;
plot(cumsum(lambda)/sum(lambda),'.:')
%figure 4
figure;
plot(Wpc(:,1),Wpc(:,2),'.')