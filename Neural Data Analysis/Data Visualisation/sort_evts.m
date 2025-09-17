function[evt_loom,evt_flash] = sort_evts(evt,Ntrial,onset,which_epoch)

Ltrial = numel(evt)/Ntrial;
evt_mat = zeros(Ntrial,Ltrial);
for n = 1:Ntrial
    evt_mat(n,:) = evt((n-1)*Ltrial+1:n*Ltrial);
end
if which_epoch == 2
   stim = [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1];
elseif which_epoch == 3
   stim = [0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1];
end
evt_loom = evt_mat(stim==1,onset);
evt_flash = evt_mat(stim==0,onset);


