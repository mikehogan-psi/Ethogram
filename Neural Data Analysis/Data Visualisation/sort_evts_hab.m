function evt_hab = sort_evts_hab(evt, tpre, tpost)

Ltrial = ((tpre + tpost) * 15);
Ntrial = numel(evt)/ Ltrial ;

evt_mat = zeros(Ntrial,Ltrial);

for n = 1:Ntrial
    evt_mat(n,:) = evt((n-1)*Ltrial+1:n*Ltrial);
end

evt_hab = evt_mat(:,(tpre * 15));


