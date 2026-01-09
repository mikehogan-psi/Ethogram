function[good_trials] = extract_good_trials(X)
[Np,Nd,Ltrial,Ntrial] = size(X);
good_trials = false(1,Ntrial);
for t = 1:Ntrial
    temp = squeeze(X(1,1,:,t));
    if sum(isnan(temp))==0
        good_trials(t) = true;
    end
end