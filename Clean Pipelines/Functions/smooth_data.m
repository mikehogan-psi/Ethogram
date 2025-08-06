function[R1,T1,b1] = smooth_data(R1,T1,b1,good_trials,wind)
[Ltrial,Nd,Ntrial] = size(T1);
Nshape = size(b1,2);
for t = 1:Ntrial
    if good_trials(t)
        for n = 1:3
            for m = 1:3
                R1(n,m,:,t) = filtfilt(wind, 1, squeeze(R1(n,m,:,t)));
            end
        end
        for l = 1:Ltrial
            [Ur,Sr,Vr] = svd(squeeze(R1(:,:,l,t)));
            R1(:,:,l,t) = Ur*Vr';
        end
        for n = 1:3
            T1(:,n,t) = filtfilt(wind,1,squeeze(T1(:,n,t)));
        end
        for n = 1:Nshape
            b1(:,n,t) = filtfilt(wind,1,squeeze(b1(:,n,t)));
        end
    end
end