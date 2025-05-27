function[x1] = extrapolate_missing(x,Lwind,TH)
for n = 1:numel(Lwind)
    x1(n,:) = extrapolate_missing_sub(x,Lwind(n),TH);
end
x1 = nanmean(x1);


function[x1] = extrapolate_missing_sub(x,Lwind,TH)
L = numel(x);
ind_nan = find(isnan(x));
ind_nan = ind_nan((ind_nan>Lwind)&(ind_nan<L-Lwind));
ind_num = find(~isnan(x));
x1 = x;
for n = 1:numel(ind_nan)
    temp = x(ind_nan(n)-Lwind:ind_nan(n)+Lwind);
    ind_nan_temp = find(isnan(temp));
    ind_num_temp = find(~isnan(temp));
    if ((numel(ind_num_temp)/numel(temp)>TH)&(~isnan(temp(1)))&(~isnan(temp(end))))
        temp(ind_nan_temp) = pchip(ind_num_temp,temp(ind_num_temp),ind_nan_temp);
    end
    x1(ind_nan(n)-Lwind:ind_nan(n)+Lwind) = temp;
end