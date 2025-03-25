function[D_train,lik] = load_coordinates(filename)
Nfile = numel(filename);
Maindata = [];
for n = 1:Nfile
    Maindata = cat(2,Maindata,csvread(filename{n},3,1)');
end
[Nbp, Nframe] = size(Maindata);
Nbp = size(Maindata,1)/3;
lik = Maindata(3:3:end,:);
D_train = zeros(Nbp,2,Nframe);
for n = 1:Nframe
    D_train(:,1,n) = Maindata(1:3:end,n);
    D_train(:,2,n) = Maindata(2:3:end,n);
end
clear Maindata;