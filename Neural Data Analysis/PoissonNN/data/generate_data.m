function[] = generate_da0t0
a(which_test)
Nsamples_cell = 1500;
Ncell = 4;
N = Nsamples_cell*Ncell;
Ndim = 5;
%init randn data
X0 = randn(Nsamples_cell,Ndim);
%perform some nonlinear operations
A = randn(Ndim);
H = tanh(X0*A);
%init linear heads
W = randn(Ndim, Ncell);

%generate the data
cell_ids = [];
X = [];
y = [];
for n = 1:Ncell
    X = cat(1, X, X0);
    yc = poissrnd(exp(H*W(:,n)));
    y = cat(1, y, yc);
    cell_ids = cat(1, cell_ids, (n-1)*ones(Nsamples_cell,1));
end
y(y>200) = 200;
cell_ids = cell_ids + (which_test-1)*Ncell;
rec_id = which_test*ones(numel(cell_ids),1);

%transpose
X = X';
y = y';
cell_ids = cell_ids';
rec_id = rec_id';

figure; 
plot(y);

save(['test' num2str(which_test) '.mat'],'X','cell_ids','y','rec_id','-v7.3');
