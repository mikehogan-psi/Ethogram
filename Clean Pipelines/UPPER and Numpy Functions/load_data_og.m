function[x,L] = load_data_og(camera_files)

Nc = numel(camera_files);

for c = 1:Nc
data(:,:,c) = csvread(camera_files{c},3,1);
disp(sprintf('loaded camera %s',num2str(c)));
end

%extract x (homogeneous coord) and L
[Nt,Np,Nc] = size(data);
Np = Np/3;
x = ones(Nt,3*Np,Nc); %1 is
for p = 1:Np
for c = 1:Nc
x(:,3*(p-1)+1:3*(p-1)+2,c) = data(:,3*(p-1)+1:3*(p-1)+2,c);
end
end
L = data(:,3:3:end,:);