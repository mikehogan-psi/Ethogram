function[transform] = triangulate_new_final(video_path, camID)

numCams = numel(camID);
allPosesArray = []; 
allIntrinsics = [];
load([video_path '\Pcal.mat']);
P = P(camID);

for n = 1:numCams
    temp = load([video_path '\cam' num2str(camID(n)) '.mat']);
    CoordsCell{n} = temp.data.Coords(:,[1 2]);
    if n == 1
       Ncoord = size(temp.data.Coords,1);
       Xw = temp.X;
    end
end

X = zeros(Ncoord,3);
for n = 1:Ncoord
    x = [];
    for m = 1:numCams
        x = cat(2,x,[CoordsCell{m}(n,1); CoordsCell{m}(n,2); 1]);
    end
    indok = find(~isnan(x(1,:)));
    temp = ls_triangulate(x(:,indok),P(indok));
    X(n,:) = temp(1:3)'/temp(4);
end

%calculate Procrustes
% Perform Procrustes analysis
[~, ~, transform] = procrustes(Xw, X);
X = transform.b*X*transform.T + transform.c;

figure; hold on;
plot3(Xw(:,1),Xw(:,2),Xw(:,3),'k.','MarkerSize',12);
plot3(X(:,1),X(:,2),X(:,3),'r.','MarkerSize',12);

save('transform','transform');



