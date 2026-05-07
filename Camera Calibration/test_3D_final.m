function[] = test_3D_final(video_path, camID)

% Init data
global data;
data.camID = camID;
data.path = video_path;

% Load intrinsics7
load([video_path '\camera_intrinsics.mat'],'cameraParams');
data.cameraParams = cameraParams{camID};

% 3D world coordinates
studH = 0.17; pitchH = 0.96; pitchXY = 0.8; 
Xval = cumsum([0, 6*pitchXY, 6*pitchXY, 7*pitchXY, 6*pitchXY, 6*pitchXY]);
Yval = cumsum([0, 6*pitchXY, 6*pitchXY, 7*pitchXY, 6*pitchXY, 6*pitchXY]);
Zval = [studH, 5*pitchH+studH, 8*pitchH+studH];
X = []; Y = []; Z = [];
for n = 1:6
    X = [X; Xval'];
    Y = [Y; repmat(Yval(n),6,1)];
    if n == 1 | n == 4
        Z = [Z; [Zval(2) Zval(3) Zval(2) Zval(3) Zval(2) Zval(3)]'];
    elseif n == 3 | n == 6
        Z = [Z; [Zval(3) Zval(2) Zval(3) Zval(2) Zval(3) Zval(2)]'];
    else
        Z = [Z; Zval(1)*ones(6,1)];
    end
end
numPoints = numel(X);
data.X = [Y, X, Z];

% Init 2D coordinates
data.Coords = NaN*zeros(numPoints,2);

% Load images
f = dir([video_path '\*.jpg']);
data.img = imread([video_path '\' f(camID).name]);

% Init points count
data.pointsCount = 0;

% Generate the ui
fig = uifigure('Name', 'Image Cycle', 'Position', [100 100 1200 600]);
fig.HandleVisibility = 'on'; 

% Add the panels
tlo = tiledlayout(fig, 1, 2); 
axL = nexttile(tlo); % Left tile
axR = nexttile(tlo); 

% Add the buttons
labelBtn = uibutton(fig, 'push', 'Text', 'Label', ...
    'Position', [25 500 100 40], 'UserData', false, ...
    'ButtonPushedFcn', {@labelFcn, fig, tlo, axL, axR});
backBtn = uibutton(fig, 'push', 'Text', 'Back', ...
    'Position', [25 20 100 40], 'UserData', false, ...
    'ButtonPushedFcn', {@backFcn});
nextBtn = uibutton(fig, 'push', 'Text', 'Next', ...
    'Position', [1000 20 100 40], ...
    'UserData', false, ... 
    'ButtonPushedFcn', {@nextFcn}); 
saveBtn = uibutton(fig, 'push', 'Text', 'Save', ...
    'Position', [1000 500 100 40], ...
    'UserData', false, ... 
    'ButtonPushedFcn', {@saveFcn}); 


function[] = labelFcn(src, event, fig, tlo, axL, axR)
global data;
cla(axL); cla(axR); 
data.pointsCount = data.pointsCount+1;
stem3(data.X(:,1),data.X(:,2),data.X(:,3),'.k','MarkerSize',12,'Parent', axL);
hold(axL, 'on');
stem3(data.X(data.pointsCount,1),data.X(data.pointsCount,2), ...
    data.X(data.pointsCount,3),'.r','MarkerSize',20,'Parent', axL);
view(axL, -20,45);
xlabel('X','Parent',axL);
ylabel('Y','Parent',axL);
zlabel('Z','Parent',axL);
title('Proceed from 1 to 6','Parent', axL);
axes(axR);
imshow(data.img,'Parent', axR);
hold(axR, 'on');
data.Coords(data.pointsCount,:) = ginput(1); 
plot(data.Coords(:,1), data.Coords(:,2),'.r', 'MarkerSize', 12);

function[] = backFcn(src, event)
global data; 
data.pointsCount = data.pointsCount-1;

function[] = nextFcn(src, event)
global data; 
data.pointsCount = data.pointsCount+1;

function[] = saveFcn(src, event)
global data; 
X = data.X;
Coords = data.Coords;
camID = data.camID;
cameraParams = data.cameraParams;
indOK = find(~isnan(Coords(:,1)));
[worldPose, inlierIdx, status] = estworldpose(Coords(indOK,:),X(indOK,:),...
    cameraParams.Intrinsics,'MaxNumTrials',10000,'Confidence',0.99, ...
    'MaxReprojectionError',5);
cameraPose = invert(worldPose);
K = cameraParams.Intrinsics.K;
R = cameraPose.R;
T = cameraPose.Translation;
P = K*[R, T'];

%reprojec test
xtest = P*[X'; ones(1,size(X,1))];
xtest(1,:)  = xtest(1,:)./xtest(3,:);
xtest(2,:)  = xtest(2,:)./xtest(3,:);
figure; 
imshow(data.img)
hold on;
plot(xtest(1,:),xtest(2,:),'.r','MarkerSize',12);

%save
save([data.path '\cam' num2str(camID)],'X','data','camID','cameraPose',...
    'cameraParams','K','R','T','P');
