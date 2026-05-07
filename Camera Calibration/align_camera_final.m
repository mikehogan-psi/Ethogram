function[] = align_camera_final(video_path, cameraSet)
% Here load intrinsics
% find images that are shared across all cameras
% use one image as reference
% realign all cameras to this
% load 3D images (only need 2 or 3 cameras)
% triangulate 3D images
% compute alignment transformation between reference frames
% apply transformation to get the new R, T for all cameras
% recalculate P matrices for all cameras

% I should be able to calculate the transformation direcly right?

warning('off');

load([video_path '\camera_intrinsics.mat']);
cameraParams = cameraParams(cameraSet);
numCams = numel(cameraParams);

histValid = false(numCams, Nimgs);
for n = 1:numCams
    histValid(n,:) = hist(imageValidPerCam{n},[1:Nimgs]) == 1;
end

%pick one random img that is correctly detected by all cameras
indAllValid = find(sum(histValid)==numCams);
if numel(indAllValid) == 0
   error(['Could not find an image seen by all cameras. ' ...
       'Change the camera set of re-acquire the checkerboard video.']);
end
indValid = randsample(indAllValid,1);

%register all cameras to the same reference image
K = cell(1,numCams);
T = cell(1,numCams);
R = cell(1,numCams);
P = cell(1,numCams);
for n = 1:numCams
    indImg = find(imageValidPerCam{n} == indValid);
    K{n} = cameraParams{n}.Intrinsics.K;
    Rc = cameraParams{n}.PatternExtrinsics(indImg).R;
    Tc = cameraParams{n}.PatternExtrinsics(indImg).Translation'; 
    T_cam_world = [Rc, Tc; 0 0 0 1];
    T_world_cam = inv(T_cam_world);
    R{n} = Rc; %T_world_cam(1:3,1:3); 
    T{n} = Tc; %T_world_cam(1:3,4);
    P{n} = K{n}*[R{n} T{n}];
end

save([video_path '\Pcal.mat'],'P','K','R','T');
