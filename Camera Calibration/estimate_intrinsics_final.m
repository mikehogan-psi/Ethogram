function[cameraParams] = estimate_intrinsics_final(video_path)

%video_path = 'C:\Users\mqbphrs3\OneDrive - The University of Manchester\code_cam_cal\calib280226';
warning('off');

%loading images
[imgs, worldPoints, boardSize] = load_video_set_params(video_path);
disp('loaded images')

Nimgs = numel(imgs{1});
Npoints = size(worldPoints,1);
numCams = numel(imgs);
cameraParams = cell(1,numCams);
imagePointsPerCam = cell(numCams,1);
for n = 1:numCams
    imagePointsPerCam{n} = zeros(Npoints,2,Nimgs);
end
imageSizePerCam = zeros(numCams,2);
imageValidPerCam = cell(numCams,1); %%%%%%%%%%%

makefig = true;
if makefig == true
   fig = figure; hold on;
end

for c = 1:numCams
    % Collect detections across the sequence for camera c
    imagePoints = {};
    imageSize = [];
    for k = 1:Nimgs %
        I = imgs{c}{k};
        [imagePoints_k, boardSizeFound] = detectCheckerboardPoints(I, ...
            'PartialDetections', false, 'MinCornerMetric', 0.15);
        
        if makefig == true    
            figure(fig); 
            imshow(I); hold on;
            try
                plot(imagePoints_k(:,1),imagePoints_k(:,2),'.r','MarkerSize',16);
            end
            pause(0.01); clf;
        end
        if ~isempty(imagePoints_k) && all(boardSizeFound == boardSize)
            imagePointsPerCam{c}(:,:,k) = imagePoints_k; 
            imageValidPerCam{c}(end+1) = k;
            if isempty(imageSize)
                imageSize = [size(I,2) size(I,1)];
            end
        end
        disp(sprintf('cam %s img %s',num2str(c), num2str(k)));
    end
    imageSizePerCam(c,:) = imageSize;
    
    % Estimate intrinsics for camera c
    cameraParams{c} = estimateCameraParameters(imagePointsPerCam{c}, worldPoints, ...
        'ImageSize', imageSizePerCam(c,:));

    %useful
    figure; showExtrinsics(cameraParams{c}); title(['Camera ' num2str(c)]);

end

save([video_path '\camera_intrinsics']);

%Alternatively I can also use this: 
%F = estimateFundamentalMatrix(c1, c2);
%Then perform self-calibration to estimate K
%Note: we assume the same camera took all images

