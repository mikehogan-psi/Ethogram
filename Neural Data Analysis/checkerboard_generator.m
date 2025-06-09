% Target image resolution (16:9)
img_width = 1920;
img_height = 1080;

% Desired number of checks (change to match spatial frequency)
n_checks_x = 24;
n_checks_y = 13;

% Integer check size (ensures clean tiling)
checker_size_x = floor(img_width / n_checks_x);
checker_size_y = floor(img_height / n_checks_y);

% Recalculate image size that will be created
final_width = checker_size_x * n_checks_x;
final_height = checker_size_y * n_checks_y;

% Generate base checkerboard pattern
[X, Y] = meshgrid(1:n_checks_x, 1:n_checks_y);
pattern = mod(X + Y, 2);

% Upsample to pixel resolution
checkerboard_normal = kron(pattern, ones(checker_size_y, checker_size_x));

% Pad to target size if needed
padded_checkerboard = zeros(img_height, img_width);
padded_checkerboard(1:final_height, 1:final_width) = checkerboard_normal(1:final_height, 1:final_width);

% Create inverted version
checkerboard_inverted = 1 - padded_checkerboard;

% Save images
imwrite(padded_checkerboard, 'checkerboard_normal.png');
imwrite(checkerboard_inverted, 'checkerboard_inverted.png');

% Display
figure;
subplot(1,2,1);
imshow(padded_checkerboard, []);
title('Normal Checkerboard');

subplot(1,2,2);
imshow(checkerboard_inverted, []);
title('Inverted Checkerboard');
