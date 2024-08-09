function [ROI_center_x,ROI_center_y,ROI_center_z,min_radius] = MaxCubeSize_CentreLoc(middleLayer,nbLayers)

[y, x] = find(middleLayer); % Extract the white pixels from the middle layer
distances = pdist2([x, y], [x, y]); % Compute the distance between all the white pixels 

%% Find the maximum distance and its indices
[maxDistance, maxIndices] = max(distances(:)); % Find the maximum distance between 2 pixels and the pixels indices
[row, col] = ind2sub(size(distances), maxIndices); % Extract the row and column of the maxDistance

%% Calculate the center of the maximum diameter
ROI_center_x = floor((x(row) + x(col)) / 2); % Compute the x centre of the binary image 
ROI_center_y = floor((y(row) + y(col)) / 2); % Compute the y centre of the binary image 
ROI_center_z = floor(nbLayers/2); % z centre already knowns 
fprintf('Center of maximum diameter: (%.2f, %.2f)\n', ROI_center_x, ROI_center_y);
fprintf('Maximum diameter: %.2f pixels\n', maxDistance);

%% Find the maximum and minimum radius from the centre 
distances_from_center = pdist2([x, y], [ROI_center_x, ROI_center_y]);% Calculate distances of all pixels from the center
pixel_data_matrix_3d = cat(3, x, y, distances_from_center);% Create a 3D matrix to store x, y, and distances_from_center

max_radius = max(distances_from_center);
fprintf('Maximum rayon: %.2f pixels\n', max_radius);

max_distance_per_column = zeros(length(unique(x)), 1);% Create a matrix to store only one maximum distance per column

% Iterate over unique x coordinates and find the maximum distance
for i = 1:length(unique(x))
    current_x = unique(x(i));
    indices = find(x == current_x);
    max_distance_per_column(i) = max(distances_from_center(indices));
end

min_radius = min(max_distance_per_column);
fprintf('Minimum rayon: %.2f pixels\n', min_radius);
min_distance_pixel = pixel_data_matrix_3d(pixel_data_matrix_3d(:,:,3) == min_radius, :);% Extract pixel with minimum distance
max_distance_pixels = pixel_data_matrix_3d(pixel_data_matrix_3d(:,:,3) == max_radius, :);% Extract pixels with maximum distance

%% %%%%%%%%%%%%%%%%%%%     Visualization code     %%%%%%%%%%%%%%%%%%%%%%%%%
figure; imshow(middleLayer); hold on;
plot(ROI_center_x, ROI_center_y, 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot([x(row), x(col)], [y(row), y(col)], 'r-', 'LineWidth', 2);
plot(max_distance_pixels(:, 1), max_distance_pixels(:, 2), 'bo', 'MarkerSize', 5, 'LineWidth', 2);
plot(min_distance_pixel(1, 1), min_distance_pixel(1, 2), 'mo', 'MarkerSize', 5, 'LineWidth', 2);
plot([ROI_center_x, min_distance_pixel(1, 1)], [ROI_center_y, min_distance_pixel(1, 2)], 'm-', 'LineWidth', 2);
plot([ROI_center_x, max_distance_pixels(:, 1)], [ROI_center_y, max_distance_pixels(:, 2)], 'b-', 'LineWidth', 2);
title('Binary Image with All White Pixels, Center, Max Distance Pixels, Min Distance Pixel, and Lines');
hold off

end