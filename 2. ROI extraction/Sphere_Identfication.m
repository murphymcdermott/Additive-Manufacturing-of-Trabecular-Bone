function [Sphere,sphereRadius,Voxels,sphereCenter] = Sphere_Identfication(Full_Binary_sample,optimalCubeSize)

%% %%%%%%%%%%%%%%%%%%     PARAMETERS & VARIABLES     %%%%%%%%%%%%%%%%%%%%%%
Full_Binary_sample = logical(Full_Binary_sample); % Convert to logical to analyze
Matrixsize = size(Full_Binary_sample); % Get the 3 dimensions sizes
nbLayers = Matrixsize(3); % Extract the number of z layers
middleLayer = Full_Binary_sample(:, :, floor(nbLayers/2) + 1); % Extract the middle layer
[xIndices, yIndices, zIndices] = ind2sub(size(Full_Binary_sample), find(Full_Binary_sample == 1)); %Identification of wite pixels in full sample


%% %%%%%%%%%%%%%%%%%%%     ROI_1 IDENTIFICATION     %%%%%%%%%%%%%%%%%%%%%%%
[ROI_center_x,ROI_center_y,ROI_center_z,~] = MaxCubeSize_CentreLoc(middleLayer,nbLayers); % Identify the centre of the Full Sample 
cubecenter = [ROI_center_x,ROI_center_y,ROI_center_z]; % Store the coordinates of the centre 
halfSize = floor(optimalCubeSize / 2); % Calculate half size, assuming cubeSize is odd

% Define the sphere's center and radius
sphereCenter = cubecenter; sphereRadius = sqrt(3) * halfSize;

% Preallocate ROI_half for efficiency, assuming you know the size, e.g., Matrixsize
Sphere = zeros(Matrixsize); % Assuming Matrixsize is [xSize, ySize, zSize]
Voxels = [xIndices, yIndices, zIndices];

% Iterate over each voxel to see if it's within the sphere
for j = 1:size(Voxels, 1)
    voxel = Voxels(j, :);
    distance = sqrt(sum((voxel - sphereCenter).^2)); % Euclidean distance from the sphere's center
    
    if distance <= sphereRadius
        x = Voxels(j, 1);
        y = Voxels(j, 2);
        z = Voxels(j, 3);
        Sphere(x, y, z) = 1; % Mark the voxel as inside the sphere
    end
end
sphereRadius1 = floor(sphereRadius);
Sphere = Sphere(floor(ROI_center_x-sphereRadius1):floor(ROI_center_x+sphereRadius1), ...
                                floor(ROI_center_y-sphereRadius1):floor(ROI_center_y+sphereRadius1), ...
                                floor(ROI_center_z-sphereRadius1):floor(ROI_center_z+sphereRadius1));
end