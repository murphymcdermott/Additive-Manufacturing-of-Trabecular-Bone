%% ******************     GENERAL INFORMATIONS     ************************

% This function will go through the extraction of the CT-scan images and
% then use Otsu Method to do its binarization. For this stage, we perform
% two separate otsu filtering to ensure that we only keep the trabecular
% bone. Finally the images are binarized and stacked to build a 3D Matrix
% that will represent the full bone sample. 

% The second stage is to extract from this large volume an ROI that will be
% the final volume that we will analyse for the abaqus FEM tests. Some
% filtering is applied to obtain the best ratio between BVTV and ROI size.

% Once the ROI has been determined, we rotate it to align is eigen vectors
% with its principal axis to allow a simpler interpretation of the FEM
% resutls. 

%% ************************************************************************

%% %%%%%%%%%%%%%%%%%%%%%% INITIAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath("C:\Users\huppe\OneDrive\Documents\BioFabLab\Functions")); % path to function folder
CT_data_folder = 'C:\Users\huppe\OneDrive\Documents\BioFabLab\Images_Matrices\CT-data\175';
binary_folder = 'C:\Users\huppe\OneDrive\Documents\BioFabLab\Images_Matrices\Binary_Folder\175'; 
Matrix_Folder = 'C:\Users\huppe\OneDrive\Documents\BioFabLab\Images_Matrices\Full_Sample_Matrix';
ROI_Folder = 'C:\Users\huppe\OneDrive\Documents\BioFabLab\Images_Matrices\ROI_Matrix';
MatrixNb = 38;
Data = struct();
%% %%%%%%%%%%%%%%%%%%%%%%%% OTSU BINARIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imageList = dir(fullfile(CT_data_folder, '*.tif')); % get the list of all the images
for j = 1:length(imageList)
    imagePath = fullfile(CT_data_folder, imageList(j).name); % extract a single image
    rawImage = imread(imagePath); % read the selected image
    binary_Image = function_Otsu(rawImage); % binarize using Otsu method the selected image
    outputFileName = fullfile(binary_folder, ['binarized_', imageList(j).name]);
    imwrite(binary_Image, outputFileName); % Save the image in the desired folder
end
fprintf('Otsu thresholding and saving of all images completed.\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% STACKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binarizedList = dir(fullfile(binary_folder, '*.tif')); % get the list of binarized images
first_image = imread(fullfile(binary_folder, binarizedList(1).name)); 
[W, H] = size(first_image); % store the size of the initial image
D = numel(binarizedList); % number of layers for the stacking
binaryStack = zeros(W, H, D); % 3D stack of the sample 
binaryStack(:, :, 1) = first_image; % first layer of the stack 

for i = 2:D % stacking all the layers 
    binaryStack(:, :, i) = imread(fullfile(binary_folder, binarizedList(i).name)); % stacking the layers to build the 3D stack
end
Data.Binary_Stack = binaryStack;
Matrix_Name = sprintf('Binarized_sample_%d.mat', MatrixNb);
outputFilePath_Matrix = fullfile(Matrix_Folder, Matrix_Name);
save(outputFilePath_Matrix,'binaryStack','-mat');

binaryStack = logical(binaryStack); % Matrix to analyze
binaryStack_size = size(binaryStack); % Get the 3 dimensions sizes
Data.Binary_Stack_logical = binaryStack; Data.Binary_Stack_size = binaryStack_size;

%% %%%%%%%%%%%  Find the biggest diagonal of the cube  %%%%%%%%%%%%%%%%%%%%
nbLayers = binaryStack_size(3); Data.Layer_number = nbLayers;% Extract the number of z layers 
middleLayer = binaryStack(:, :, floor(nbLayers/2) + 1); Data.Middle_layer = middleLayer;% Adjusted to get the middle layer
[ROI_center_x,ROI_center_y,ROI_center_z,minDiam] = MaxCubeSize_CentreLoc(middleLayer,nbLayers);
Data.ROI_centre.ROI_centre_X = ROI_center_x; Data.ROI_centre.ROI_centre_Y = ROI_center_y; Data.ROI_centre.ROI_centre_Z = ROI_center_z;

%% %%%%%%%%%%%%%%%%%%%%%%M  Compute BVTV  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxCubeSize = floor(minDiam*cos(pi/4));
cubeSize = 1:1:floor(maxCubeSize); % Size of the cube
BVTV_values = zeros(length(cubeSize), 1);

for k = 1:length(cubeSize)
    % Define the ranges for each dimension
    xRange = (ROI_center_x - floor(cubeSize(k)/2)):(ROI_center_x + floor(cubeSize(k)/2)); xRange = max(1, min(xRange, binaryStack_size(1)));
    yRange = (ROI_center_y - floor(cubeSize(k)/2)):(ROI_center_y + floor(cubeSize(k)/2)); yRange = max(1, min(yRange, binaryStack_size(2)));
    zRange = (ROI_center_z - floor(cubeSize(k)/2)):(ROI_center_z + floor(cubeSize(k)/2)); zRange = max(1, min(zRange, binaryStack_size(3)));
    
    cube = binaryStack(xRange, yRange, zRange);    % Extract the cube from the original matrix
    BVTV_values(k) = (sum(cube(:)>0) / numel(cube))*100;    % BV/TV in %
end

%% %%%%%%%%%%%%%%%%%%  ROI Identification  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cube_BVTV = [(cubeSize)',BVTV_values]; Data.Cube_BVTV = Cube_BVTV;
last10Percent_size = round(0.1 * length(BVTV_values));
averageLast10Percent = mean(BVTV_values(end - last10Percent_size + 1:end));
variance = var(BVTV_values(end - last10Percent_size + 1:end));

% Define the range for filtering
lowerBound = averageLast10Percent - variance*averageLast10Percent;
upperBound = averageLast10Percent + variance*averageLast10Percent;

% Filter BVTV values within the specified range
filteredBVTV = BVTV_values(BVTV_values >= lowerBound & BVTV_values <= upperBound);
filteredCube_BVTV = (Cube_BVTV((end-length(filteredBVTV)+1:end),:));
meanFiltered = mean(filteredBVTV);

diffBVTV = [(filteredCube_BVTV),abs(filteredCube_BVTV(:,2)-meanFiltered)];
Threshold = 0.05*meanFiltered;

for l = 1:length(diffBVTV)

    if size(diffBVTV, 1) >= l && diffBVTV(l, 3) < Threshold
    disp(" SAMPLE OPTIMAL ROI PARAMETERS : ")
    disp("-------------------------------------")
    fprintf('Cube size :%.2f pixels\n', diffBVTV(l, 1));
    fprintf('BVTV :%.2f %%\n', diffBVTV(l, 2));
    disp("-------------------------------------")
    optimalCube_Size = diffBVTV(l, 1); 
    optimal_BVTV = diffBVTV(l, 2); 
    OptimalCube = binaryStack(floor(ROI_center_x - optimalCube_Size / 2):floor(ROI_center_x + optimalCube_Size / 2),floor(ROI_center_y - optimalCube_Size / 2):floor(ROI_center_y + optimalCube_Size / 2),floor(ROI_center_z - optimalCube_Size / 2):floor(ROI_center_z + optimalCube_Size / 2));
    end        
    break;
end
Data.Optimal_cube_size = optimalCube_Size; Data.Optimal_BVTV = optimal_BVTV; Data.ROI1 = OptimalCube;

%% %%%%%%%%%%%%%%%%%%%%%%%%%% ROI Saving  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROI_Name = sprintf('ROI1_%d.mat', MatrixNb);
outputFilePath_ROI = fullfile(ROI_Folder, ROI_Name);
save(outputFilePath_ROI,'OptimalCube','-mat');

%% %%%%%%%%%%%%%%%%%%%%%  VISULIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upperBound_array = ones(size(Cube_BVTV(:,2)))*upperBound;
lowerBound_array = ones(size(Cube_BVTV(:,2)))*lowerBound;

figure
scatter(Cube_BVTV(:,1), Cube_BVTV(:,2),'Color','b','LineWidth',1)
hold on 
scatter(filteredCube_BVTV(:,1), filteredCube_BVTV(:,2),'MarkerEdgeColor','g','MarkerFaceColor','g','LineWidth',1)
hold on 
scatter(Cube_BVTV(optimalCube_Size,1), Cube_BVTV(optimalCube_Size,2),'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1)
hold on
plot(Cube_BVTV(:,1),upperBound_array(:,1),'color','k','LineStyle','-')
hold on 
plot(Cube_BVTV(:,1),lowerBound_array(:,1),'Color','k','LineStyle','-')
legend('BVTV vs cube size', 'Within bound cube size', 'Optimal Cube size and BVTV', 'Upper Bound', 'Lower Bound'); % corrected legend
xlabel("Cube size [pixel]")
ylabel(" BVTV value [%]")
title("BVTV vs Cube size")

%% %%%%%%%%%%%%%%%%%%%%  ROI RE-ORIENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binarized_sample = double(Data.Binary_Stack);
[Rotm,Faces,Vertices,Eigen_Values_ROI_2,Eigen_Vectors_ROI_2,ROI_2] = ROI2_cubedsphere_sphere(binarized_sample,Data.Optimal_cube_size,1);
origin = min(Vertices);  Vertices = Vertices - origin+1;
Data.Faces = Faces; Data.Vertices = Vertices ; Data.Rotm = Rotm; Data.Eigen_values_ROI2 = Eigen_Values_ROI_2; Data.Eigen_vectors_ROI2 = Eigen_Vectors_ROI_2;
Data.ROI_2 = ROI_2; 

%% %%%%%%%%%%%%%%%%%%% SAMPLE CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.05; % threshold of acceptance 
v1 = Eigen_Values_ROI_2(1,1); v2 = Eigen_Values_ROI_2(2,1); v3 = Eigen_Values_ROI_2(3,1);% Extract norms for the current set of vectors
norms = sort([v1, v2, v3], 'descend');% Sort the norms to simplify comparison logic
DoA = norms(1)/norms(3); Data.DoA = DoA;% Calculate the degree of anisotropy (DoA)
fprintf("Degree of Anisotropy : %.2f \n" , DoA);

% Check the classification based on the tolerance
if isApproxEqual(v1, v2, alpha) && isApproxEqual(v2, v3, alpha)
disp("Sample is Isotropic") ; Data.Type = 'Isotropic';
elseif isApproxEqual(v1, v2, alpha) || isApproxEqual(v2, v3, alpha) || isApproxEqual(v1, v3, alpha)
disp("Sample is Transverse"); Data.Type = 'Transverse';
else
disp("Sample is Orthotropic"); Data.Type = 'Orthotropic';
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE AREA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Resolution = 1;
Data.Area = fun_Area(Data.Vertices,Data.Faces,Resolution);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% COMPUTE ROI STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ++++++++++  HEX MODEL

Resolution = 1; 
Young_Modulus_Bone = 600e6; Young_Modulus_Marrow = 10e0; Poisson_Ratio = 0.3; ROI = Data.ROI_2;
E_app_Struct_Hex = fun_Hex_StiffnessMat(Resolution, Young_Modulus_Bone, Young_Modulus_Marrow,Poisson_Ratio, ROI,1,6);

%% ++++++++++  TET MODEL

% Create the mesh 
Mesh_Tet = fun_FV_Tet(Data.ROI_2); % This function provide the triangular connectivity of both the Bone and the Marrow 

%Create the INP's --> Import the STL in ABAQUS and save the INP in the desired folder for the next steps 

% Store all the necessary info
FEM_Tet.test_Name = {'CompressionX','CompressionY','CompressionZ','ShearXY','ShearXZ','ShearYZ'};
FEM_Tet.ROI_inp = 'C:\\Users\\mh3407\\Documents\\BioFabLab\\Bone-38.inp'; 
FEM_Tet.Complement_inp = 'C:\\Users\\mh3407\\Documents\\BioFabLab\\Marrow-38.inp'; 
FEM_Tet.Sample_Name = 'Sample-38';
FEM_Tet.Material1 = 'Bone'; FEM_Tet.Matrial1_Youngs_Modulus = 10e9; FEM_Tet.Matrial1_Poisson_Coeff = 0.3; 
FEM_Tet.Material2 = 'Void'; FEM_Tet.Matrial2_Youngs_Modulus = 10e1; FEM_Tet.Matrial2_Poisson_Coeff = 0.3; 
FEM_Tet.Sets_name = {'Left','Right','Bottom','Tom','Front','Back'}; Data.Vertices
FEM_Tet.Set_range = [min(size(Data.ROI_2(:,1))), max(size(Data.ROI_2(:,1))); min(size(Data.ROI_2(:,2))), max(size(Data.ROI_2(:,2))); min(size(Data.ROI_2(:,3))), max(size(Data.ROI_2(:,3)))];
FEM_Tet.Resolution = 1;
FEM_Tet.Displacement = 1; 
FEM_Tet.Job_names = {'Job-CompressionX','Job-CompressionY','Job-CompressionZ','Job-ShearXY','Job-ShearXZ','Job-ShearYZ'};

E_app_Struct_Tet = fun_Tet_StiffnessMat(FEM_Tet);

%% %%%%%%%%%%%%%%%%%%%%%%%%% Clear the variables %%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Data E_app_Struct_Tet Mesh_Tet FEM_Tet

%% %%%%%%%%%%%%%%%%%%%%%  SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isEqual = isApproxEqual(v1, v2, alpha)
% Helper function to check if two values are approximately equal based on alpha
lowerBound = (1-alpha) * v1;
upperBound = (1+alpha) * v1;
isEqual = (v2 >= lowerBound) && (v2 <= upperBound);
end
