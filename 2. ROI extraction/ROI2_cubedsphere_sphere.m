function [Rotm,Faces,Vertices,Eigen_Values_ROI_2,Eigen_Vectors_ROI_2,ROI_2] = ROI2_cubedsphere_sphere(Full_Binary_sample,optimalCubeSize,Resolution)
[x,y,z] = ndgrid(-1:2:1) ;
ROI_Nodes = [x(:),y(:),z(:)]*optimalCubeSize/2 ;

[x_sample,y_sample,z_sample] = ind2sub(size(Full_Binary_sample),find(Full_Binary_sample)) ; % get the x,y,z matrices of the full binary sample 
Voxels_sample = [x_sample,y_sample,z_sample] ; % Store the location of the white pixels in a n*3 matrix of voxels

[Sphere,sphereRadius,Voxels,sphereCenter] = Sphere_Identfication(Full_Binary_sample,optimalCubeSize);
FV_ROI_1 = fun_FV(Sphere,0) ;
[Eigen_Values_ROI_1,Eigen_Vectors_ROI_1] = fun_TRI(FV_ROI_1) ;

%% RE-ORIENTATION OF THE FIRST EXTRACTED ROI TO ALIGN IT WITH ITS PRINCIPAL AXIS 
Rotm = [Eigen_Vectors_ROI_1(1:2,:) ; cross(Eigen_Vectors_ROI_1(1,:),Eigen_Vectors_ROI_1(2,:))] ; % Extract the Rotation Matrix of the ROI within the sample 
if ~(isequal(round(Rotm*Rotm',10),eye(3)) & isequal(round(Rotm'*Rotm,10),eye(3)) & round(det(Rotm),10)==1)
    display('Mauvaise matrice de rotation')
end  

[Axis,Angle] = fun_Rotm2Axang(Rotm) ; % Get the angle and the axis around which to rotate the ROI matrix 
Binary_Image_Rot = imrotate3(Full_Binary_sample,Angle*180/pi,Axis,'nearest') ; % Rotate the Full sample around the axis by the define angle 
[x,y,z] = ind2sub(size(Binary_Image_Rot),find(Binary_Image_Rot)) ; % get the x,y,z matrices of the Rotated full binary sample 
Voxels = [x,y,z] ; % Store the location of the white pixels in a n*3 matrix of voxels
Voxels_ROI_2 = Voxels(sum(abs((Voxels-mean(Voxels))*Resolution)<=optimalCubeSize/2,2)==3,:) ; % Extract the voxels of the oriented ROI from the rotated sample 
Voxels_ROI_2 = Voxels_ROI_2-min(Voxels_ROI_2)+1 ; % Centre the voxels 
ROI_2 = zeros([1,1,1]*max(Voxels_ROI_2(:))) ; % Create the matrix that will store the ROI voxels 
ROI_2(sub2ind(size(ROI_2),Voxels_ROI_2(:,1),Voxels_ROI_2(:,2),Voxels_ROI_2(:,3))) = 1 ; % Allocate  each voxels the value of 1 in the n*n*n matrix of ROI_2

FV_ROI_2 = fun_FV(ROI_2,0) ; % Extract the vertices coordinates and the connectivity of ROI_2 
[ROI_2] = ROI2_croping(ROI_2); % Crop the ROI_2 volume 
FV_ROI_2.vertices = (FV_ROI_2.vertices-mean(Voxels)) ; % Centre the ROI_2 matrix in the full sample 
Vertices = FV_ROI_2.vertices; % Store the re-oriented ROI vertices 
Faces = FV_ROI_2.faces; % Store the re-oriented ROI faces connectivity (triangles) 
[Eigen_Values_ROI_2,Eigen_Vectors_ROI_2] = fun_TRI(FV_ROI_2) ; % Compute the new Eigen Vectors of ROI_2 to verify its orientation 

% Calcul des angles avec le repère de référence
[~,idx] = max(abs(Eigen_Vectors_ROI_2),[],1) ;
Angles = diag(acos(Eigen_Vectors_ROI_2(idx,:)*eye(3))*180/pi) ;
Angles(Angles>90) = 180-Angles(Angles>90); 
Mean_angles = mean(Angles);

%% STL SAVING OF THE RE-ORIENTED ROI 
T = triangulation(Faces,Vertices); FileName = 'oriented_ROI2.stl';
stlwrite(T,FileName);

%% %%%%%%%%%%%%%%%%%%%     Visualization code     %%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the non-oriented ROI (ROI_1)
halfSize = floor(optimalCubeSize / 2); % Calculate half size, assuming cubeSize is odd

vertices1 = [-halfSize,-halfSize,-halfSize; -halfSize,-halfSize,halfSize; -halfSize,halfSize,-halfSize;
            -halfSize,halfSize,halfSize; halfSize,-halfSize,-halfSize; halfSize,-halfSize,halfSize;
             halfSize,halfSize,-halfSize; halfSize,halfSize,halfSize];
             vertices = vertices1 + sphereCenter;
faces = [ 1, 2, 6, 5; 2, 4, 8, 6; 3, 4, 8, 7; 1, 2, 4, 3;  1, 3, 7, 5; 5, 6, 8, 7];

% Plot the oriented ROI (ROI_2)
EigenVector = Eigen_Vectors_ROI_1;
orientedVertices = (EigenVector*vertices1'+sphereCenter')';
orientedFaces = [ 1, 2, 6, 5; 2, 4, 8, 6; 3, 4, 8, 7; 1, 2, 4, 3;  1, 3, 7, 5; 5, 6, 8, 7];

% Plot the Sphere 
[X, Y, Z] = sphere;
X = X * sphereRadius + sphereCenter(1);
Y = Y * sphereRadius + sphereCenter(2);
Z = Z * sphereRadius + sphereCenter(3);

figure; 
axis equal; hold on;  grid on; view (45,10);
patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'g', 'FaceAlpha', 0.5); % Patch of the vertices and faced of ROI_1
patch('Vertices', orientedVertices, 'Faces', orientedFaces, 'FaceColor', 'red', 'FaceAlpha', 0.8); % Patch of the vertices and faced of ROI_2 
surf(X, Y, Z, 'EdgeColor', 'none','FaceColor', 'c'); colormap('default');alpha(0.3); % Plot the surface of the sphere

% Plot the eigenvectors 
quiver3(sphereCenter(1),sphereCenter(2),sphereCenter(3),Eigen_Vectors_ROI_1(1,1),Eigen_Vectors_ROI_1(2,1),Eigen_Vectors_ROI_1(3,1), 70, 'LineWidth', 2,'Color','black'); 
quiver3(sphereCenter(1),sphereCenter(2),sphereCenter(3),Eigen_Vectors_ROI_1(1,2),Eigen_Vectors_ROI_1(2,2),Eigen_Vectors_ROI_1(3,2), 70, 'LineWidth', 2,'Color','black');
quiver3(sphereCenter(1),sphereCenter(2),sphereCenter(3),Eigen_Vectors_ROI_1(1,3),Eigen_Vectors_ROI_1(2,3),Eigen_Vectors_ROI_1(3,3), 70, 'LineWidth', 2,'Color','black');
scatter3(Voxels_sample(:,1),Voxels_sample(:,2),Voxels_sample(:,3),'o','MarkerEdgeAlpha',0.01,'MarkerEdgeColor','b'); % Full sample
xlabel('X'); ylabel('Y');  zlabel('Z'); 

%% %%%%%%%%%%%%%%%%%%%%%     SUB FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FV = fun_FV(Binary_Image,Smoothing)
[x,y,z] = ndgrid(1:size(Binary_Image,1),1:size(Binary_Image,2),1:size(Binary_Image,3)) ;
[Faces,Vertices] = MarchingCubes(x,y,z,Binary_Image,0.5) ;
if Smoothing
FV_temp = smoothSurfaceMesh(surfaceMesh(Vertices,Faces),1) ;
Faces = double(FV_temp.Faces) ; Vertices = double(FV_temp.Vertices) ;
end
FV.faces = Faces ; FV.vertices = Vertices ;
end

function [Eigen_Values,Eigen_Vectors] = fun_TRI(FV)
[Vcb,Fcb] = cubedsphere(10,'equiangular') ;
Normals = faceNormal(triangulation(double(FV.faces),double(FV.vertices))) ;
Surfaces = 0.5*sqrt(sum(cross(FV.vertices(FV.faces(:,2),:)-FV.vertices(FV.faces(:,1),:),FV.vertices(FV.faces(:,3),:)-FV.vertices(FV.faces(:,1),:),2).^2,2)) ;
Norms = sum(abs(Vcb*Normals').*repmat(Surfaces',size(Vcb,1),1),2) ;
[~,Eigen_Values,Eigen_Vectors] = ellipsoid_fit(Vcb.*Norms) ; 
Eigen_Vectors = Eigen_Vectors' ;
end
 
function [ROI_2] = ROI2_croping(ROI_2)
[nonzero_rows, nonzero_cols, nonzero_layers] = ind2sub(size(ROI_2), find(ROI_2)); %Find non-zero elements indices along each dimension
min_x = min(nonzero_cols); max_x = max(nonzero_cols); min_y = min(nonzero_rows); %Find the minimum and maximum indices along each dimension
max_y = max(nonzero_rows); min_z = min(nonzero_layers); max_z = max(nonzero_layers);
cube_size = max([max_x - min_x, max_y - min_y, max_z - min_z]) + 1; %Find the size of the smallest cube containing all non-zero elements
pad_x = (cube_size - (max_x - min_x + 1)) / 2; %Calculate the padding needed for each dimension
pad_y = (cube_size - (max_y - min_y + 1)) / 2;
pad_z = (cube_size - (max_z - min_z + 1)) / 2;
smallest_cube = zeros(cube_size, cube_size, cube_size); %Create the smallest cube containing all non-zero elements
smallest_cube(pad_y+1:pad_y+max_y-min_y+1, pad_x+1:pad_x+max_x-min_x+1, pad_z+1:pad_z+max_z-min_z+1) = ROI_2(min_y:max_y, min_x:max_x, min_z:max_z);
ROI_2 = smallest_cube;
end

    function [Axis,Angle] = fun_Rotm2Axang(R)
Trace = trace(R) ;
Angle = acos((Trace-1)/2) ;
if Angle > 0
    Axis = (1/(2*sin(Angle)))*[R(3,2)-R(2,3) ; R(1,3)-R(3,1) ; R(2,1)-R(1,2)]' ;
else
    Axis = [1,0,0] ;
end
Axis = Axis/norm(Axis) ;
    end

end
