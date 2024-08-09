function Area = fun_Area(Nodes,Faces,Resolution)
Area = struct();
%% Identify nodes of interest for each side
nodeIndicesXMax = find(Nodes(:,1) <= max(Nodes(:,1)) & Nodes(:,1) >= max(Nodes(:,1)-Resolution)); %Area.Vertices.x_max = nodeIndicesXMax;
nodeIndicesXMin = find(Nodes(:,1) >= min(Nodes(:,1)) & Nodes(:,1) <= min(Nodes(:,1)+Resolution)); %Area.Vertices.x_min = nodeIndicesXMin;
nodeIndicesYMax = find(Nodes(:,2) <= max(Nodes(:,2)) & Nodes(:,2) >= max(Nodes(:,2)-Resolution)); %Area.Vertices.y_max = nodeIndicesYMax;
nodeIndicesYMin = find(Nodes(:,2) >= min(Nodes(:,2)) & Nodes(:,2) <= min(Nodes(:,2)+Resolution)); %Area.Vertices.y_min = nodeIndicesYMin;
nodeIndicesZMax = find(Nodes(:,3) <= max(Nodes(:,3)) & Nodes(:,3) >= max(Nodes(:,3)-Resolution)); %Area.Vertices.z_max = nodeIndicesZMax;
nodeIndicesZMin = find(Nodes(:,3) >= min(Nodes(:,3)) & Nodes(:,3) <= min(Nodes(:,3)+Resolution)); %Area.Vertices.z_min = nodeIndicesZMin;

%% PARAMETERS AND OUTPUTS
allNodeIndices = {nodeIndicesXMax, nodeIndicesXMin, nodeIndicesYMax, nodeIndicesYMin, nodeIndicesZMax, nodeIndicesZMin}; % Indices of the nodes of interest
facesForSides = struct('Faces', [], 'FaceVertices', []);% Initialize a structure array to store face information for each side
totalAreas = struct('Side', [], 'Area', []); % Store the areas per side 
sideColors = {'r', 'g', 'b', 'c', 'm', 'y'}; Area.side_colors = sideColors;% Colors for each side
Side_name = {'Right', 'Left', 'Top', 'Bottom', 'Front', 'Back'}; % Side names 

%% %%%%%%%%%%%%%%%     COMPUTE AREA OF EACH SIDE     %%%%%%%%%%%%%%%%%%%%%%
for sideIndex = 1:numel(allNodeIndices)
    sideFacesIndices = [];    % Initialize an empty array to store the indices of faces containing nodes on the side
    
    % Iterate over each face
    for i = 1:size(Faces, 1)
        faceIndices = Faces(i, :);    % Get the indices of the nodes for the current face
   
        % Check if all vertices are on the current side
        if all(ismember(faceIndices, allNodeIndices{sideIndex}))
            sideFacesIndices = [sideFacesIndices, i];
        end
    end
    
    sideFacesIndices = sideFacesIndices';
    facesForSides(sideIndex).Faces = Faces(sideFacesIndices, :);    % Store face information for the side
    sideVertices = struct('Node1', [], 'Node2', [], 'Node3', []);  

    for i = 1:length(sideFacesIndices)
        for j = 1:3
            vertexIndex = Faces(sideFacesIndices(i), j);
            sideVertices(i).(['Node' num2str(j)]) = Nodes(vertexIndex, :);
        end
    end
    facesForSides(sideIndex).FaceVertices = sideVertices;
    
    % Compute and display the total area for the current side
    totalArea = ComputeArea(sideVertices, length(sideFacesIndices));
    facesForSides(sideIndex).TotalArea = totalArea;

    vertex_coordinates = Nodes(allNodeIndices{sideIndex}, :);
    Area.Vertices(sideIndex).Name = Side_name{sideIndex};
    Area.Vertices(sideIndex).Location = allNodeIndices{sideIndex};
    Area.Vertices(sideIndex).Coordinates = vertex_coordinates;

    % Store total area in the structure
    totalAreas(sideIndex).Side = sideIndex;
    totalAreas(sideIndex).Area = totalArea;
    disp(['Total area of ' Side_name{sideIndex} ' side : ' num2str(totalArea)]);
    Area.Face(sideIndex).Side_Name = Side_name{sideIndex};
    Area.Face(sideIndex).Total_Areas = totalArea;
    Area.Face(sideIndex).Faces_indices = sideFacesIndices;
    Area.Face(sideIndex).Face_vertices = sideVertices;
    
end

%% %%%%%%%%%%%%%%%%%%%     Visualization code     %%%%%%%%%%%%%%%%%%%%%%%%%

% Plot each side with a different color
figure(); view(45, 10); axis equal ;
patch('Faces', Faces, 'Vertices', Nodes, 'facecolor', 'cyan')
hold on 
for sideIndex = 1:numel(facesForSides)
    sideColor = sideColors{sideIndex};

    % Plot faces on the side
    for i = 1:length(facesForSides(sideIndex).Faces)
        face = facesForSides(sideIndex).Faces(i, :);
        x = Nodes(face, 1);
        y = Nodes(face, 2);
        z = Nodes(face, 3);
        fill3(x, y, z, sideColor, 'EdgeColor', sideColor,'FaceAlpha',0.7);
        hold on;
    end

    % Plot nodes on the side
    scatter3(Nodes(allNodeIndices{sideIndex}, 1), ...
             Nodes(allNodeIndices{sideIndex}, 2), ...
             Nodes(allNodeIndices{sideIndex}, 3), ...
             50, sideColor, 'filled');
    hold on;
end

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
end

%% %%%%%%%%%%%%%%%%%%      SUB FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function totalArea = ComputeArea(face_vertices,numFaces)
totalArea = 0;
for i = 1:numFaces
    % Extract coordinates for each vertex of the triangle
    vertex1 = face_vertices(i).Node1;
    vertex2 = face_vertices(i).Node2;
    vertex3 = face_vertices(i).Node3;
    
    area = 0.5 * norm(cross(vertex2 - vertex1, vertex3 - vertex1)); % Compute the area of the triangle formed by the three vertices
    totalArea = totalArea + area;   % Add the area to the total
end

end