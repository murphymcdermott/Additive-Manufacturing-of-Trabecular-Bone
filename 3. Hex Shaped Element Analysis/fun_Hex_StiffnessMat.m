function [Mesh] = fun_Hex_StiffnessMat(Resolution, Young_Modulus_Bone, Young_Modulus_Marrow,Poisson_Ratio, ROI,Deformation,nbTests)


%% %%%%%%%%%%%%%%%%%%%%%%%     VARIABLE SET UP     %%%%%%%%%%%%%%%%%%%%%%%%
Young_Modulus_List = [Young_Modulus_Bone,Young_Modulus_Marrow] ; % Young's modulus of both materials 
[x,y,z] = ind2sub(size(ROI),find(ROI)) ; % Create x,y,z matrices that store the location of ROI white pixels 
Voxels = [x,y,z] ; Mesh.Voxels = Voxels; % Extract from x,y,z the coordiantes of the voxels and store them in a n*3 matrix of coordinates 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mesh Creation 
Size = max(Voxels(:))+1 ; % Mesh size 
Nodes = [repmat(1:Size,1,Size^2) ; reshape(repmat(repmat(1:Size,1,Size),Size,1),1,[]) ; reshape(repmat(1:Size,Size^2,1),1,[])]' ;
Nodes = (Nodes-mean(Nodes))*Resolution ; Mesh.Nodes = Nodes;
Connectivity = [1,2,2+Size,1+Size,1+Size^2,2+Size^2,2+Size+Size^2,1+Size+Size^2] ; % Set up the hex element shaped connectivity 
Elements = repmat(repmat(repmat(Connectivity,Size-1,1)+(0:Size-2)'*1,Size-1,1)+reshape(repmat((0:Size-2),Size-1,1),1,[])'*Size,Size-1,1)+reshape(repmat((0:Size-2),(Size-1)^2,1),1,[])'*Size^2 ; % apply the connectivity to the elements 
[~,Order] = sort([3,7,8,4,2,6,5,1]) ; % Sorts the nodes in the proper order within the elements 
Elements = Elements(:,Order) ; Mesh.Elements = Elements; 

% Node Sets creation 
idx_Nodes = unique(Elements) ; % Extract the unique nodes contained in the element matrix  
Set_1 = idx_Nodes(Nodes(idx_Nodes,1)==min(Nodes(idx_Nodes,1))) ; % Left
Set_2 = idx_Nodes(Nodes(idx_Nodes,1)==max(Nodes(idx_Nodes,1))) ; % Right 
Set_3 = idx_Nodes(Nodes(idx_Nodes,2)==min(Nodes(idx_Nodes,2))) ; % Bottom 
Set_4 = idx_Nodes(Nodes(idx_Nodes,2)==max(Nodes(idx_Nodes,2))) ; % Top 
Set_5 = idx_Nodes(Nodes(idx_Nodes,3)==min(Nodes(idx_Nodes,3))) ; % Front 
Set_6 = idx_Nodes(Nodes(idx_Nodes,3)==max(Nodes(idx_Nodes,3))) ; % Back 
Nodes_Sets = {Set_1,Set_2,Set_3,Set_4,Set_5,Set_6} ;
Mesh.Node_Sets = Nodes_Sets;

% Create element sets 
Element_Set_1 = find(logical(ROI(:))) ; % Extract the nodes indices contained in the ROI (Bone)
Element_Set_2 = setdiff(1:size(Elements,1),find(logical(ROI(:))))' ; % Extract the nodes indices contained in the Complementary volume (Bone marrow)
Elements_Sets = {Element_Set_1,Element_Set_2} ; % Store the element sets in a single cell structure 
Mesh.Element_Sets = Elements_Sets;
 
%% %%%%%%%%%%%%%%%%%%%     Visualization code     %%%%%%%%%%%%%%%%%%%%%%%%%

figure(1) ; set(gcf,'color','w') ; axis equal ; hold on ; view(25,10)
Faces = [1,2,3,4 ; 5,6,7,8 ; 1,2,6,5 ; 3,4,8,7 ; 2,3,7,6 ; 4,1,5,8] ;
for i_Element = 1:size(Element_Set_1,1)
    patch('faces',Faces,'vertices',Nodes(Elements(Element_Set_1(i_Element),:),:),'facecolor','cyan','facealpha',1)
end
for i_Element = 1:size(Element_Set_2,1)
    patch('faces',Faces,'vertices',Nodes(Elements(Element_Set_2(i_Element),:),:),'facecolor','y','facealpha',1)
end


%% %%%%%%%%%%%%%%%%%%%%%     Boundary Conditions     %%%%%%%%%%%%%%%%%%%%%%

Displacement = Deformation*(max(Nodes(unique(Elements),3))-min(Nodes(unique(Elements),3))) ;

% Format 1 --> Fixed face ; Format 2 --> moving face 
Format_1 = ['** Name: BC-%1.f Step: Initial Type: Displacement/Rotation',char(10),'*Boundary',char(10),'SET-%1.f, %1.f, %1.f',char(10)] ;
Format_2 = ['** Name: BC-%1.f Step: Step-1 Type: Displacement/Rotation',char(10),'*Boundary',char(10),'SET-%1.f, %1.f, %1.f, %d',char(10)] ;

% Compression xx
BC_1 = sprintf(Format_1,1,1,1,1) ; % Fix the left side  along x
BC_2 = sprintf(Format_2,2,2,1,1,-Displacement) ; % Displacement of the right side along x
BC{1} = {BC_1,BC_2} ;

% Compression yy
BC_1 = sprintf(Format_1,1,3,2,2) ; % Fix the Bottom side along y 
BC_2 = sprintf(Format_2,2,4,2,2,-Displacement) ; % Displacement of the top side along y 
BC{2} = {BC_1,BC_2} ;

% Compression zz
BC_1 = sprintf(Format_1,1,5,3,3) ; % Fix the front side  along z
BC_2 = sprintf(Format_2,2,6,3,3,-Displacement) ; % Displacement of the back side along z
BC{3} = {BC_1,BC_2} ;

% Cisaillement xy
BC_1 = sprintf(Format_1,1,3,1,1) ; % Fix the Bottom side  along x
BC_2 = sprintf(Format_1,2,3,2,2) ; % Fix the Bottom side  along y
BC_3 = sprintf(Format_2,3,4,2,2,0) ; % Fix the top side  along y
BC_4 = sprintf(Format_2,4,4,1,1,Displacement) ; % Displacement of the top side along x 
BC{4} = {BC_1,BC_2,BC_3,BC_4} ;

% Cisaillement xz
BC_1 = sprintf(Format_1,1,5,1,1) ; % Fix the Front side  along x
BC_2 = sprintf(Format_1,2,5,3,3) ; % Fix the Front side  along z
BC_3 = sprintf(Format_2,3,6,3,3,0) ; % Fix the Back side  along z
BC_4 = sprintf(Format_2,4,6,1,1,Displacement) ; % Displacement of the Back side along x 
BC{5} = {BC_1,BC_2,BC_3,BC_4} ;

% Cisaillement yz
BC_1 = sprintf(Format_1,1,3,3,3) ; % Fix the Bottom side  along z
BC_2 = sprintf(Format_1,2,3,2,2) ; % Fix the Bottom side  along y
BC_3 = sprintf(Format_2,3,4,2,2,0) ; % Fix the Top side  along y
BC_4 = sprintf(Format_2,4,4,3,3,Displacement) ; % Displacement of the Top side along z
BC{6} = {BC_1,BC_2,BC_3,BC_4} ;

%% %%%%%%%%%%%%%%%%%%%%%%%%% LAUNCH THE PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%
initialFiles = {dir(pwd).name}; 

for i_Test = 1:nbTests
    Inp_File_Name = 'Input_File.inp' ;     % Create INP file 
    fun_Inp(Inp_File_Name,Nodes,Elements,Elements_Sets,Nodes_Sets,Young_Modulus_List,Poisson_Ratio,BC{i_Test})
    Py_File_Name = 'Script.py' ; % Create python file 
    fun_Py(Py_File_Name,Inp_File_Name)
    
    % Launch abaqus calculations 
    system(strrep('abaqus cae noGUI=Script.py','Script.py',Py_File_Name))
    
    % Results extraction 
    system('abaqus cae noGUI=Outputs.py')
    Stress{i_Test} = table2array(readtable('Stress.csv')) ;  % Store the Stress in the stress cell array
    Strain{i_Test} = table2array(readtable('Strain.csv')) ;  % Store the Strain in the stress cell array
    RF{i_Test} = table2array(readtable('RF.csv')) ; % Store the Reaction forces in the stress cell array
    Volume = table2array(readtable('Volume.csv')) ; % Store the Volumes in the stress cell array
    
    % Delete temporary files 
    newFiles = setdiff({dir(pwd).name}, initialFiles);
    FilesToDelete = newFiles(~endsWith(newFiles, '.py') & ~endsWith(newFiles, '.odb'));
    delete(FilesToDelete{:});
end

%% %%%%%%%%%%%%%%%     STIFFNESS MATRIX Calculation      %%%%%%%%%%%%%%%%%%

Sigma = cell2mat(cellfun(@(x) (sum(x.*repmat(Volume,1,6),1)./sum(Volume))',Stress,'un',0)) ; 
Epsilon = cell2mat(cellfun(@(x) (sum(x.*repmat(Volume,1,6),1)./sum(Volume))',Strain,'un',0)) ;
det(Epsilon)
E_App = Sigma/Epsilon ;
E_App_Sym = 0.5*(E_App+E_App'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%% Save the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Symetric Stiffness Matrix : ');disp(E_App_Sym);
disp('---------------------------');Mesh.E_app_Struct.E_App_Sym = E_App_Sym;Mesh.E_app_Struct.E_app = E_App;
disp('-------------------------');

end