function Mesh_Tet = fun_FV_Tet(ROI)
%ROI = Data.ROI_2;

Bins = bwconncomp(ROI,6) ;
[~,idx] = max(cell2mat(cellfun(@(x) size(x,1),Bins.PixelIdxList,'un',0))) ;
ROI(setdiff(1:numel(ROI),Bins.PixelIdxList{idx})) = 0 ;
[x,y,z] = ind2sub(size(ROI),find(ones(size(ROI)))) ;
Grid = [x,y,z] ;
% Création de la surface triangulée
[x,y,z] = ndgrid(1:size(ROI,1),1:size(ROI,2),1:size(ROI,3)) ;
[Faces,Vertices] = MarchingCubes(x,y,z,ROI,0.5) ;
Faces = double(Faces) ; Vertices = double(Vertices) ;
% Extraction des bords libres
Free_Edges = fun_Free_Edges(Faces) ;

% Affichage
figure(1) ; set(gcf,'color','w') ; axis equal ; hold on ; view(25,10)
patch('faces',Faces,'vertices',Vertices,'facecolor','cyan')
scatter3(Grid(:,1),Grid(:,2),Grid(:,3),'filled','k')
line([Vertices(Free_Edges(:,1),1),Vertices(Free_Edges(:,2),1)]',[Vertices(Free_Edges(:,1),2),Vertices(Free_Edges(:,2),2)]',[Vertices(Free_Edges(:,1),3),Vertices(Free_Edges(:,2),3)]','linewidth',2,'color','r')

%% Triangularisation des bords
[Nodes,Connectivity_Bone,Connectivity_Marrow] = deal(cell(6,1)) ;
for i_Face = 1:6
    if mod(i_Face,2)
    idx_Vertices = find(Vertices(:,ceil(i_Face/2))==min(Vertices(:,ceil(i_Face/2)))) ;
    idx_Grid = find(Grid(:,ceil(i_Face/2))==min(Grid(:,ceil(i_Face/2)))) ;
    else
    idx_Vertices = find(Vertices(:,ceil(i_Face/2))==max(Vertices(:,ceil(i_Face/2)))) ;
    idx_Grid = find(Grid(:,ceil(i_Face/2))==max(Grid(:,ceil(i_Face/2)))) ; 
    end
    Points = [Vertices(idx_Vertices,:) ; Grid(idx_Grid,:)] ;
    % Bords libres appartenant à la face traitée
    Free_Edges_Face = Free_Edges(sum(ismember(Free_Edges,idx_Vertices),2)==2,:) ; 
    Free_Edges_temp = zeros(size(Free_Edges_Face)) ;
    for i_Edge = 1:size(Free_Edges_Face,1)
        [~,b] = ismember(Vertices(Free_Edges_Face(i_Edge,:),:),Points,'rows') ;
        Free_Edges_temp(i_Edge,:) = b' ;
    end
    Free_Edges_Face = Free_Edges_temp ;
    % Bords correspondants aux arrêtes du ROI cubique
    Points_temp = Points(:,setdiff(1:3,ceil(i_Face/2))) ;
    Boundary = unique([find(sum(Points_temp==min(Points_temp(:)),2)) ; find(sum(Points_temp==max(Points_temp(:)),2))],'rows') ;
    [~,Order] = sort(atan2(Points_temp(Boundary,2)-mean(Points_temp(Boundary,2)),Points_temp(Boundary,1)-mean(Points_temp(Boundary,1)))) ;
    Boundary = Boundary(Order) ;  
    Free_Edges_Boundary = Boundary([1:size(Boundary,1) ; [2:size(Boundary,1),1]]') ;
    % Définition des polygones
    [a,b] = groupcounts(Free_Edges_Face(:)) ;
    idx_Extremities = find(ismember(Boundary,b(a==1))) ;
    idx_Extremities = [idx_Extremities ; idx_Extremities(1)] ;
    Edges_Polygone = cell(size(idx_Extremities,1)-1,1) ;
    for i_Polygon = 1:size(idx_Extremities,1)-1
        if i_Polygon==size(idx_Extremities,1)-1
        Edges_Polygone{i_Polygon} = [Free_Edges_Boundary(idx_Extremities(i_Polygon):size(Free_Edges_Boundary,1),:) ; Free_Edges_Boundary(1:idx_Extremities(i_Polygon+1)-1,:)] ;
        else
        Edges_Polygone{i_Polygon} = Free_Edges_Boundary(idx_Extremities(i_Polygon):idx_Extremities(i_Polygon+1)-1,:) ;
        end
    end     
    % Sélection d'un segment sur deux et ajout des bords libres de la face
    if  i_Face == 1 || i_Face == 3 || i_Face == 6
        Edges_Polygone = [Free_Edges_Face ; cat(1,Edges_Polygone{2:2:end})] ;
    else 
        Edges_Polygone = [Free_Edges_Face ; cat(1,Edges_Polygone{1:2:end})] ;
    end
    % Triangularisation sous contraintes
    DT = delaunayTriangulation(Points_temp,Edges_Polygone) ;
    Connectivity_Bone{i_Face} = DT.ConnectivityList(isInterior(DT),:) ;
    Connectivity_Marrow{i_Face} = DT.ConnectivityList(~isInterior(DT),:) ;
    Nodes{i_Face} = Points ;
    % % Affichage
    figure(1) ; set(gcf,'color','w') ; axis equal off ; hold on ; view(-90,0)
    patch('faces',Faces,'vertices',Vertices,'facecolor','cyan')
    scatter3(Points(:,1),Points(:,2),Points(:,3),'filled')
    line([Points(Free_Edges_Face(:,1),1),Points(Free_Edges_Face(:,2),1)]',[Points(Free_Edges_Face(:,1),2),Points(Free_Edges_Face(:,2),2)]',[Points(Free_Edges_Face(:,1),3),Points(Free_Edges_Face(:,2),3)]','linewidth',2,'color','r')
    scatter3(Points(Boundary,1),Points(Boundary,2),Points(Boundary,3),'filled','k')
    line([Points(Free_Edges_Boundary(:,1),1),Points(Free_Edges_Boundary(:,2),1)]',[Points(Free_Edges_Boundary(:,1),2),Points(Free_Edges_Boundary(:,2),2)]',[Points(Free_Edges_Boundary(:,1),3),Points(Free_Edges_Boundary(:,2),3)]','linewidth',2,'color','k')
    line([Points(Edges_Polygone(:,1),1),Points(Edges_Polygone(:,2),1)]',[Points(Edges_Polygone(:,1),2),Points(Edges_Polygone(:,2),2)]',[Points(Edges_Polygone(:,1),3),Points(Edges_Polygone(:,2),3)]','linewidth',2,'color','b')
    patch('faces',Connectivity_Bone{i_Face},'Vertices',Points,'facecolor','cyan')
    patch('faces',Connectivity_Marrow{i_Face},'Vertices',Points,'facecolor','y')  
end

%%
Connectivity_Bone([1,4,5]) = cellfun(@(x) flip(x,2),Connectivity_Bone([1,4,5]),'un',0) ;
Connectivity_Marrow([1,4,5]) = cellfun(@(x) flip(x,2),Connectivity_Marrow([1,4,5]),'un',0) ;
% Assemblage des tables de connectivité
Sizes = cell2mat(cellfun(@(x) size(x,1),Nodes,'un',0)) ; Nodes = cat(1,Nodes{:}) ;
Sizes = num2cell(cumsum([0 ; Sizes(1:end-1)])) ;
Connectivity_Bone = cellfun(@(x,y) x+y,Connectivity_Bone,Sizes,'un',0) ; Connectivity_Bone = cat(1,Connectivity_Bone{:}) ;
Connectivity_Marrow = cellfun(@(x,y) x+y,Connectivity_Marrow,Sizes,'un',0) ; Connectivity_Marrow = cat(1,Connectivity_Marrow{:}) ;
% Assemblage avec la surface osseuse
Faces_Bone = [flip(Faces,2) ; Connectivity_Bone+size(Vertices,1)] ; Vertices_Bone = [Vertices ; Nodes] ;
Faces_Marrow = [Faces ; Connectivity_Marrow+size(Vertices,1)] ; Vertices_Marrow = [Vertices ; Nodes] ;

% Affichage
% figure(1) ; set(gcf,'color','w') ; axis equal ; hold on ; view(45,25)
% patch('faces',Faces_Bone,'vertices',Vertices_Bone,'facecolor','cyan','facealpha',1)
% patch('faces',Connectivity_Bone,'vertices',Nodes,'facecolor','cyan','facealpha',1)
% Normals_Bone = faceNormal(triangulation(Faces_Bone,Vertices_Bone)) ;
% Centers_Bone = incenter(triangulation(Faces_Bone,Vertices_Bone)) ;
% quiver3(Centers_Bone(:,1),Centers_Bone(:,2),Centers_Bone(:,3),Normals_Bone(:,1),Normals_Bone(:,2),Normals_Bone(:,3),'color','r') ;
% patch('faces',Faces_Marrow,'vertices',Vertices_Marrow,'facecolor','y','facealpha',1)
% Normals_Marrow = faceNormal(triangulation(Faces_Marrow,Vertices_Marrow)) ;
% Centers_Marrow = incenter(triangulation(Faces_Marrow,Vertices_Marrow)) ;
% quiver3(Centers_Marrow(:,1),Centers_Marrow(:,2),Centers_Marrow(:,3),Normals_Marrow(:,1),Normals_Marrow(:,2),Normals_Marrow(:,3),'color','r') ;

%% %%%%%%%%%%%%%%%%% FILTER NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bone_Nodes_Index = unique(Faces_Bone(:));
Bone_Nodes = Vertices_Bone(Bone_Nodes_Index,:);

Marrow_Nodes_Index = unique(Faces_Marrow(:));
Marrow_Nodes = Vertices_Bone(Marrow_Nodes_Index,:);

Mesh_Tet.Nodes = Vertices_Bone;
Mesh_Tet.Bone_Nodes = Bone_Nodes;
Mesh_Tet.Marrow_Nodes = Marrow_Nodes;
Mesh_Tet.Bone_Faces = Faces_Bone;
Mesh_Tet.Marrow_Faces = Faces_Marrow;

stlwrite(triangulation(Mesh_Tet.Bone_Faces,Mesh_Tet.Nodes),'Bone-175.stl')
stlwrite(triangulation(Mesh_Tet.Marrow_Faces,Mesh_Tet.Nodes),'Marrow-175.stl')

end

