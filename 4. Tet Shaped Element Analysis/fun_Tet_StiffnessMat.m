function Mesh = fun_Tet_StiffnessMat(FEM_Tet)
%% %%%%%%%%%%%%%%%%%%%%%%%    VARIABLES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FEM_Tet.test_Name = {'CompressionX','CompressionY','CompressionZ','ShearXY','ShearXZ','ShearYZ'};
FEM_Tet.ROI_inp = 'C:\\Users\\mh3407\\OneDrive\\Documents\\BioFabLab\\Tet_analysis\\Bone-175.inp'; 
FEM_Tet.Complement_inp = 'C:\\Users\\mh3407\\OneDrive\\Documents\\BioFabLab\\Tet_analysis\\Marrow-175.inp'; 
FEM_Tet.Sample_Name = 'Sample-175';
FEM_Tet.Material1 = 'Bone'; FEM_Tet.Matrial1_Youngs_Modulus = 600e6; FEM_Tet.Matrial1_Poisson_Coeff = 0.3; 
FEM_Tet.Material2 = 'Void'; FEM_Tet.Matrial2_Youngs_Modulus = 10e0; FEM_Tet.Matrial2_Poisson_Coeff = 0.3; 
FEM_Tet.Sets_name = {'Left','Right','Bottom','Tom','Front','Back'};
FEM_Tet.Set_range = [min(size(Data.ROI_2(:,1))), max(size(Data.ROI_2(:,1))); min(size(Data.ROI_2(:,2))), max(size(Data.ROI_2(:,2))); min(size(Data.ROI_2(:,3))), max(size(Data.ROI_2(:,3)))];
FEM_Tet.Resolution = 1;
FEM_Tet.Displacement = 1; 
FEM_Tet.Job_names = {'Job-CompressionX','Job-CompressionY','Job-CompressionZ','Job-ShearXY','Job-ShearXZ','Job-ShearYZ'};

%% %%%%%%%%%%%%%%%%% LAUNCH THE CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialFiles = {dir(pwd).name};

for i_Test = 4:numel(FEM_Tet.test_Name)
    Test = FEM_Tet.test_Name{i_Test};
    Py_File_Name = [Test,'.py'] ;
    file_id = fopen(Py_File_Name,'wt') ;
    fun_py_Tet(FEM_Tet,Test);
    system(strrep('abaqus cae noGUI=Script.py','Script.py',Py_File_Name))
    %Modify the output 
    fun_output(FEM_Tet.Job_names{i_Test},FEM_Tet.test_Name{i_Test})
    system('abaqus cae noGUI=Output.py')
    Stress{i_Test} = table2array(readtable('Stress.csv')) ; 
    Strain{i_Test} = table2array(readtable('Strain.csv')) ; 
    RF{i_Test} = table2array(readtable('RF.csv')) ;
    Volume = table2array(readtable('Volume.csv')) ;
    % Delete temporary files 
    newFiles = setdiff({dir(pwd).name}, initialFiles);
    FilesToDelete = newFiles(~endsWith(newFiles, '.py') & ~endsWith(newFiles, '.odb'));
    delete(FilesToDelete{:});
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% E_App Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%

Sigma = cell2mat(cellfun(@(x) (sum(x.*repmat(Volume,1,6),1)./sum(Volume))',Stress,'un',0)) ;
Epsilon = cell2mat(cellfun(@(x) (sum(x.*repmat(Volume,1,6),1)./sum(Volume))',Strain,'un',0)) ;
det(Epsilon)
E_App = Sigma/Epsilon ;
E_App_Sym = 0.5*(E_App+E_App'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%% Save the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Symetric Stiffness Matrix : ');
disp(E_App_Sym);
disp('---------------------------');
Mesh.E_app_Struct.E_App_Sym = E_App_Sym;
Mesh.E_app_Struct.E_app = E_App;
disp('-------------------------');
Mesh.Stresses = Stress; Mesh.Strain = Strain; Mesh.Volumes = Volume; Mesh.RF = RF;

end