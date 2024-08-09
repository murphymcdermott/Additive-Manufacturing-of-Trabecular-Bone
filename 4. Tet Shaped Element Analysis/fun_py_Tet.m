function fun_py_Tet(FEM_Tet,Test)

%% %%%%%%%%%%%%%%%%%%%%%%%   VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROI_inp = FEM_Tet.ROI_inp;
Complement_inp = FEM_Tet.Complement_inp;
Sample_Name = FEM_Tet.Sample_Name;
Material1 = FEM_Tet.Material1; Matrial1_Youngs_Modulus = FEM_Tet.Matrial1_Youngs_Modulus; Matrial1_Poisson_Coeff = FEM_Tet.Matrial1_Poisson_Coeff;
Material2 = FEM_Tet.Material2; Matrial2_Youngs_Modulus = FEM_Tet.Matrial2_Youngs_Modulus; Matrial2_Poisson_Coeff = FEM_Tet.Matrial2_Poisson_Coeff;
Sets_name = FEM_Tet.Sets_name;
Set_range = FEM_Tet.Set_range; %minX, maxX, minY, maxY, minZ, maxZ
Resolution = FEM_Tet.Resolution;
Displacement = FEM_Tet.Displacement;

%% %%%%%%%%%%%%%%%% CREATE THE PYTHON FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Py_File_Name = [Test,'.py'] ;
file_id = fopen(Py_File_Name,'wt') ;

%% %%%%%%%%%%%%%%%%%%%%% EDIT PYTHON FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Headers 
fprintf(file_id,"from part import *\nfrom material import *\nfrom section import *\nfrom assembly import *\nfrom step import *\nfrom interaction import *\nfrom load import *\nfrom mesh import *\nfrom optimization import *\nfrom job import *\nfrom sketch import *\nfrom visualization import *\nfrom connectorBehavior import *\n\n");

% Import models and parts 
fprintf(file_id,sprintf("## Import models and parts ##\nmdb.Model(modelType=STANDARD_EXPLICIT, name='%s')\n",Sample_Name));
fprintf(file_id,sprintf("mdb.ModelFromInputFile(inputFileName='%s', name='Complement-model')\n",Complement_inp));
fprintf(file_id,"mdb.models['Complement-model'].parts['PART-1'].mergeNodes(nodes=mdb.models['Complement-model'].parts['PART-1'].nodes, tolerance=1e-06)\nmdb.models['Complement-model'].parts.changeKey(fromName='PART-1', toName='Complement')\n");
fprintf(file_id,sprintf("mdb.ModelFromInputFile(inputFileName='%s', name='ROI-model')\n",ROI_inp));
fprintf(file_id,"mdb.models['ROI-model'].parts.changeKey(fromName='PART-1', toName='Bone')\nmdb.models['ROI-model'].parts['Bone'].mergeNodes(nodes=mdb.models['ROI-model'].parts['Bone'].nodes, tolerance=1e-06)\n");
fprintf(file_id,sprintf("mdb.models['%s'].Part('Bone', mdb.models['ROI-model'].parts['Bone'])\nmdb.models['%s'].Part('Complement', mdb.models['Complement-model'].parts['Complement'])\n",Sample_Name,Sample_Name));

% Delete other parts
fprintf(file_id,"\n## Delete other parts ##\ndel mdb.models['ROI-model']\ndel mdb.models['Model-1']\ndel mdb.models['Complement-model']\n");
% Tri to tet
fprintf(file_id,sprintf("\n## Tri to tet ##\nmdb.meshEditOptions.setValues(enableUndo=True, maxUndoCacheElements=0.5)\nmdb.models['%s'].parts['Bone'].generateMesh(elemShape=TET)\nmdb.models['%s'].parts['Complement'].generateMesh(elemShape=TET)\n",Sample_Name,Sample_Name));

% Properties and Section assignement
fprintf(file_id,"\n## Properties and Section assignement ##\n");
fprintf(file_id,sprintf("mdb.models['%s'].Material(name='%s')\nmdb.models['%s'].materials['%s'].Elastic(table=((%2.f, %f), ))\n",Sample_Name,Material1,Sample_Name,Material1,Matrial1_Youngs_Modulus,Matrial1_Poisson_Coeff));
fprintf(file_id,sprintf("mdb.models['%s'].Material(name='%s')\nmdb.models['%s'].materials['%s'].Elastic(table=((%2.f, %f), ))\n",Sample_Name,Material2,Sample_Name,Material2,Matrial2_Youngs_Modulus,Matrial2_Poisson_Coeff));
fprintf(file_id,sprintf("mdb.models['%s'].HomogeneousSolidSection(material='%s', name='%s-section', thickness=None)\nmdb.models['%s'].HomogeneousSolidSection(material='%s', name='%s-section', thickness=None)\n",Sample_Name,Material1,Material1,Sample_Name,Material2,Material2));

% Assign sets to sections 
fprintf(file_id,"\n## Assign sets to sections  ##\n");
fprintf(file_id,sprintf("mdb.models['%s'].parts['Complement'].Set(elements=mdb.models['%s'].parts['Complement'].elements[:],name='%s-Nodes')\nmdb.models['%s'].parts['Complement'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=mdb.models['%s'].parts['Complement'].sets['%s-Nodes'], sectionName='%s-section', thicknessAssignment=FROM_SECTION)\n",Sample_Name,Sample_Name,Material2,Sample_Name,Sample_Name,Material2,Material2));
fprintf(file_id,sprintf("mdb.models['%s'].parts['Bone'].Set(elements=mdb.models['%s'].parts['Bone'].elements[:],name='%s-Nodes')\nmdb.models['%s'].parts['Bone'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=mdb.models['%s'].parts['Bone'].sets['%s-Nodes'], sectionName='%s-section', thicknessAssignment=FROM_SECTION)\n",Sample_Name,Sample_Name,Material1,Sample_Name,Sample_Name,Material1,Material1));

% Assembly
fprintf(file_id,"\n## Assembly  ##\n");
fprintf(file_id,strrep("mdb.models['%s'].rootAssembly.DatumCsysByDefault(CARTESIAN)\nmdb.models['%s'].rootAssembly.Instance(dependent=ON, name='Bone-1', part=mdb.models['%s'].parts['Bone'])\nmdb.models['%s'].rootAssembly.Instance(dependent=ON, name='Complement-1', part=mdb.models['%s'].parts['Complement'])\nmdb.models['%s'].rootAssembly._previewMergeMeshes(instances=(mdb.models['%s'].rootAssembly.instances['Bone-1'], mdb.models['%s'].rootAssembly.instances['Complement-1']), nodeMergingTolerance=1e-06)\n",'%s',Sample_Name));
fprintf(file_id,strrep("mdb.models['%s'].rootAssembly.InstanceFromBooleanMerge(domain=MESH, instances=(mdb.models['%s'].rootAssembly.instances['Bone-1'], mdb.models['%s'].rootAssembly.instances['Complement-1']), mergeNodes=BOUNDARY_ONLY, name='Combined-part', nodeMergingTolerance=1e-06, originalInstances=SUPPRESS)\n",'%s',Sample_Name));

% Sets
fprintf(file_id,"\n## Sets  ##\n");
set_format = "mdb.models['%s'].rootAssembly.Set(name='%s', nodes=mdb.models['%s'].rootAssembly.instances['Combined-part-1'].nodes.getByBoundingBox(xMin=%2.f, xMax=%2.f, yMin=%2.f, yMax=%2.f, zMin=%2.f, zMax=%2.f))\n";
fprintf(file_id,sprintf((set_format),Sample_Name,Sets_name{1},Sample_Name,Set_range(1,1),Set_range(1,1)+Resolution,Set_range(2,1),Set_range(2,2),Set_range(3,1),Set_range(3,2)));
fprintf(file_id,sprintf((set_format),Sample_Name,Sets_name{2},Sample_Name,Set_range(1,2)-Resolution,Set_range(1,2),Set_range(2,1),Set_range(2,2),Set_range(3,1),Set_range(3,2)));
fprintf(file_id,sprintf((set_format),Sample_Name,Sets_name{3},Sample_Name,Set_range(1,1),Set_range(1,2),Set_range(2,1),Set_range(2,1)+Resolution,Set_range(3,1),Set_range(3,2)));
fprintf(file_id,sprintf((set_format),Sample_Name,Sets_name{4},Sample_Name,Set_range(1,1),Set_range(1,2),Set_range(2,2)-Resolution,Set_range(2,2),Set_range(3,1),Set_range(3,2)));
fprintf(file_id,sprintf((set_format),Sample_Name,Sets_name{5},Sample_Name,Set_range(1,1),Set_range(1,2),Set_range(2,1),Set_range(2,2),Set_range(3,1),Set_range(3,1)+Resolution));
fprintf(file_id,sprintf((set_format),Sample_Name,Sets_name{6},Sample_Name,Set_range(1,1),Set_range(1,2),Set_range(2,1),Set_range(2,2),Set_range(3,2)-Resolution,Set_range(3,2)));

% Step
fprintf(file_id,"\n## Step  ##\n");
fprintf(file_id,sprintf("mdb.models['%s'].StaticStep(name='%s', previous='Initial')\nmdb.models['%s'].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 'IVOL'))\n",Sample_Name, Test,Sample_Name));


% Boundary-Conditions 
fprintf(file_id,"\n## Boundary Condition  ##\n");
BC_Format = "mdb.models['%s'].DisplacementBC(amplitude=UNSET, createStepName='%s', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-%s', region=mdb.models['%s'].rootAssembly.sets['%s'], u%d=%d)\n";


if strcmp(Test, 'CompressionX')
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{1},Sample_Name,Sets_name{1},1,Displacement));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{2},Sample_Name,Sets_name{2},1,0));
    
elseif strcmp(Test, 'CompressionY')
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{3},Sample_Name,Sets_name{3},2,Displacement));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{4},Sample_Name,Sets_name{4},2,0));
   
elseif strcmp(Test, 'CompressionZ')
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{5},Sample_Name,Sets_name{5},3,Displacement));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{6},Sample_Name,Sets_name{6},3,0));
   
elseif strcmp(Test, 'ShearXY')
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{3},Sample_Name,Sets_name{3},1,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{3},Sample_Name,Sets_name{3},2,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{4},Sample_Name,Sets_name{4},2,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{4},Sample_Name,Sets_name{4},1,Displacement));
   
elseif strcmp(Test, 'ShearXZ')
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{5},Sample_Name,Sets_name{5},1,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{5},Sample_Name,Sets_name{5},3,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{6},Sample_Name,Sets_name{6},3,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{6},Sample_Name,Sets_name{6},1,Displacement));
   
elseif strcmp(Test, 'ShearYZ')
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{3},Sample_Name,Sets_name{3},3,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{3},Sample_Name,Sets_name{3},2,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{4},Sample_Name,Sets_name{4},2,0));
    fprintf(file_id,sprintf((BC_Format),Sample_Name, Test, Sets_name{4},Sample_Name,Sets_name{4},3,Displacement));

end


% Jobs 
fprintf(file_id,"\n## Job  ##\n");
fprintf(file_id,sprintf("mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model='%s', modelPrint=OFF, multiprocessingMode=DEFAULT, name='Job-%s', nodalOutputPrecision=SINGLE, numCpus=1, numGPUs=0, numThreadsPerMpiProcess=1, queue=None, resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)\n",Sample_Name,Test));
fprintf(file_id,sprintf("mdb.jobs['Job-%s'].submit(consistencyChecking=OFF)",Test));
