
function fun_Py(File_Name,Inp_File_Name)
% Python file creation 
fid = fopen(File_Name,'wt') ;

%%  File Modifications 
% Headers
fprintf(fid,['from part import *',char(10),'from material import *',char(10),'from section import *',char(10),'from assembly import *',char(10),'from step import *',char(10),'from interaction import *',char(10),'from load import *',char(10),'from mesh import *',char(10),'from optimization import *',char(10),'from job import *',char(10),'from sketch import *',char(10),'from visualization import *',char(10),'from connectorBehavior import *',char(10)]) ;

% Import Inp file
fprintf(fid,strrep([strrep('mdb.ModelFromInputFile(inputFileName=Input_File.inp','Input_File.inp',['''',pwd,'\',Inp_File_Name,'''']),', name=''Input_File'')',char(10)],'\','/')) ; 

% Outputs request
fprintf(fid,['mdb.models[''Input_File''].fieldOutputRequests[''F-Output-1''].setValues(variables=(''S'', ''PE'', ''PEEQ'', ''PEMAG'', ''LE'', ''U'', ''RF'', ''CF'', ''CSTRESS'', ''CDISP'', ''IVOL''))',char(10)]) ;

% Execute job
fprintf(fid,['mdb.Job(atTime=None, contactPrint=OFF, description='''', echoPrint=OFF, explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model=''Input_File'', modelPrint=OFF, multiprocessingMode=DEFAULT, name=''Job-1'', nodalOutputPrecision=SINGLE, numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='''', type=ANALYSIS, userSubroutine='''', waitHours=0, waitMinutes=0)',char(10)]) ; 
fprintf(fid,['mdb.jobs[''Job-1''].submit(consistencyChecking=OFF)',char(10)]) ; 
end