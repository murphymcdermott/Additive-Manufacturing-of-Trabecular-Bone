function fun_output(Job_name,Step_name)

Py_File_Name = 'Output.py' ;
file_id = fopen(Py_File_Name,'wt') ;

fprintf(file_id,"from odbAccess import *\n");
fprintf(file_id,"import numpy as np\n\n");

fprintf(file_id,"# Access odb\n");
fprintf(file_id,sprintf("odb = openOdb('%s.odb')\n",Job_name));
fprintf(file_id,"# Access data\n");
fprintf(file_id,sprintf("Stress_temp = odb.steps['%s'].frames[-1].fieldOutputs['S'].values\n",Step_name));
fprintf(file_id,sprintf("Strain_temp = odb.steps['%s'].frames[-1].fieldOutputs['E'].values\n",Step_name));
fprintf(file_id,sprintf("Volume_temp = odb.steps['%s'].frames[-1].fieldOutputs['IVOL'].values\n\n",Step_name));

fprintf(file_id,"Stress = np.zeros([len(Stress_temp),6])\n");
fprintf(file_id,"Strain = np.zeros([len(Strain_temp),6])\n");
fprintf(file_id,"Volume = np.zeros([len(Volume_temp),6])\n");

fprintf(file_id,"for i in range (len(Stress_temp)):\n");
fprintf(file_id,"    Stress[i,:] = Stress_temp[i].data\n");
fprintf(file_id,"    Strain[i,:] = Strain_temp[i].data\n");
fprintf(file_id,"    Volume[i] = Volume_temp[i].data\n\n");

fprintf(file_id,"# Save data\n");
fprintf(file_id,"Volume = Volume[:,0] \n");
fprintf(file_id,"np.savetxt('Stress.csv',Stress,fmt ='%%s')\n");
fprintf(file_id,"np.savetxt('Strain.csv',Strain,fmt ='%%s')\n");
fprintf(file_id,"np.savetxt('Volume.csv',Volume,fmt ='%%s')\n\n");

fprintf(file_id,"# Reaction forces\n");
fprintf(file_id,sprintf("RF_temp = odb.steps['%s'].frames[-1].fieldOutputs['RF'].values\n",Step_name));
fprintf(file_id,"RF = np.zeros([len(RF_temp),3])\n");
fprintf(file_id,"for i in range (len(RF_temp)):\n");
fprintf(file_id,"    RF[i,:] = RF_temp[i].data\n\n");

fprintf(file_id,"# Save data\n");
fprintf(file_id,"np.savetxt('RF.csv',RF,fmt ='%%s')\n");

end