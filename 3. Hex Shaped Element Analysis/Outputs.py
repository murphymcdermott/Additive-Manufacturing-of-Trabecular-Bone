from odbAccess import *
import numpy as np

# Access odb
odb = openOdb('Job-1.odb') ;
# Access data
Stress_temp = odb.steps['Step-1'].frames[-1].fieldOutputs['S'].values
Strain_temp = odb.steps['Step-1'].frames[-1].fieldOutputs['E'].values
Volume_temp = odb.steps['Step-1'].frames[-1].fieldOutputs['IVOL'].values
Stress = np.zeros([len(Stress_temp),6])
Strain = np.zeros([len(Strain_temp),6])
Volume = np.zeros([len(Volume_temp),6])
for i in range (len(Stress_temp)):
    Stress[i,:] = Stress_temp[i].data
    Strain[i,:] = Strain_temp[i].data
    Volume[i] = Volume_temp[i].data
# Save data
Volume = Volume[:,0] 
np.savetxt('Stress.csv',Stress,fmt ='%s') 
np.savetxt('Strain.csv',Strain,fmt ='%s') 
np.savetxt('Volume.csv',Volume,fmt ='%s')

# Reaction forces
RF_temp = odb.steps['Step-1'].frames[-1].fieldOutputs['RF'].values
RF = np.zeros([len(RF_temp),3])
for i in range (len(RF_temp)):
    RF[i,:] = RF_temp[i].data
# Save data
np.savetxt('RF.csv',RF,fmt ='%s')
