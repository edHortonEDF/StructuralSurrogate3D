# -*- coding: utf-8 -*-
"""
This is the primary script that should be executed to run the surrogate model.

All the inputs required are located in this file unless otherwise stated.

Created on Fri Jun 20 08:14:55 2025

@author: horto1e
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rom_utilities as util
import stiffnessMatrix as sm
import MaterialModel as matmod

import NewtonRhapson_generic as NR

# ========= INPUTS ============
#Open props file

props_fname = 'Units.csv'
#applied conditions
load_fname = 'Loading.csv'
constraint_fname = 'Constraint.csv'
temp_fname = 'UnitTemps_0.csv'
path = 'C:\\Temp\\HLoop\\3D\\SimpleTests\\shear\\'
#path = ''


model_type = "3D" # Options: "1D", "3D"
#This removes some components from the elastic stiffness matrix. 
#Only add constraints/applied alods in Fxx,Mxx. Loads in other directions will cause model to not converge.

#constraint type - determines which D_e and B matrix to use.
#PE - plane strain - also need to set uXzz to 0 in constraints file.
#PS - plane stress (not implemented)
#Axi - axisymmetric - (not implemented)
constraint_type = 'PE'

#Testing
Testing = 0
results_fname = "expected_results.txt"

#Jacobian Method - Currently not implemented. CHange in NR_generic method
NumericalJ = 0
pert = 1e-8

 

# ========= EDND OF INPUTS ============
nodal_dofs = 8 # number of degrees of freedom at each node. This should only change is the model is fundamentally altered to have more
#dofs per node. Note that each element will have nodal_dofs*2 degrees of freedom.
# uxx,uxy,uxz,rxz,uyy,uzz,r_yy,r_zz

if Testing == 1:
    fea = pd.read_csv(path+results_fname,delimiter='\t')
    fea = pd.DataFrame.to_numpy(fea)  
    
#Open properties file
properties = util.openPropsFile(path,props_fname,nodal_dofs,model_type,constraint_type)

#Open loading file
load_dict = util.openLoadingFile(properties,path,load_fname,constraint_fname,temp_fname)




#generate stiffness matrix.
stiff = sm.generateStiffnessMatrix(properties)


n_t_steps = len(load_dict["t"])
n_units = len(properties)-1
n_nodes = int(len(stiff)/nodal_dofs)
total_dofs = len(stiff)

#Initialise stuff
u = np.zeros([n_t_steps,n_nodes*nodal_dofs])
F = np.zeros([n_t_steps,n_nodes*nodal_dofs])

#This should be hardcoded as 6. (number of stress components.)
stress_u = np.zeros([n_t_steps,6,n_units])
strain_u = np.zeros([n_t_steps,6,n_units])
strain_u_th = np.zeros([n_t_steps,6,n_units])
strain_u_el = np.zeros([n_t_steps,6,n_units])



#step through time increments

for t_step in range(0,n_t_steps):
    
    
    
    #Calculate new dn (distance of each unit's connections from the neutral axis of the node.) based on shift of neutral axis. neutral axis shifts due to changing D_e due to plasticity.
    
    #Calculate a and b matrix for each unit. - do it here so that if tn shifts a and b must be recalculated (dn changes).
    if constraint_type != 'Axi':
        properties = sm.aAndbMatrixCalc(properties,nodal_dofs) #No yet implemented correctly
    else:
        properties = sm.aAndbMatrixCalc_axi(properties,nodal_dofs)
    
    print('Step:',t_step+1,'/',n_t_steps)
    #current applied load/constraint
    F_app =load_dict["F_app"][t_step,:]
    u_app =load_dict["u_app"][t_step,:]
    
    
    
    #previous applied load/constraint - assume start at 0.
    if t_step == 0:
        F_app_prev = np.zeros([total_dofs])
        u_app_prev = np.zeros([total_dofs])
        u_prev = np.zeros([n_nodes*nodal_dofs])
    else:
        F_app_prev = load_dict["F_app"][t_step-1,:]
        u_app_prev = load_dict["u_app"][t_step-1,:]
        u_prev = u[t_step-1,:]
        
    #update T and dT in state variables
    #T is temperature at the previous increment (same as everything else in state var.)
    for unit in range(0,n_units):
        
        if t_step == 0:
            properties[str(unit)]["StateVariables"]["T"] = load_dict["T_app"][t_step,unit]
            properties[str(unit)]["StateVariables"]["dT"] = 0
        else:
            properties[str(unit)]["StateVariables"]["T"] = load_dict["T_app"][t_step-1,unit]
            properties[str(unit)]["StateVariables"]["dT"] = load_dict["T_app"][t_step,unit] - load_dict["T_app"][t_step-1,unit]
        
    
        
    #==== Load increment ====
    
    #Applied load increment:
    dF_app = F_app - F_app_prev
    
    #Applied displacements
    du_app = [None]*total_dofs
    for u_it in range(0,total_dofs):
        #If current step is N then N
        if u_app[u_it] == 'N':
            du_app[u_it] = 'N'
        #If previous step was 'N' and current is a number then just that u. This is not added to currently calculated stuff, it is absolute.
        elif u_app[u_it] != 'N' and u_app_prev[u_it] == 'N':
            du_app[u_it] = u_app[u_it]
        #Both prev and current are number then u1-u0
        elif u_app[u_it] != 'N' and u_app_prev[u_it] != 'N':
            du_app[u_it] = u_app[u_it]-u_app_prev[u_it]
    
    
    dF_rf_app = np.zeros([total_dofs])
    # Target reaction force forces
    for i in range(0,total_dofs):
        if du_app[i] == 'N':
            dF_rf_app[i] = dF_app[i]
        else:
            dF_rf_app[i] = 0
    
    
    #Initial Guess accounting for constraints (applied displacements and rotations) - u_0
    du_g = np.zeros([total_dofs])
    for i in range(0,total_dofs):
        if du_app[i] == 'N':
            du_g[i] = 0
        else:
            du_g[i] = du_app[i]
    
    
    
    #warning about applied force + constraint on same dof
    for dof in range(0,len(du_app)):
        if du_app[dof] != 'N' and dF_app[dof] != 0:
            print('WARNING: Both force and constraint applied to node '+str(dof)+'! Only constraint will be applied')
    
    #Input = du_g,dF_rf_app
    #output = du_g_converged
    
    #Newton Raphson for full bar config
    du_fin,properties  = NR.NewtonRaphson(du_g,dF_rf_app,du_app,properties) 
   
    #
    #If you want the final forces. Not sure this is correct...
    dF_rf_final = np.dot(stiff,du_fin)
    if t_step == 0:
        u[t_step,:] = du_fin
        F[t_step,:] = dF_rf_final*-1 # THis needs to change later
        
    else:
        u[t_step,:] =  u[t_step-1,:]+du_fin
        F[t_step,:] = F[t_step-1,:] + dF_rf_final*-1
        
    

        
    #update state variables
    #will need to be model specific. - currentl assumes a conventional material model that uses PEEQ and backstress to track plasticity.
    for unit in range(0,n_units):
        #total strain
        properties[str(unit)]["StateVariables"]["eto"] +=  properties[str(unit)]["StateVariables"]["deto"]
        properties[str(unit)]["StateVariables"]["deto"] = np.zeros([6])
        #elastic strain
        properties[str(unit)]["StateVariables"]["eel"] +=  properties[str(unit)]["StateVariables"]["deel"]
        properties[str(unit)]["StateVariables"]["deel"] = np.zeros([6])
        #plastic strain
        properties[str(unit)]["StateVariables"]["ep"] +=  properties[str(unit)]["StateVariables"]["dep"]
        properties[str(unit)]["StateVariables"]["dep"] = np.zeros([6])
        #thermal strain
        properties[str(unit)]["StateVariables"]["eth"] +=  properties[str(unit)]["StateVariables"]["deth"]
        properties[str(unit)]["StateVariables"]["deth"] = np.zeros([6])
        #stress
        properties[str(unit)]["StateVariables"]["stress"] +=  properties[str(unit)]["StateVariables"]["dstress"]
        properties[str(unit)]["StateVariables"]["dstress"] = np.zeros([6])
        #alpha
        properties[str(unit)]["StateVariables"]["alpha"] +=  properties[str(unit)]["StateVariables"]["dalpha"]
        properties[str(unit)]["StateVariables"]["dalpha"] = np.zeros([6])
           
    #stress evolution in each unit
    for unit in range(0,n_units):
        stress_u[t_step,:,unit] = properties[str(unit)]["StateVariables"]["stress"]
        strain_u[t_step,:,unit] = properties[str(unit)]["StateVariables"]["eto"]
        strain_u_th[t_step,:,unit] = properties[str(unit)]["StateVariables"]["eth"]
        strain_u_el[t_step,:,unit] = properties[str(unit)]["StateVariables"]["eel"]
    
    
# END OF SIMULATION LOOP 

#Post processing:
    
#mises stress
mises = np.zeros([n_t_steps,n_units])
hyd =np.zeros([n_t_steps,n_units])

for t in range(0,n_t_steps):
    for unit in range(0,n_units):
        s_vec = stress_u[t,:,unit]
        mises[t,unit],hyd[t,unit] = util.misescalc(s_vec)

print('==============')
print('Check the following before analysing results:')
print('1. additional constraints are on/off as required')
print('2. Check material models are correct in file')
print('3. Check applied loads and constraints are correct in files')
print('4. Check initial neutral axis values are correct. (This will be automated eventually)')
print('==============')

print('Model Type: ',model_type)

print('STRESS')
print(stress_u[-1,:])
print('MISES')
print(mises[-1,:])
print('STRAIN')
print(strain_u[-1,:])

plt.figure(dpi=200)
plt.rcParams["figure.figsize"] = (5,5)
if Testing == 1:
    plt.plot(fea[:,0],fea[:,1])
for i in range(0,n_units):
    #Sxx (direct)
    
    plt.plot(strain_u[:,0,i],stress_u[:,0,i],'x')
    
plt.ylabel('Stress [MPa]')
plt.xlabel('Strain')
#plt.legend(['Target','A','B','C','D','E','F'])

marker = ['.','+','x','1','2','3']

plt.figure(dpi=200)
#stress-strain in unit 1 - all components
for i in range(0,6):
    #plt.plot(strain_u[:,i,0],stress_u[:,i,0],'-')
    plt.plot(load_dict["t"],stress_u[:,i,0],marker[i])
    
plt.plot(load_dict["t"],mises[:,0])
plt.legend(['S11','S22','S33','S12','S13','S23','Mises'])
plt.title('Unit 1')
plt.xlabel('time')
plt.ylabel('Stress [MPa]')

    

plt.figure(dpi=200)
plt.plot(load_dict["t"],mises[:,0])
plt.plot(load_dict["t"],hyd[:,0])

plt.ylabel('Stress')
plt.xlabel('time')  
plt.legend(['Mises','Hydrostatic'])
plt.title('Mises and hydrostatic stress in unit 1') 


marker = ['.-','-+','-x','1','2','3']

broken = 0
if broken == 1:
    #for every dof: - onyl last value of stress
    nodes = np.unique(properties["Raw"][:,1:3])
    for n in range(0,n_nodes):
        node = nodes[n]
        distance = []
        
        units_node = 0
        for unit in range(0,n_units):
            if properties["Raw"][unit,1] == node :
                distance.append(properties["Raw"][unit,3])
                units_node += 1
                
            elif properties["Raw"][unit,2] == node:
                distance.append(properties["Raw"][unit,4])
                units_node += 1
                 
        stress = np.zeros([6,units_node])
        strain = np.zeros([6,units_node])    
        
        u_1 = 0
        for unit in range(0,n_units):
            if properties["Raw"][unit,1] == node or properties["Raw"][unit,2] == node:
                stress[:,u_1] = stress_u[-1,:,unit]
                strain[:,u_1] = strain_u[-1,:,unit]
                u_1 += 1
    
        plt.figure(dpi=200)
        for i in range(0,1):
            plt.plot(distance,stress[i,:],marker[i])
        plt.title('Node:' + str(node))
        plt.legend(['S11','S22','S33','S12','S13','S23','Mises'])    
        plt.ylabel('Stress')
        plt.xlabel('distance from loading position')         
           
        plt.figure(dpi=200)
        for i in range(0,1):
            plt.plot(distance,strain[i,:],marker[i])
        plt.title('Node:' + str(node))
        plt.legend(['E11','E22','E33','E12','E13','E23','Mises'])    
        plt.ylabel('Strain')
        plt.xlabel('distance from loading position')         
            
#stress

plt.figure(dpi=200)
for i in range(0,6):
    plt.plot(properties["Raw"][:,3],stress_u[-1,i,:],marker[i])
plt.plot(properties["Raw"][:,3],mises[-1,:],'k-')
plt.legend(['S11','S22','S33','S12','S13','S23','Mises'])
plt.ylabel('Stress')
plt.xlabel('distance from loading position')

strains_on = 1
if strains_on == 1:
    #strain
    marker = ['.-','-+','-x','1','2','3']
    plt.figure(dpi=200)
    for i in range(0,6):
        plt.plot(properties["Raw"][:,3],strain_u[-1,i,:],marker[i])
    #plt.plot(properties["Raw"][:,3],mises[-1,:],'k-')
    plt.legend(['E11','E22','E33','E12','E13','E23'])
    plt.ylabel('Strain')
    plt.xlabel('distance from loading position') 
 
    
 
strains_therm_on = 0
if strains_therm_on == 1:    
    
    leg_list = [['E11 - total','E11 - elastic','E11 - thermal'],
                ['E22 - total','E22 - elastic','E22 - thermal'],
                ['E33 - total','E33 - elastic','E33 - thermal']]
    
    
    #strain components at t=1 for E11
    marker = ['.-','-+','-x','1','2','3']
    for f in range(0,3):
        plt.figure(dpi=200)
    
        plt.plot(properties["Raw"][:,3],strain_u[-1,f,:],marker[0])
        plt.plot(properties["Raw"][:,3],strain_u_el[-1,f,:],marker[1])
        plt.plot(properties["Raw"][:,3],strain_u_th[-1,f,:],marker[2])
        #plt.plot(properties["Raw"][:,3],mises[-1,:],'k-')
        plt.legend(leg_list[f])
        plt.ylabel('Strain')
        plt.xlabel('distance from loading position') 



   
