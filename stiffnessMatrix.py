# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 16:47:45 2025

@author: horto1e

"""

import numpy as np
import rom_utilities as util


#input:

def generateStiffnessMatrix(properties):
    
    # 3D method - 
    #1. construct local B matrix
    #2. construct local A matrix (A^T = B^T ?)
    #3. construct /unit stiffness 
    
    units_raw = properties["Raw"]    
    n_units = len(properties)-1
    
    nodes = np.unique(properties["Raw"][:,1:3])
    n_nodes = len(nodes)
    
    n_dof_unit = np.size(properties["0"]["a_mat"],axis=0) # each unit is connected to two nodes thereforce 8.
    nodal_dofs = int(n_dof_unit/2)
    n_dof_global = int(n_nodes*nodal_dofs) #ux,uy,uz,rz for each node
    

    K_global = np.zeros([n_dof_global,n_dof_global])
    
    #unit stiffnesses
    K_unit = np.zeros([n_dof_unit,n_dof_unit,n_units])    
    
#Calculate stiffnesses in every unit
    for unit in range(0,n_units):
        A = properties[str(unit)]["Area"]
        L = properties[str(unit)]["Length"]
        C_ijkl = properties[str(unit)]["C_ijkl"]
        a_mat = properties[str(unit)]["a_mat"]
        b_mat = properties[str(unit)]["b_mat"]
        
        #
        unit_connections = units_raw[unit,1:3]
        unit_dn = units_raw[unit,3:5]
        #E = units[unit,6]
        
        #a_mat and b_mat calculation here.
        
        K_unit[:,:,unit] = a_mat @ C_ijkl @ b_mat
        
        

            
    #Note that this assmbly method assumes all dofs associated with a node are grouped together.
    #If not this does not work!! Be careful!
            
    # Next assemble all the stiffnesses into the correct places in the global stiffness matrix
    #1. split K_unit into 4 parts:
    #   F1 => u1 K_unit[:4,:4]
    #   F1 => u2 K_unit[:4,4:8]
    #   F2 => u1 K_unit[4:8,:4]
    #   F2 => u2 K_unit[4:,4:]
    #1. look at the connecting nodes [bot,top]
    #2. use that information to figure out where 
    for unit in range(0,n_units):
        
        unit_connections = units_raw[unit,1:3]
        bot_con = unit_connections[0]
        top_con = unit_connections[1]
        #Find index in nodes (sorted list of nodes. THis is the order they appear in the stiffness matrix.)
        #sorted in ascending order. top_con does not necessarily == top_ind but it can
        top_ind = int(list(np.where(nodes == unit_connections[1]))[0])
        bot_ind = int(list(np.where(nodes == unit_connections[0]))[0])
        
        k11 = K_unit[:nodal_dofs,:nodal_dofs,unit]
        k12 = K_unit[:nodal_dofs,nodal_dofs:,unit]
        k21 = K_unit[nodal_dofs:,:nodal_dofs,unit]
        k22 = K_unit[nodal_dofs:,nodal_dofs:,unit]
        
        top_loc = top_ind*nodal_dofs
        bot_loc = bot_ind*nodal_dofs
        #top to top (1:1)
        K_global[top_loc:top_loc+nodal_dofs,top_loc:top_loc+nodal_dofs] += k11
        
        #top to bottom (1:2)
        K_global[top_loc:top_loc+nodal_dofs,bot_loc:bot_loc+nodal_dofs] += k12
        #bottom to top (2:1)
        K_global[bot_loc:bot_loc+nodal_dofs,top_loc:top_loc+nodal_dofs] += k21
        #bottom to bottom (2:2)
        K_global[bot_loc:bot_loc+nodal_dofs,bot_loc:bot_loc+nodal_dofs] += k22
        a=1
    
    
    
    

    
    return K_global







def aAndbMatrixCalc(properties,nodal_dofs):
    
    #Calculates the a and b matrix for each unit.
    
    units_raw = properties["Raw"]    
    n_units = len(properties)-1
    
    nodes = np.unique(properties["Raw"][:,1:3])
    n_nodes = len(nodes)
    
    n_dof_unit = nodal_dofs * 2 # each unit is connected to two nodes so number of dofs at each node * 2
    n_dof_global = n_nodes*nodal_dofs #ux,uy,uz,rz,uyy,uzz for each node  
    tn,dn_range = util.neutralAxisCalc(properties)#neutral axis calc at each node. NOT per unit
#Calculate stiffnesses in every unit
    for unit in range(0,n_units):
        
        Ax = properties[str(unit)]["Area"]
        L = properties[str(unit)]["Length"]
        
        #temp
        Lx = L[0]
        Ly = L[1]
        Lz = L[2]
        Ay = Ax
        Az = Ax
        
        nu = properties[str(unit)]["Parameters"][1]
        
        
        nu_e = nu/(1-nu)
        
        
        #distance of unit from neutral axis.
        unit_dn = units_raw[unit,3:5]
        #Construct B matrix - 6 stresses, 8 dofs per unit
        b_mat = np.zeros([6,n_dof_unit]) 
        
        
        #part of b_mat for each bar in element
        bx = np.zeros([6,n_dof_unit]) 
        by = np.zeros([6,n_dof_unit]) 
        bz = np.zeros([6,n_dof_unit]) 
        
        #change this for altering dof connections with stress
        #TOp left referes to the TOP of the unit
        
        #distance from neutral axis
        #1 find distance of connection from neutral axis
        #needs to be done separately for each node.
        #d_lhs = unit_dn + 0.75 #unit_dn + connection point
        #d_tn = d_lhs - 0.5 #d_lhs - tn
        
        
        #Bx (axial bar)
        bx[0,0] = 1
        bx[0,3] = unit_dn[1] #
        
        bx[4,2] = 1
        
        bx[:,8:12] = bx[:,:4]*-1
        bx[0,11] = unit_dn[0]*-1
        
        #shear - xy - wx/2*(L-x?)
        #Here L = total area of all bars in unit? ()
        #x = distance from left hand side
        ax = np.copy(bx) 
        
        
        kap = 2/3 # I think I need to include this because im not accounting for the fact tau_xy = 3V/2A and how this interacts with the strain?
        #tau_xy = G.gamma_xy
        #3V/2A = G.gamma_xy
        #V = G.A.2/3.gamma_xy
        
        
        
        
        h_0 = dn_range[units_raw[unit,1]]#total 'height' of node. This is the total area of all bars on that node? (Or is it the total area covered (so dn_max - dn_min)?)
        h_1 = dn_range[units_raw[unit,1]]
        if tn[units_raw[unit,1]] - unit_dn[0] > 0:
            y_0 =  tn[units_raw[unit,1]] - unit_dn[0]-Ax/2
        else:
            y_0 =  tn[units_raw[unit,1]] - unit_dn[0]+Ax/2
        
        if tn[units_raw[unit,2]] - unit_dn[1] > 0:
            y_1 =  tn[units_raw[unit,2]] - unit_dn[1] - Ax/2
        else:
            y_1 =  tn[units_raw[unit,2]] - unit_dn[1] + Ax/2
        
        width = 1#out-of-plane width of beam (b)
        

        #find Ay for unit
        if y_0 <= 0:
            A_0 = y_0 + h_0/2
            y_0_bar = abs(y_0-A_0/2)
            
        else:
            A_0 = -y_0 + h_0/2
            y_0_bar = abs(y_0+A_0/2)
            
            
        if y_1 <= 0:
            A_1 = y_1 + h_1/2
            y_1_bar = abs(y_1-A_1/2)
            
        else:
            A_1 = -y_1 + h_1/2
            y_1_bar = abs(y_1+A_1/2)
        
        #Ay for unit
        Ay_0 = A_0*y_0_bar
        Ay_1 = A_1*y_1_bar
        
        
        #Ay @ mid of node. 
        A_0_mid = h_0/2
        A_1_mid = h_1/2
        
        y_0_mid = A_0_mid/2
        y_1_mid = A_1_mid/2
        
        Ay_0_mid = A_0_mid*y_0_mid
        Ay_1_mid = A_1_mid*y_1_mid
        

        Ay_Ay_0 = Ay_0/Ay_0_mid
        Ay_Ay_1 = Ay_1/Ay_1_mid
        
        

            #r_xz
        bx[3,3] = Ay_Ay_0*Lx*kap
        bx[3,11] = -Ay_Ay_1*Lx*kap
        
        #u_xy
        bx[3,1] = Ay_Ay_0*kap
        bx[3,9] = -Ay_Ay_1*kap
        #copy unit layout to a before *1/L
        #print(y_0,Q_0)
        
        bx *=1/Lx
        
        
        ax[3,3] = 1
        ax[3,11] = -1
        
        ax[3,1] = 1
        ax[3,9] = -1
        ax = ax * Ax
        
        
        #By (transverse bar)
        by[1,4] = 1
        by[1,12] = -1
        
        by[1,6] = unit_dn[1]
        by[1,14] = -unit_dn[0]
        
        #copy unit layout to a before *1/L
        ay = np.copy(by) * Ay
        
        by *=1/Ly
        
        #Bz (out of plane bar)
        bz[2,5] = 1
        bz[2,13] = -1
        
        bz[2,7] = unit_dn[1]
        bz[2,15] = -unit_dn[0]
        
        #copy unit layout to a before *1/L
        az = np.copy(bz) * Az
        
        bz *=1/Lz
        
        b_mat = bx + by + bz

        #A matrix construction
        a_mat = ax + ay + az
  
        a_mat = np.transpose(a_mat)

                
        
        
        
        #assign 1/L and Area to B and A respectively 
        #b_mat *= 1/L
        #a_mat *= A
        
        properties[str(unit)]["a_mat"] = a_mat
        properties[str(unit)]["b_mat"] = b_mat


    return properties





def aAndbMatrixCalc_axi(properties,nodal_dofs):
    
    #Calculates the a and b matrix for each unit.
    
    units_raw = properties["Raw"]    
    n_units = len(properties)-1
    
    nodes = np.unique(properties["Raw"][:,1:3])
    n_nodes = len(nodes)
    
    n_dof_unit = nodal_dofs * 2 # each unit is connected to two nodes so number of dofs at each node * 2
    n_dof_global = n_nodes*nodal_dofs #ux,uy,uz,rz,uyy,uzz for each node  
    tn,dn_range = util.neutralAxisCalc(properties)#neutral axis calc at each node. NOT per unit
#Calculate stiffnesses in every unit
    for unit in range(0,n_units):
        
        Ax = properties[str(unit)]["Area"]
        L = properties[str(unit)]["Length"]
        
        
        r_inner = L[2]
        #temp
        Lx = L[0]
        Ly = L[1]
        
        t_lhs = [0.5+units_raw[unit,3],0.5+units_raw[unit,4]]
        r_mid = r_inner+t_lhs[0]
        r_outer = r_mid + Ax/2
        
        r_inner = r_outer - 0.1
        Ax = np.pi*(r_outer**2 - r_inner**2)#top of cylinder (pi.(r_o^2-r_i^2))
        Ay = 2*np.pi*(r_outer)*Ax #side of cylinder (2.pi.r.h)
        Az = Ax
        
        A_tot = sum(units_raw[:,5])
        
        #print('=====')
        #print(unit)
        #print(Ax)
        #print(Ay)
        
        
        
        nu = properties[str(unit)]["Parameters"][1]
        
        
        nu_e = nu/(1-nu)
        
        
        #distance of unit from neutral axis.
        unit_dn = units_raw[unit,3:5]
        #Construct B matrix - 6 stresses, 8 dofs per unit
        b_mat = np.zeros([6,n_dof_unit]) 
        
        
        #part of b_mat for each bar in element
        bx = np.zeros([6,n_dof_unit]) 
        by = np.zeros([6,n_dof_unit]) 
        bz = np.zeros([6,n_dof_unit]) 
        
        #change this for altering dof connections with stress
        #TOp left referes to the TOP of the unit
        
        #B matrix:
        
        #strain_xx (row 0)
        bx[0,0] = 1
        bx[0,3] = unit_dn[1] #
        bx[0,8] = -1
        bx[0,11] = -unit_dn[0]
        
        #strain_yy (row 1)
        bx[1,4] = 1
        bx[1,12] = -1
        
        #strain_theta (row 2) - should be different for every unit
        bx[2,4] = 1
        bx[2,12] = 1
        
        #shear_theta_y (row 3)
        #bx[3,1] = 1
        #bx[3,9] = -1
        
        
        
        
        
        #copy unit layout to a before *1/L
        ax = np.copy(bx)
        
        
        bx[0,:] *= 1/Lx
        bx[1,:] *= 1/Ly
        bx[2,:] *= 1/(r_mid)
        bx[3,:] *= 1/Lx

        
        b_mat = bx

        #A matrix construction
        a_mat = ax 
  
        a_mat = np.transpose(a_mat)

        ax[:,0] *= Ax
        ax[:,1] *= Ay
        ax[:,2] *= Ax
        ax[:,3] *= Ax
        
        
        

        
        properties[str(unit)]["a_mat"] = a_mat
        properties[str(unit)]["b_mat"] = b_mat


    return properties