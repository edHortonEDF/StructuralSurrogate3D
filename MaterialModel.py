# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 09:09:51 2025

@author: horto1e
"""


import numpy as np
import random as rd
import sys

#convert stress or strain vector to matrix
def vec2mat(s_vec):
    s_mat = np.zeros([3,3])
    for i in range(0,3):
        s_mat[i,i] = s_vec[i]
        
    s_mat[0,1] = s_vec[3]
    s_mat[0,2] = s_vec[4]
    s_mat[1,2] = s_vec[5]
    
    s_mat[1,0] = s_vec[3]
    s_mat[2,0] = s_vec[4]
    s_mat[2,1] = s_vec[5]
    return s_mat


#convert stress or strain matrix to vecotr.
def mat2vec(s_mat):
    s_vec = np.zeros([6])
    s_vec[0] = s_mat[0,0]
    s_vec[1] = s_mat[1,1]
    s_vec[2] = s_mat[2,2]
    s_vec[3] = s_mat[0,1]
    s_vec[4] = s_mat[0,2]
    s_vec[5] = s_mat[1,2]
    
    return s_vec



def LinearElastic(unit_props,du_g, node_list):

    
    #

    
    a_mat = unit_props["a_mat"]
    b_mat = unit_props["b_mat"]
    D_e = unit_props["D_e"]
    #Material parameters
    E = unit_props["Parameters"][0]
    nu = unit_props["Parameters"][1]
    #sig_y = unit_props["Parameters"][1]
    #H = unit_props["Parameters"][2]
    
    #State Variables
    eto_bts= unit_props["StateVariables"]["eto"]
    eel_bts= unit_props["StateVariables"]["eel"]
    p_bts  = unit_props["StateVariables"]["ep"]
    stress_bts = unit_props["StateVariables"]["stress"]
    alpha_bts = unit_props["StateVariables"]["alpha"]
    
    T_bts = unit_props["StateVariables"]["T"]
    dT = unit_props["StateVariables"]["dT"]
    #deto
    #deel
    #dp
    

    bot_con = unit_props["Connections"][0]
    top_con = unit_props["Connections"][1]
    #Find index in nodes (sorted list of nodes. THis is the order they appear in the stiffness matrix.)
    #sorted in ascending order. top_con does not necessarily == top_ind but it can
    top_ind = int(list(np.where(node_list == top_con))[0])
    bot_ind = int(list(np.where(node_list == bot_con))[0])
    
    
    nodal_dofs = int(np.size(a_mat,axis=0)/2)
    
    top_loc = top_ind*nodal_dofs
    bot_loc = bot_ind*nodal_dofs
    
    
    du_unit = np.zeros([nodal_dofs*2]) #x2 because du_units has top and bottom of bar in

    du_unit[:nodal_dofs] = du_g[top_loc:top_loc+nodal_dofs]
    du_unit[nodal_dofs:] = du_g[bot_loc:bot_loc+nodal_dofs]
    
    
    #add shear components

    
    #dstrain - e_t  = e_e + e_p + e_th
    deto = b_mat @ du_unit
    #we are interested in CONSTRAINED thermal behaviour so thermal strain is negative.
    deth = np.zeros([6])
    deth[:3] = dT*16e-6 #16e-6 is the thermal expansion. Add this to the material properties please.
    #deth[3] = -deth[0]*1

    deel = deto - deth
    #deel[3] = -deth[0]*0.1
    
    dstress_tr = D_e @ deel
    
    stress_tr = stress_bts + dstress_tr

    #All linear elastic so trial = actual
    dstress = dstress_tr
    
    unit_props["C_ijkl"] = D_e
    #print('total strain:', deto)
    #print('thermal strain:', deth)
    #print('elastic strain:',deel)
    #print('=====')
    

    #internal forces - element form
    df_el = a_mat @ dstress
    print(df_el)
    #move internal forces to global nodes
    df_global = np.zeros([len(du_g)])
    
    df_global[top_loc:top_loc+nodal_dofs] = df_el[:nodal_dofs]
    df_global[bot_loc:bot_loc+nodal_dofs] = df_el[nodal_dofs:]
    

    #COnsistent tangent modulus
    C_ctm = D_e
    
    
    stress = stress_tr
    dstress = stress - stress_bts
    #deel = deto-deth
    dalpha = np.zeros([6])
    #update state variables 
    #deto = deel - deth
    
    unit_props["StateVariables"]["deto"] = deto
    unit_props["StateVariables"]["deel"] = deel
    unit_props["StateVariables"]["dep"] = 0
    unit_props["StateVariables"]["dstress"] = dstress
    unit_props["StateVariables"]["dalpha"] = dalpha
    unit_props["StateVariables"]["deth"] = deth
    
    unit_props["C_ijkl"] = C_ctm
    return df_global,unit_props
    




# ========== Chaboche hardening

def Chaboche(unit_props,du_g, node_list):
    n_nodes = int(len(du_g)/2)
    #du_unit = np.zeros([8])
    #
    Area = unit_props["Area"]
    Length = unit_props["Length"]
    
    a_mat = unit_props["a_mat"]
    b_mat = unit_props["b_mat"]
    D_e = unit_props["D_e"]
    #Material parameters
    E = unit_props["Parameters"][0]
    nu = unit_props["Parameters"][1]
    sig_y = unit_props["Parameters"][2]
    C = unit_props["Parameters"][3]
    D = unit_props["Parameters"][4]
    Q = unit_props["Parameters"][5]
    b = unit_props["Parameters"][6]
    H = C
    #sig_y = unit_props["Parameters"][1]
    #H = unit_props["Parameters"][2]
    
    

    
    #State Variables
    eto_bts= unit_props["StateVariables"]["eto"]
    eel_bts= unit_props["StateVariables"]["eel"]
    p_bts  = unit_props["StateVariables"]["ep"]
    stress_bts = unit_props["StateVariables"]["stress"]
    alpha_bts = unit_props["StateVariables"]["alpha"]
    #deto
    T_bts = unit_props["StateVariables"]["T"]
    dT = unit_props["StateVariables"]["dT"]
    

    bot_con = unit_props["Connections"][0]
    top_con = unit_props["Connections"][1]
    #Find index in nodes (sorted list of nodes. THis is the order they appear in the stiffness matrix.)
    #sorted in ascending order. top_con does not necessarily == top_ind but it can
    top_ind = int(list(np.where(node_list == top_con))[0])
    bot_ind = int(list(np.where(node_list == bot_con))[0])
    
    
    nodal_dofs = int(np.size(a_mat,axis=0)/2)
    
    top_loc = top_ind*nodal_dofs
    bot_loc = bot_ind*nodal_dofs
    
    
    du_unit = np.zeros([nodal_dofs*2]) #x2 because du_units has top and bottom of bar in

    du_unit[:nodal_dofs] = du_g[top_loc:top_loc+nodal_dofs]
    du_unit[nodal_dofs:] = du_g[bot_loc:bot_loc+nodal_dofs]
    
    EG = E/(2*(1+nu))

    
    #equivalent plastic strain - before time step
    p_bts_contract = np.sum(p_bts**2)
    p_bts_eq = np.sqrt(2/3*p_bts_contract)
    
    #strain = (du_top-du_bot)/Length
    #need displacement AND rotation.
    deto = b_mat @ du_unit
    
    
    #we are interested in CONSTRAINED thermal behaviour so thermal strain is negative.
    deth = np.zeros([6])
    deth[:3] = dT*16e-6
    #deto = deto + deth
    
    deto_tr = np.copy(deto) - deth
    

    
    M1 = np.zeros([6,6])
    for i in range(0,3):
        M1[i,i] = 1
    for i in range(3,6):
        M1[i,i] = 2
    #deto_mat = vec2mat(deto)
    
    #calculate predictor stress
    dstress_tr = D_e @ deto_tr
    
    stress_tr = stress_bts + dstress_tr
    
    s_hyd_tr = (stress_tr[0] + stress_tr[1] + stress_tr[2])/3
    
    
    #deviatoric
    s_tr = np.copy(stress_tr)
    for i in range(0,3):
        s_tr[i] = s_tr[i] - s_hyd_tr
    
    #stress - alpha
    eta_tr = s_tr - alpha_bts
    
    
    eta_tr_contract = np.sum(eta_tr**2)
    
    stress_tr_eq = np.sqrt(3/2*eta_tr_contract)
    #stress_tr_eq = eta_tr[0]
    
    
    #Isotropic
    
    sig_0_tr = sig_y + Q*(1-np.exp(-b*p_bts_eq))
    
    
    f_tr = stress_tr_eq - sig_0_tr
    
    if f_tr <= 0:
        deel = np.copy(deto_tr) - deth
        
        dstress = D_e @ deel
        
        deto = deel
        
        dep = np.zeros([6])
        dalpha = np.zeros([6])
        C_ctm = D_e
    else:       
        
        flow = eta_tr/stress_tr_eq
        
        dep_eq = 0.0#f_tr/(2*EG+C+b*Q)
        
        tol = 1e-8
        conv = 0
        it = 0
        while conv == 0:
            it += 1
            if it > 1000:
                print('Material model failed to converge')
                print(dep_eq)
                sys.exit()
            
            
            # Trial updates of internal variables
            dep = 1.5 * dep_eq * flow  # plastic strain increment
            #dalpha = dep_eq * (C * flow - D * alpha_bts)
            #alpha_ets = alpha_bts + dalpha  # backstress
            p_eq = p_bts_eq + dep_eq  # accumulated plastic strain
            
            K = Q * (1 - np.exp(-b * p_eq))
            
            sig_0_ets = sig_y + K  # isotropic hardening

            if D == 0:
                dalpha = flow * H * dep_eq - D*alpha_bts*dep_eq
                alpha_ets = alpha_bts + dalpha 
            else:
                alpha_ets = flow*(C/D*(1-np.exp(-D*dep_eq)))+alpha_bts*np.exp(-D*dep_eq)
                dalpha = alpha_ets -alpha_bts
            
            #deviatoric stress -
            deel = deto - dep - deth
            
            #deel[1] = (nu*dep[2]-nu*deel[0])/(1-nu) # plane strain correction
            dstress = D_e @ deel
            
            
            stress = stress_bts + dstress
            
            s_hyd = (stress[0] + stress[1] + stress[2])/3 
            
            s_ets = np.copy(stress)
            for i in range(0,3):
                s_ets[i] = s_ets[i] - s_hyd
            
            #s = s_tr - 2*EG*dep
            
            eta = s_ets-alpha_ets
            
            eta_contract = np.sum(eta**2)
            
            stress_eq = np.sqrt(3/2*eta_contract)
            #stress_eq = eta[0]
            
 
            f = stress_eq - sig_0_ets
            
            if abs(f) <= tol:
                 conv = 1
                 #print('converged in ', it, ' iterations')
                 #find full stress tensor

                 deel = deto - dep - deth
                 #deel[1] = (nu*dep[2]-nu*deel[0])/(1-nu)
                   
                 deto = dep + deel
                 
                 dstress = D_e @ deel
                 
                 
                 
            else:
                numerical = 1
                if numerical == 0:
                    dY_dep_eq = b*Q*np.exp(-b*p_eq)
                    h = C-np.sqrt(3/2)*np.dot(flow,D*alpha_bts)
                    df_dep_eq = -2*EG*1.5+h+dY_dep_eq
                    
                    dep_eq -= f/df_dep_eq
                else:
                #Numerical J for material model... maybe not the best
                    pert = 1e-10
                    p_vals = np.zeros([2])
                    f_vals = np.zeros([2])
                    
                    for i in range(0,2):
                        if i == 0:
                            if dep_eq < 1e-10:
                                dep_eq_pert = dep_eq
                            else:
                                dep_eq_pert = dep_eq - pert
                        elif i == 1:
                            dep_eq_pert = dep_eq + pert
                        
                        if D == 0:
                            dalpha = flow * H * dep_eq_pert - D*alpha_bts*dep_eq_pert
                            alpha_ets = alpha_bts + dalpha 
                        else:
                            alpha_ets = flow*(C/D*(1-np.exp(-D*dep_eq_pert)))+alpha_bts*np.exp(-D*dep_eq_pert)
                            dalpha = alpha_ets -alpha_bts
                        
                        sig_0 = sig_y + Q*(1-np.exp(-b*(p_bts_eq+dep_eq_pert)))
                        
                        dep_pert = 3/2 * flow * dep_eq_pert
                    
                    
                        deel = deto - dep_pert - deth
                    

             
                        dstress = D_e @ deel
                        
                        
                        stress = stress_bts + dstress
                        
                        s_hyd = (stress[0] + stress[1] + stress[2])/3 
                        
                        s = np.copy(stress)
                        for j in range(0,3):
                            s[j] = s[j] - s_hyd
                        
                        eta = s-alpha_ets
                        
                        eta_contract = np.sum(eta**2)
                        
                        stress_eq = np.sqrt(3/2*eta_contract)
                        #stress_eq = eta[0]
                    
                        f_vals[i] = stress_eq - sig_0
                        p_vals[i] = dep_eq_pert
                    
                    
                    
                    df_dep_eq = (f_vals[1]-f_vals[0])/(p_vals[1]-p_vals[0])
                    
                    d_g = -1/df_dep_eq*f
                    dep_eq = dep_eq + d_g
                
       
        #Consistent tangent modulus
        
        norm_eta  = np.sqrt(3/2*np.dot(eta,eta))
        dflow_dsig = (np.identity(6) - np.outer(flow, flow)) / norm_eta


        dK_dp = b*Q*np.exp(-b*p_bts_eq+dep_eq)
        
        dalpha_dp = C*np.exp(-D*dep_eq)*np.identity(6)
        
        H_prime = np.dot(flow,np.dot(D_e,flow)) + dK_dp + np.dot(flow,np.dot(dalpha_dp,flow))
        
        C_ctm = D_e - np.outer(np.dot(D_e, flow), np.dot(flow, D_e)) / H_prime
        #project = np.sqrt(3/2)*np.outer(flow,flow)


        #dalpha_dp = C*np.exp(-D*dep_eq) - D*alpha_bts*np.exp(-D*dep_eq)        

        
        
        #h = C-np.sqrt(3/2)*np.transpose(flow)*M1*D*alpha_ets#+dY_dp
        
        
        #C_ctm = D_e - 6*EG**2*((flow*np.transpose(flow))/(3*EG+h))
        
        #C_ctm = D_e - 6*(EG)**2/(3*EG+dalpha_dp+dY_dp)*project

    
    #internal forces - element form
    df_el = a_mat @ dstress

    #move internal forces to global nodes
    df_global = np.zeros([len(du_g)])
    
    df_global[top_loc:top_loc+nodal_dofs] = df_el[:nodal_dofs]
    df_global[bot_loc:bot_loc+nodal_dofs] = df_el[nodal_dofs:]
    


    
    

    #update state variables 
    unit_props["StateVariables"]["deto"] = deto 
    unit_props["StateVariables"]["deel"] = deel
    unit_props["StateVariables"]["dep"] = dep
    unit_props["StateVariables"]["dstress"] = dstress
    unit_props["StateVariables"]["dalpha"] = dalpha
    unit_props["StateVariables"]["deth"] = deth
    
    unit_props["C_ijkl"] = C_ctm
    return df_global,unit_props
    

 
