# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:02:18 2025

@author: horto1e
"""


import numpy as np
import pandas as pd
import math
import stiffnessMatrix as sm
import sys


def Convert2Float(value):
        try:
            value = float(value)
        except ValueError:
            pass
    
        return value
    
# open the properties file
def openLoadingFile(properties,path,load_fname,constraint_fname,temp_fname):
    
    units = properties["Raw"]
    
    n_units = len(units)
    
    load = pd.read_csv(path+load_fname)
    load = pd.DataFrame.to_numpy(load)
    
    constraint = pd.read_csv(path+constraint_fname)
    constraint = pd.DataFrame.to_numpy(constraint)
    
    temp = pd.read_csv(path+temp_fname)
    temp = pd.DataFrame.to_numpy(temp) 
    
    #find number of nodes
    nodes = np.unique(units[:,1:3])
    n_nodes = len(nodes)
    
    #n_t_steps = np.size(load,0)
    
    #load_only = load[:,1:np.size(load,1)]
    
    #n_dofs = np.size(load_only,1)
    #Initialize dofs
    #for node in range(0,n_nodes):
    #    F_app = np.zeros([n_t_steps,n_nodes])
    #    M_app = np.zeros([n_t_steps,n_nodes])
    #    u_app = np.zeros([n_t_steps,n_nodes])
    #    r_app = np.zeros([n_t_steps,n_nodes])
    
    F_app = load[:,1:np.size(load,1)]

    u_app = constraint[:,1:np.size(load,1)]

    T_app = temp[:,1:]
    
    load_dict = {
        "t":load[:,0],
        "F_app":F_app,
        "u_app":u_app,
        "T_app":T_app
        }
    
    return load_dict



# open the properties file
def openPropsFile(path,fname,nodal_dofs,model_type,constraint_type):
    units = pd.read_csv(path+fname)
    units = pd.DataFrame.to_numpy(units)
    
    
    
    n_units = np.size(units,0)
    #create dictionary entry for every unit
    properties = {}
    for unit in range(0,n_units):
        # 1. construct elastic stiffness matrix
        
        E = units[unit,10]
        nu =units[unit,11]
        
        
        #type of simulation - 3D,PS,PE
        
        #model_type = "1D"
        
        
        D_e = np.zeros([6,6])
        
        if model_type == "3D":
            if constraint_type == 'PE':
                #3D elasticity
                lam = E/((1+nu)*(1-2*nu))
                
                
                D_e[:3,:3] = nu
                for i in range(0,3):
                    D_e[i,i] = 1 - nu
                for i in range(3,6):
                    D_e[i,i] = (1-2*nu)/2
                    
                D_e = D_e*lam
            elif constraint_type == 'Axi':
                #3D elasticity
                lam = E/((1+nu)*(1-2*nu))
                
                
                D_e[:2,:2] = nu
                for i in range(0,2):
                    D_e[i,i] = 1 - nu
                    
                D_e[2,2] = 1-nu
                D_e[3,3] = (1-2*nu)/2
                    
                D_e = D_e*lam
                
                
            else:
                print("Plane stress no yet implemented")
                sys.exit()
        elif model_type == "1D":
            D_e[0,0] = E
        else:
            print('Invalid model type (1D or 3D)')
            sys.exit()
            

        
        state_variables = {
            "eto" : np.zeros([6]),
            "deto" : np.zeros([6]),
            "eel" : np.zeros([6]),
            "deel" : np.zeros([6]),
            "ep" : np.zeros([6]),
            "dep" : np.zeros([6]),
            "eth": np.zeros([6]),
            "deth": np.zeros([6]),
            "stress" : np.zeros([6]),
            "dstress" : np.zeros([6]),
            "alpha" : np.zeros([6]),
            "dalpha" : np.zeros([6]),
            "T": 0,
            "dT": 0
            
            }
        
        
        #2. construct unit dictionary entry
        unit_dict = {"Name": units[unit,0],
                    "Connections": units[unit,1:3],
                    "Position": units[unit,3:5],
                    "Area": units[unit,5],
                    "Length": units[unit,6:9], #change to have 3 lengths associated with it
                    "C_ijkl":  D_e,
                    "C_f": 0,
                    "D_e": D_e,
                    "Material Model": units[unit,9],
                    "Parameters": units[unit,10:np.size(units,1)],
                    "StateVariables": state_variables,    #This number is defined by the material model. 
                    "a_mat": 0,
                    "b_mat": 0,
                         }
        
        properties.update({
            str(unit):unit_dict
            })
        properties.update({
            "Raw":units
            })
        
        if constraint_type == "PE":
            properties = sm.aAndbMatrixCalc(properties,nodal_dofs)
        elif constraint_type == "Axi":
            properties = sm.aAndbMatrixCalc_axi(properties,nodal_dofs)
        else:
            print("Plane stress not here yet, sorry :/")
    return properties



def misescalc(s_vec):
    #s_vec - full 3D stress tensor
    s_hyd = (s_vec[0] + s_vec[1] + s_vec[2])/3
    
    
    #deviatoric
    s = np.copy(s_vec)
    for i in range(0,3):
        s[i] = s[i] - s_hyd
    
    #stress - alpha
    eta = s
    
    
    eta_contract = np.sum(eta**2)
    
    stress_eq = np.sqrt(3/2*eta_contract)

    return stress_eq,s_hyd



def neutralAxisCalc(properties):
    #ony gives neutral axis in the x direction currently.
    props = properties["Raw"]
    nodes = np.unique(props[:,1:3])
    n_units = np.size(props,0)
    dn_range = np.zeros([max(nodes)+1])
    
    tn = np.zeros([max(nodes)+1])
    
    #Always use props. the properties dictionary may not be fully filled out when this is run
    k=0
    for node in nodes:
        node_units = []
        top = 0
        bot = 0
        for unit in range(0,n_units):
            if props[unit,1] == node or props[unit,2] == node:
                node_units.append(unit)
                
                
        dn_list = []
        A_list = []
        for unit in node_units:
            #find most negative
            if props[unit,1] == node:
                dn_list.append(props[unit,3])
                A_list.append(props[unit,5])
            if props[unit,2] == node:  
                dn_list.append(props[unit,4])
                A_list.append(props[unit,5])
            
        dn_min = np.min(dn_list)
        dn_max = np.max(dn_list)
        min_ind = np.argmin(dn_list)
        max_ind = np.argmax(dn_list)
        A_min = A_list[min_ind]
        A_max = A_list[max_ind]
        dn_range[node] = (dn_max + A_max/2) - (dn_min - A_min/2) 
            
        for unit in node_units:
            A = props[unit,5]
            L = props[unit,6]
            E = props[unit,10]
            if props[unit,1] == node:
                dn = props[unit,3]
            if props[unit,2] == node:  
                dn = props[unit,4]
                
                
            t_LHS = dn# - dn_min + A/2
            bot += E*A/L
            
            #find the most negatibe dn. 
            #t_LHS = current_dn - neg_dn + A/2
            
            top += E*A/L*t_LHS
        
            
        tn[node] = top/bot
        #dn_range[node] = dn_max-dn_min
        k += 1
        
        
    return tn,dn_range # tn is relavtive to coordinates of node (so y = dn -A/2- tn)
    
    
    
    
    
    
    
    
    
    
    
    
    
    