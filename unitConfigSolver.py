# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 08:26:14 2025

@author: horto1e
"""
import numpy as np
import MaterialModel as matmod
import stiffnessMatrix as sm
import sys

def duTodF(du,du_old,properties,du_app):
    #Relation between x (du) and y (dF_rf)
    #1.input: du_guess
    #2. run du_g through material model to get C_ep for each unit
    #3. assemble global stiffness matrix, K
    #4. calculate dF_rf = K.du_g
    #5. output: dF_rf,K
    
    #May be able to use this with material model as well:
    #du = plastic strain increment (always >0) and y = yield function (or (sig-du)-sig_0)
    
    #This is where you add your relation between the inputs (du) and outputs (y)
    #Don't think there is any limitations as long as it is a valid relation.
    
    n_units = len(properties)-1
    
    dF = np.zeros([len(du)])
    
    #This needs to be adjusted to give the gradient at the point not dsig/deto.
    #recommend we start with numerical perturbation (+-1e-8)
    node_list = np.unique(properties["Raw"][:,1:3])
    for unit in range(0,n_units):
        #print('Unit:',unit)
        unit_props = properties[str(unit)]
        material_model = unit_props["Material Model"]
        #Select correect material model
        if material_model == 'LinearElastic':
            df_unit,unit_props = matmod.LinearElastic(unit_props, du, node_list)
        elif material_model == 'LinearK':
            df_unit,unit_props = matmod.LinearKH(unit_props, du, node_list)
        elif material_model == 'Chaboche':
            df_unit,unit_props = matmod.Chaboche(unit_props, du, node_list)
        else:
            print('No such material model. Options are: LinearElastic, LinearK, LinearK_LinearISO, Voce, Chaboche')
            sys.exit()
        #Options: LinearElastic,LinearK,LinearK_LinearISO,Voce,Combined
        #C_ep,unit_props = matmod.Voce(unit_props, du)
        #C_ep,unit_props  = matmod.LinearKH(unit_props, du)
        #C_ep,unit_props  = matmod.Combined(unit_props, du)
        #C_ep,unit_props  =  matmod.Combined(unit_props, du)

        dF += df_unit
        
    
    
    stiff = sm.generateStiffnessMatrix(properties)
    

    

        
    #update state variables in properties.
    
    
    #Outputs
    #y - function output (so dF_rf in surrogate)
    #K - gradient matridu at the point we are considering.
    #In the surrogate model case, this is the stiffness matridu
    return dF,stiff#,stiff