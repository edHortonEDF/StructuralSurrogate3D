# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 20:44:35 2025

@author: horto1e
"""

import numpy as np
import unitConfigSolver as unitSolver
import stiffnessMatrix as sm
import sys
    
def NewtonRaphson(x_0,y_target,du_app,properties):
    
    #Inputs:
    #1. Initial guess for x with n dimensions
    #2. Target for y with n dimensions
    #3. method relating x to y
    
    #x_0 = initial du_g
    #x_guess = du_g guess
    #K = constraint modified stiffness matrix
    #y_guess = dF
    #y_target = array of applied forces and constraints.
    
    #I For this script input and output need to have same dimensions (need square K for inversion)
    
    n = len(x_0)
    n_dofs_global = n
    
    
    
    tolerance = 1e-6
    pert = 1e-8
    NumericalJ = 1 #This only is needed for use with the dR_ddu nr_method. It does nothing with Kmod.
    nr_method = 'dR_ddu' #Kmod or dR_ddu.
    #speed rankings
    #1. dR_ddu w/ analytical Jacobian
    #2. Kmod 
    #3. dR_ddu with numerical Jacobian
    
    #m = len(y_target)
    
    #Initialisation
    x_guess = x_0
    y_guess = np.zeros([n])
    
    is_constrained = np.array([du != 'N' for du in du_app])
    
    #Guess 0
    #print('iteration 0') #Yes this should be x_guess,x_guess. x_old is not created yet.
    y_guess,K = unitSolver.duTodF(x_guess,x_guess,properties,du_app)
    
    y_adjusted,stiff_adjusted = constraintAdjustments(x_guess,y_guess,K,du_app)  
    y_target_adj = np.copy(y_target)
    y_target_adj[is_constrained] = 0

    
    #y_guess = y_target             
    #convergence?
    #residual = max(abs(y_target_adj-y_adjusted))
    residual = np.max(np.abs((y_target - y_guess)[~is_constrained]))

    #guess i
    if residual <= tolerance:
       conv = 1
    else:
       conv = 0
    
    #K = sm.generateStiffnessMatrix(properties)
    it_count = 0
    dR_ddu_fix = stiff_adjusted
    #K_elastic = K_mod
    while conv == 0:
        
        it_count += 1
        #print('iteration', it_count)
        if it_count > 10000:
            print('max iterations reached')
            sys.exit()
            
        x_old = x_guess
        y_old = y_guess
        
        
        
        
        
        #Adjust x and y for           
        y_adjusted,stiff_adjusted = constraintAdjustments(x_guess,y_guess,K,du_app)
        y_target_adj = np.copy(y_target)
        y_target_adj[is_constrained] = 0


        if nr_method == 'dR_ddu':
            x_inc = np.linalg.solve(dR_ddu_fix,y_target-y_adjusted) 
            if NumericalJ == 0 and it_count == 1:
                print('Analytical CTM')
            if NumericalJ == 1 and it_count == 1:
                print('Numerical CTM')   
        else:
            print('Newton Raphson Method Invalid. Check spelling')
        
        x_guess = x_inc + x_old
        
        y_guess,K = unitSolver.duTodF(x_guess,x_old,properties,du_app)

        
        y_adjusted,stiff_adjusted = constraintAdjustments(x_guess,y_guess,K,du_app)
        y_target_adj = np.copy(y_target)
        y_target_adj[is_constrained] = 0
        #R = dF_target - dF_calc(du)
        #residual = max(abs(y_target_adj-y_adjusted))
        residual = np.max(np.abs((y_target - y_guess)[~is_constrained]))
        #print('residual:', residual)
        
        
        if abs(residual) < tolerance:
            conv = 1
            print('converged in', it_count,' iterations' )
            #print('du: ',x_guess)
            
            
                
            
        else:
            #calculate consistent tangent modulus (dR/ddu)
            #for a small change is ddu, see how the residuals change.
            
            #For analytical method:
            #construct stiffness matrix out of C_ctm property of units and use that as dR_ddu.
            #still need to fix it for constraint
            #dR_ddu_ana = np.zeros([int(n_nodes*2),int(n_nodes*2)])
            #dR_ddu_ana = sm.generateStiffnessMatrixCTM(properties)
            y_old = y_guess
            #Numerical Jacobian - uses central difference to find updated residuals.
            if NumericalJ == 1:
                dR_ddu = np.zeros([len(K),len(K)])
                
                for dof in range(0,int(n_dofs_global)):

                    #backward du
                    x_back = x_old + x_inc 
                    x_back[dof] = x_back[dof] - pert
                    y_back,dummy = unitSolver.duTodF(x_back,x_old,properties,du_app)
                    R_back = y_target-y_back #backward residual 
                    
                    
                    #Forward du
                    x_forw = x_old + x_inc
                    x_forw[dof] = x_forw[dof] + pert
                    y_forw,dummy = unitSolver.duTodF(x_forw,x_old,properties,du_app)
                    R_forw = y_target - y_forw # forward residual
                    
                    #dR_ddu[i,j] = dR[i]/ddu[j]
                    for col in range(0,int(n_dofs_global)):
                            dR_ddu[dof,col] = -(R_forw[col]-R_back[col])/(x_forw[dof]-x_back[dof])
                            #print(dR_ddu[dof,col])
                #modify dR_ddu for constraint.
                #3 === modify stiffness matrix due to constraints ===
                dR_ddu_fix = np.copy(dR_ddu)
            
                #I.e. apply boundary conditions - i.e. modify for constraint
                dummy = np.zeros([len(du_app)])
                dummy,dR_ddu_fix = constraintAdjustments(x_guess,dummy,dR_ddu,du_app)  
            else:
                #ANALYTICAL CONSTRAINT
                dR_ddu_ana = np.zeros([int(n_dofs_global),int(n_dofs_global)])
                dR_ddu_ana = sm.generateStiffnessMatrix(properties)
                #modify dR_ddu for constraint.
                #3 === modify stiffness matrix due to constraints ===
                dR_ddu_fix = np.copy(dR_ddu_ana)
                
                #I.e. apply boundary conditions - i.e. modify for constraint
                dummy = np.zeros([len(du_app)])
                dummy,dR_ddu_fix = constraintAdjustments(x_guess,dummy,dR_ddu_ana,du_app) 
    return x_guess,properties



#inputs 
#x_initial = [0,0] # du_g_0 (so 0)
#y_target = [100,512] # -dF_app. In other words dF_rf_app

#x_conv = NewtonRaphson(x_initial,y_target)



#This method adjusts the force and stiffness array in the correct manner.
def constraintAdjustments(x_guess,y_guess,stiff,du_app):
    
    is_constrained = np.array([du != 'N' for du in du_app])

    
    #Adjust force vector
    y_adjusted = np.copy(y_guess)
    y_adjusted[is_constrained] = 0
    

    
    #Adjust stiffness matrix
    stiff_adjusted = np.copy(stiff)
    stiff_adjusted[is_constrained, :] = 0
    stiff_adjusted[:, is_constrained] = 0
    for i in range(len(stiff)):
        if is_constrained[i]:
            stiff_adjusted[i, i] = 1
    for i in range(len(stiff)):
        if stiff_adjusted[i,i] == 0:
            stiff_adjusted[i,i] = 1
    
    return y_adjusted, stiff_adjusted











