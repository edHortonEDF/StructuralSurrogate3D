
**Single Element Formulation**

Consider a single beam element with 2 nodes. This element has 8 degrees of freedom at each node which allow for a 6 dimensional stress tensor to be calculated for any given applied load or displacement. 

The three beams correspond to an axial beam with 4 degrees of freedom ($u_{xx},u_{xy},u_{xz},\theta_{xz}$), a transverse beam ($u_{yy},\theta_{yz}$), and an out-of-plane beam ($u_{zz},\theta_{zx}$). The forces applied to an element, f_{e},can be related to the resulting deflections and rotations, u_e, via the element stiffness matrix, $k_e$, as follows:

$f_e=k_e u_e$


The element stiffness matrix is constructed using a set of three governing equations:
	The degrees of freedom related to the element, $u_e$, can be related to the deformation (represented as strain), $\epsilon$ via:
	
$\epsilon=Bu_e$

Where B is a matrix that defines the relationship between $u_e$ and $\epsilon$. 

The strain is related to the internal forces (represented as stress) via:

$\sigma=C\epsilon$

Where $\sigma_{ij}$ is the stress and $C_{ijkl}$ is the 6x6 material model. 

The stress is related to the applied forces via Equation 4:

$f_e=A^T \sigma$	

Where $A^T$ is a matrix that defines the relationship between $f_e$ and $\sigma$. In this case $A^T=B^T$.

With these governing equations it is possible link an applied force to a relative change in the length of the element, assuming the problem is fully constrained (i.e. one node is pinned).


**General Solver Functionality**

Firstly, the system must be fully constrained. The simplest manner for achieve this is by forcing all the degrees of freedom at once node to be 0, i.e. pinned in place. If the system is not fully constrained the problem is unsolvable and the surrogate model will fail. 
An initial guess for the degrees of freedom, u, is chosen (usually 0 in all positions). At this point each element is run through the material model with the associated degrees of freedom to find the strain and stress.
The stress is converted into the resulting forces in the unit via:

$df_e  =Ad\sigma_e$


This gives the forces associated with the unit positions in the correct location in the context of all the degrees of freedom in the global stiffness matrix ($df_e^{global}$).
Through this method the total combination of resultant forces due to the degrees of freedom of all units can then be found by simply adding  $df_e^{global}$  from every unit together.

$F = \sum_{i=1}^{N}df_{e,i}^{global}$

This is then compared to the applied forces vector, ignoring any rows which contain Dirichlet boundary conditions. If the residual is greater than a set tolerance, then the initial guess for $u$ has not reached equilibrium and a new guess for u is found using:

$u=K^{-1}.F$

The above process is repeated iteratively using a Newton-Raphson convergence algorithm until a u that satisfies equilibrium is found.


**Material Model Functionality**

A key part of the FEA solver is the calculation of stress given a set of degrees of freedom. This occurs in the material model as follows:

The total strain increment, $d\epsilon^t$, is calculated via:

$d\epsilon^t  =Bdu_e$

The surrogate considers thermal loads so the total mechanical strain, $d\epsilon^{t,m}$,  must be separated from the thermal strain, $d\epsilon^{\theta}$, via:

$d\epsilon^{t,m}=d\epsilon^{t}-d\epsilon^{\theta}$

Where  d\epsilon^{\theta} is calculated as follows:

$d\epsilon^{\theta}=\alpha.dT.\delta$

Where  $\alpha$ is the thermal expansion coefficient, $dT$ is the temperature difference in the increment and $\delta$ is the Kronecker delta.
The trial stress, σ^tr, can then be found via:

$d\sigma^{tr}=Cd{\epsilon}^(t,m)$

as

$\sigma^{tr}_t=\sigma^tr_{t-1} + d\sigma^{tr}$

This trial stress is checked against the yield function, $f^y$, if it is less than or equal to zero then the behaviour is elastic, otherwise there is plasticity. In the case of linear elasticity $\sigma^tr_t=\sigma_t$ at all times.
The reaction force increment relative to the degrees of freedom in the element can then be found via.

$f_e=A^T.d\sigma^tr$

Finally, the element reaction force is added to the global reaction force matrix for comparison to the applied forces.


**Chaboche Material Model**

A Lemaitre-Chaboche material model has been created. It follows the same formulation as the ABAQUS implementation, allowing for simple validation between the methods. The general form is given below.
To replicate the complex behaviour seen in 316H such as cyclic hardening and ratchetting, a Lemaitre-Chaboche combined hardening material model was chosen. The built-in ABAQUS implementation of the selected model was used so details can be found in the documentation [1]. A brief outline of the approach is given below.
Initially it is assumed that each loading increment results in elastic behaviour only, i.e. σ_ij^tr=C_ijkl ϵ_ij. To test if this assumption is valid, the trial stress (σ_ij^tr) is used to check if the yield surface condition is met, Equation 21:
√(3/2 (S_ij-X_ij ):(S_ij-X_ij )  )-σ^0≤0	Equation 21

Where  S_ij is the deviatoric stress tensor, X_ij is the backstress tensor that defines the kinematic hardening behaviour and σ^0 is the isotropic hardening component.
The deviatoric stress tensor is calculated via Equation 22:
S_ij=σ_ij-1/3 tr(σ_ii )	Equation 22

The evolution of X_ij is defined through three backstress terms by Equation 23 and Equation 24.
dX_ij= ∑_(k=1)^3▒〖dX_ij^k 〗	Equation 23


dX_ij^k=C_k dϵ ̅^pl-γ_k X_ij^k |dϵ ̅^pl |	Equation 24

C_k and γ_k are material parameters, and dϵ ̅^pl is the increment of equivalent plastic strain. The isotropic hardening component, σ^0 is found via Equation 25:
σ^0=σ_ys^0+Q(1-exp⁡(-bϵ ̅^pl ) )	Equation 25

where σ_ys^0 is the yield stress, and Q and b are material parameters. A goal seek method, such as Newton-Raphson, is employed to find the value of dϵ ̅^pl that ensures Equation 21 is equal to zero.

The combined hardening material model used the linear elastic parameters given in Error! Reference source not found. combined with the combined hardening parameters given in Table 2 1. The isotropic parameters (Q and b) were taken from fits to a cyclic strain-controlled test (Δϵ=0.9%,Rϵ=-2) on 316H ELG at 550℃ termed ALM53. The methodology for the fitting can be found in [2]. 
The kinematic hardening parameters (C_i and γ_i) were also taken from [2]. It must be noted that the strain is assumed to be in percent not absolute. The fit between the Chaboche model and ALM53 is given in Error! Reference source not found.. This work is not focused on fitting methods so, for the purposes of this report the chosen parameters represent a reasonable approximation of the behaviour of 316H at 550℃.

Table 2 1: Lemaitre-Chaboche combined hardening model parameters. 
Property	Value	Unit
E	200000	MPa
ν	0.3	
σ_ys	50	MPa
C	70000	MPa
D	500	
Q	400	MPa
b	5	

