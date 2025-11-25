3D Surrogate model:
The surrogate model used in this work is designed to replicate the through-section mechanical behaviour captured by finite element analysis along a stress classification line. To achieve this replication of behaviour the local area around the SCL can be represented via a series of nodes connected by sets of ‘units’ that follow beam theory. Each unit captures the behaviour over a certain section of the line profile of surrounding material. By altering the stiffnesses of each unit in the surrogate model, it can be made to replicate the behaviour of an FEA model in areas containing stress concentrations that are induced by the geometry, such as a weld toe or notch. The surrogate has significant computational advantages over full FEA due to the fewer calculation points, needed to give an acceptable approximation of the structural response. 
The mathematical description of the ROM is given below:
The ROM functions as a series of custom 3D beam elements connected to nodes which enforce each element to share its degrees of freedom with every other element connected to that node. This means that the behaviour of each node can be changed by altering the stiffnesses of different elements along the section. Each element consists of three separate beams that are connected via the material model; an axial beam with 4 degrees of freedom (u_xx,u_xy,u_xz,θ_xz), a transverse beam (u_yy,θ_yz), and an out-of-plane beam (u_zz,θ_zx). This total of 8 degrees of freedom allows for both normal and shear stresses to be captured in the model. Figure X shows a diagram of the three beams that make up an element with their respective degrees of freedom.

The forces applied to an element,  f,can be related to the resulting deflections and rotations, u, via the element stiffness matrix, k, as follows:
f=ku
Every unit follows beam theory so to construct k, a set of three governing equations can be employed:

	The displacement, u, can be related to the deformation (represented by strain), ϵ, experienced by the material via Eq. X
ϵ_ij=〖du〗_j/(dL_i )
	Where L is the length of the beam. There is 6 strain components that are related to the 8 degrees of freedom so a 6x8 matrix termed the B matrix is created:
ϵ=Bu

	The strain can be related to the internal forces (represented at stress) via Eq. Y
σ_ij=C_ijkl ϵ_ij
Where σ_ij is the stress and C_ijkl is the 6x6 material model. 


	The stress can be related to the applied forces via Eq. Z:
f_ij=A_i.σ_ij
Where A_i is the area of the beam. There is 8 applied forces (same as degrees of freedom) related to 6 stress components so the A matrix relates these components via a 8x6 matrix.
f=A^T σ

Eq. X, Y,Z and A can be combined to give:
f=A^T.C.Bu
For a full description of how the A and B matrices are constructed please see Appendix A
Once the element stiffness matrix is found for every element the global stiffness matrix, K, which contains every degree of freedom in the model, can be constructed by simply adding the element stiffness matrices to the corresponding places in K.
Given a fully constrained system and set of known applied forces, F, the associated displacements can be found by solving the global stiffness equation:
u=K.F^(-1)
An example of a fully constrained system is a model with two nodes in which every degree of freedom at one node is fixed at 0. 
This surrogate model is designed to accept both Neumann (applied forces/moments) and Dirichlet (applied deflections/rotations) boundary conditions. To achieve this the global stiffness matrix and force vector must be adjusted. The force vector is first adjusted by setting any position in which the degree of freedom is directly defined to zero. Additionally, the columns and rows of the same degree of freedom are set to zero, except diagonals which are set as one. 

General Solver Functionality
Once the global stiffness matrix is correctly adjusted an initial guess of the displacements, u, is chosen (usually 0 in all locations). u is then fed into the material model, the strains are calculated via: 
ϵ=Bu
The resulting strains are used to find the stresses using Eq. X and finally the stresses are converted to nodal forces via Eq. Y.
The nodal force contribution from each element is added together in the global force vector. This is then compared to the adjusted applied forces vector, ignoring any rows which contain Dirichlet boundary conditions.
If the residual is greater than a set tolerance, then the initial guess for u has not reached equilibrium and a new guess for u is found using Eq. A. 
The above process is repeated iteratively using a Newton-Raphson convergence algorithm until a u that satisfies equilibrium is found.



A simple, linear elastic example of a FEA model and its associated ROM is given in figure X. Here the SCL is extracted from a notched bar. To replicate the behaviour seen in FEA, the ROM setup given in figure Xb is used. Only two nodes, one of which is completely locked

