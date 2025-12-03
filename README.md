Surrogate Model Formulation
Single Element Formulation
Consider a single beam element with 2 nodes. This element has 8 degrees of freedom at each node which allow for a 6 dimensional stress tensor to be calculated for any given applied load or displacement. The 8 degrees of freedom are linked by considering each element as three 3D beams as shown in Figure 2 1.
Figure 2 1: Figure showing one element and the beams it is made of. The dofs for one end of the beam are given.

The three beams correspond to an axial beam with 4 degrees of freedom (u_xx,u_xy,u_xz,θ_xz), a transverse beam (u_yy,θ_yz), and an out-of-plane beam (u_zz,θ_zx). The forces applied to an element, f_e,can be related to the resulting deflections and rotations, u_e, via the element stiffness matrix, k_e, as follows:

f_e=k_e u_e	Equation 1


The element stiffness matrix is constructed using a set of three governing equations:
	The degrees of freedom related to the element, u_e, can be related to the deformation (represented as strain), ϵ via Equation 2:
ϵ=Bu_e	Equation 2

Where B is a matrix that defines the relationship between u_e and ϵ. The Full B matrix used for every element is as follows:
[■(ϵ_xx@ϵ_yy@ϵ_zz  @ϵ_xy@ϵ_xz@ϵ_yz )]=1/L {■(1&0&0&dn&0&0&0&0&-1&0&0&-dn&0&0&0&0@0&0&0&0&1&0&dn&0&0&0&0&0&-1&0&-dn&0 @0&0&0&0&0&1&0&dn&0&0&0&0&0&-1&0&-dn@0&1&0&d_tn^n L_x&0&0&0&0&0&-1&0&-d_tn^n L_x&0&0&0&0@0&0&1&0&0&0&0&0&0&0&-1&0&0&0&0&0@0&0&0&0&0&0&0&0&0&0&0&0&0&0&0&0)}[■(u_(xx,i)@u_(xy,i)@u_(xz,i)@θ_(xz,i)@u_(yy,i)@u_(zz,i)@θ_(yy,i)@θ_(zz,i)@u_(xx,j)@u_(xy,j)@u_(xz,j)@θ_(xz,j)@u_(yy,j)@u_(zz,j)@θ_(yy,j)@θ_(zz,j) )]	Equation 3

Where i and j represent the top and bottom node of the element. L_x,L_y  and L_z represent the length of each beam in the element. d_th^n is the distance from the element connection point on the node to the neutral axis. Figure 2 2 gives a visual depiction of the meaning of this parameter. 
The construction of the B matrix is outlined in Appendix A

Figure 2 2: Figure describing what d_x^n is.
	The strain is related to the internal forces (represented as stress) via Equation 4
σ=Cϵ	Equation 4

Where σ_ij is the stress and C_ijkl is the 6x6 material model. 

	The stress is related to the applied forces via Equation 5:
f_e=A^T σ	Equation 5

	Where A^T is a matrix that defines the relationship between f_e and σ. In this case A^T=B^T.

With these governing equations it is possible link an applied force to a relative change in the length of the element, assuming the problem is fully constrained (i.e. one node is pinned).
