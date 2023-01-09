"""
This script uses the fenics functionality to solve the poisson equation 
- Δu = f    on Ω
u = u_d     on ∂Ω

We construct a manufactured solution u_e = 1 + x^2 + 2*y^2, and see that setting the source term f = -6 then solves the equation.
"""

from mpi4py import MPI                      #import message passage interface, allows us to define how to proscess the program
from dolfinx import mesh, fem               #import fenichs functionality
import ufl                                  #allows us to use "close to mathematics notation"
from petsc4py.PETSc import ScalarType  
import numpy as np     

#We create a simple unit square mesh, MPI.COMM_WORLD allows us to define the number of processors we would like to use when running the script. 
#e.g. mpirun -n 2 python3 diffusion.py runs the program in parallel on two processors
# (8,8) = (nx, ny)
# mesh.CellType.quadrilateral uses square cell shapes.  
domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8, mesh.CellType.quadrilateral)

#Define the function space on our domain, the type of element and degree. 
# CG denotes Lagrange type elements. Whence our solution will be a linear functions between the verticies of our square elements. 
V = fem.FunctionSpace(domain, ("CG", 1)) 

#We want to specify our boundary condition by our manufactured solution. We define a function on our grid as follows
uD = fem.Function(V)
uD.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)

# We want to enforce the values to the solution on the boundary by the function above. We then first need to extract the boundary DOFs of our domain.
tdim = domain.topology.dim 
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)          #Boundary line segments

boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets) 
bc = fem.dirichletbc(uD, boundary_dofs)                                 #Enforce uD on boundary DOF's

#We dinstinguish between trial and test spaces, but in the present 
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

f = fem.Constant(domain, ScalarType(-6)) # allows us to update the f later without recomputing our variational form, and defines it as scaler to speed up computations. 

#Weak formulation of the poisson equation
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

#We define our problem as our weak formulation with the corresponding boundary conditions. This corresponds to solving a linear system Ax = b
# which we solve using LU-factorization, using petsc we can maintain our functionality for parallel programming
problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

#We now want to compare our solution uh to our exact solution u_e
# We do this by interpolating our exact solution into the P2-function space 
V2 = fem.FunctionSpace(domain, ("CG", 2))
uex = fem.Function(V2)
uex.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)

#Compute the L2 error
L2_error = fem.form(ufl.inner(uh - uex, uh - uex) * ufl.dx)              #ufl implementation of the formula
error_local = fem.assemble_scalar(L2_error)                              #get scalar values
error_L2 = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))       #gather all proscesses 

error_max = np.max(np.abs(uD.x.array-uh.x.array))                        # maximum error, note that we compare with uD since clearly uex has more DOF's
# Only print the error on one process
if domain.comm.rank == 0:
    print(f"Error_L2 : {error_L2:.2e}")
    print(f"Error_max : {error_max:.2e}")

