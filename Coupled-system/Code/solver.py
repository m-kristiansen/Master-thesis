import meshio
import dolfinx
from mpi4py import MPI  
import ufl        

def create_mesh(mesh, cell_type, prune_z=False):
    """Converts gmsh to meshio mesh, 
    2D mesh: prune_z = True
    3D mesh: prune_z = False 
    """
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:,:2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
    return out_mesh

# Read in mesh
msh = meshio.read("straight_line.msh")

# Create and save one file for the mesh, and one file for the facets 
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
line_mesh = create_mesh(msh, "line", prune_z=True)
meshio.write("mesh.xdmf", triangle_mesh)
meshio.write("mt.xdmf", line_mesh)

with dolfinx.io.XDMFFile(MPI.COMM_WORLD, 'mesh.xdmf', "r") as xdmf:
       surface = xdmf.read_mesh(name="Grid")

with dolfinx.io.XDMFFile(MPI.COMM_WORLD, 'mt.xdmf', "r") as xdmf:
        lines = xdmf.read_mesh(name="Grid")


print(lines.ufl_cell())

dx = ufl.Measure('dx', surface)
ds = ufl.Measure('ds', lines)(subdomain_data=subdomains)
a = dolfinx.fem.form(1*ds(0))
print(dolfinx.fem.assemble_scalar(a))