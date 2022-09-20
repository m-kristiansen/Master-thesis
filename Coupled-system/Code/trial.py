import gmsh
import sys
import dolfinx
from dolfinx.io.gmshio import model_to_mesh
from mpi4py import MPI   
import ufl

gmsh.initialize()
gmsh.model.add('line')
#Add points for square box
bp1 = gmsh.model.geo.addPoint(0, 0, 0)
bp2 = gmsh.model.geo.addPoint(1, 0, 0)
bp3 = gmsh.model.geo.addPoint(1, 1, 0)
bp4 = gmsh.model.geo.addPoint(0, 1, 0)

#Add point for graph intersection with box edge
ip1 = gmsh.model.geo.addPoint(0.0, 0.5, 0)
ip2 = gmsh.model.geo.addPoint(1, 0.5, 0)

# add lines to gmsh
l1 = gmsh.model.geo.addLine(bp1, bp2) 
l2 = gmsh.model.geo.addLine(bp2, ip2) 
l3 = gmsh.model.geo.addLine(ip2, bp3) 
l4 = gmsh.model.geo.addLine(bp3, bp4) 
l5 = gmsh.model.geo.addLine(bp4, ip1) 
l6 = gmsh.model.geo.addLine(ip1, bp1)

#graph line
gl = gmsh.model.geo.addLine(ip1, ip2)

boxloop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6])
gmsh.model.geo.addPhysicalGroup(1, [l1, l2, l3, l4, l5, l6], 1)
gmsh.model.setPhysicalName(1, 1, "boundary")
boxsurface = gmsh.model.geo.addPlaneSurface([boxloop])
gmsh.model.geo.addPhysicalGroup(2, [boxsurface], 1)
gmsh.model.setPhysicalName(2, 1, "My surface") #(dim, tag, label)


gmsh.model.geo.addPhysicalGroup(1, [gl], 2)
gmsh.model.setPhysicalName(1, 2, "Graph")
gmsh.model.geo.addPhysicalGroup(0, [5], 4) #first graph node is tagged after box corners i.e. 5
gmsh.model.setPhysicalName(0, 3, "starting_point")

gmsh.model.geo.synchronize() # Sync CAD representation

#embed the graph in the mesh
gmsh.model.mesh.embed(1, [gl], 2, boxsurface) # (dim object_1, object_1, dim object_2, object_2)
gmsh.model.mesh.generate(2)

#gmsh.fltk.run()
mesh, surface_tags, line_tags = model_to_mesh(
    gmsh.model, comm=MPI.COMM_WORLD, rank=0, gdim=2)
gmsh.finalize()

dx = ufl.Measure("dx")(mesh)
dS = ufl.Measure("dS")(subdomain_data=line_tags)
a = dolfinx.fem.form(1*dS)
print(dolfinx.fem.assemble_scalar(a))