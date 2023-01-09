import gmsh
import sys


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


#bottom surface
bottom_loop = gmsh.model.geo.addCurveLoop([l1, l2, -gl, l6])
bottom_surface = gmsh.model.geo.addPlaneSurface([bottom_loop])
gmsh.model.geo.addPhysicalGroup(2, [bottom_surface], 1)
gmsh.model.setPhysicalName(2, 1, "bottom_surface") #(dim, tag, label)
gmsh.model.geo.addPhysicalGroup(1, [l1, l2, l6], 1)
gmsh.model.setPhysicalName(1, 1, "boundary_bottom_surface")
gmsh.model.geo.addPhysicalGroup(0, [bp1, bp2, ip1, ip2], 1)
gmsh.model.setPhysicalName(0, 1, "boundary_bottom_surface")


#top surface
top_loop = gmsh.model.geo.addCurveLoop([l3, l4, l5, gl])
top_surface = gmsh.model.geo.addPlaneSurface([top_loop])
gmsh.model.geo.addPhysicalGroup(2, [top_surface], 2)
gmsh.model.setPhysicalName(2, 2, "top_surface") #(dim, tag, label)
gmsh.model.geo.addPhysicalGroup(1, [l3, l4, l5], 2)
gmsh.model.setPhysicalName(1, 2, "boundary_top_surface")
gmsh.model.geo.addPhysicalGroup(0, [ip1, ip2, bp3, bp4], 2)
gmsh.model.setPhysicalName(0, 2, "boundary_top_surface")
#
# boxloop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6])
# boxsurface = gmsh.model.geo.addPlaneSurface([boxloop])
# gmsh.model.geo.addPhysicalGroup(2, [boxsurface], 3)
# gmsh.model.setPhysicalName(2, 3, "box surface") #(dim, tag, label)


gmsh.model.geo.addPhysicalGroup(1, [gl], 3)
gmsh.model.setPhysicalName(1, 3, "Graph")
gmsh.model.geo.addPhysicalGroup(0, [ip1, ip2], 3)
gmsh.model.setPhysicalName(0, 3, "Graph")

gmsh.model.geo.addPhysicalGroup(0, [ip1], 4) #first graph node is tagged after box corners i.e. 5
gmsh.model.setPhysicalName(0, 4, "starting_point")
gmsh.model.geo.addPhysicalGroup(0, [ip2], 5)
gmsh.model.setPhysicalName(0, 5, "end_point")

gmsh.model.geo.synchronize() # Sync CAD representation

#embed the graph in the mesh
gmsh.model.mesh.embed(1, [gl], 2, top_surface) # (dim object_1, object_1, dim object_2, object_2)
gmsh.model.mesh.embed(1, [gl], 2, bottom_surface)
gmsh.model.mesh.generate(2)

gmsh.fltk.run()
gmsh.write('mesh_hsplit.msh')
gmsh.finalize()
