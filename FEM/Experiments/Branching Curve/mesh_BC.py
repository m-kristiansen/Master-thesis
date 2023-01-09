import gmsh
import sys
import numpy as np
import os


def create_mesh(cls, type, name, view=True, save = True):
    gmsh.initialize(["", "-clscale", str(cls)])
    gmsh.model.add('box')
    #Add points for square box
    bp1 = gmsh.model.geo.addPoint(-1, -1, 0)
    bp2 = gmsh.model.geo.addPoint(1, -1, 0)
    bp3 = gmsh.model.geo.addPoint(1, 1, 0)
    bp4 = gmsh.model.geo.addPoint(-1, 1, 0)

    #Add point for graph intersection with box edge
    ip1 = bp1
    ip2 = gmsh.model.geo.addPoint(1, 0.0, 0)

    #graph points
    ip3 = gmsh.model.geo.addPoint(0.0, 0.0, 0)
    ip4 = gmsh.model.geo.addPoint(0.0, 0.8, 0)


    # add lines to gmsh
    l1 = gmsh.model.geo.addLine(bp1, bp2)
    l2 = gmsh.model.geo.addLine(bp2, ip2)
    l3 = gmsh.model.geo.addLine(ip2, bp3)
    l4 = gmsh.model.geo.addLine(bp3, bp4)
    l5 = gmsh.model.geo.addLine(bp4, bp1)

    #graph lines
    gl1 = gmsh.model.geo.addLine(ip1, ip3)
    gl2 = gmsh.model.geo.addLine(ip3, ip2)
    gl3 = gmsh.model.geo.addLine(ip3, ip4)



    boxloop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5])
    gmsh.model.geo.addPhysicalGroup(1, [l1, l2, l3, l4, l5], 1)
    gmsh.model.setPhysicalName(1, 1, "boundary")
    gmsh.model.geo.addPhysicalGroup(0, [bp1, bp2, bp3, bp4, ip1, ip2], 1)
    gmsh.model.setPhysicalName(0, 1, "boundary")


    boxsurface = gmsh.model.geo.addPlaneSurface([boxloop])
    gmsh.model.geo.addPhysicalGroup(2, [boxsurface], 1)
    gmsh.model.setPhysicalName(2, 1, "My surface") #(dim, tag, label)


    gmsh.model.geo.addPhysicalGroup(1, [gl1, gl2, gl3], 2)
    gmsh.model.setPhysicalName(1, 2, "Graph")

    gmsh.model.geo.addPhysicalGroup(1, [gl1], 3)
    gmsh.model.setPhysicalName(1, 3, "Graph_xy")

    gmsh.model.geo.addPhysicalGroup(1, [gl2], 4)
    gmsh.model.setPhysicalName(1, 4, "Graph_x")

    gmsh.model.geo.addPhysicalGroup(1, [gl3], 5)
    gmsh.model.setPhysicalName(1, 5, "Graph_y")

    gmsh.model.geo.addPhysicalGroup(0, [ip1, ip2, ip3, ip4], 2)
    gmsh.model.setPhysicalName(0, 2, "Graph")

    gmsh.model.geo.addPhysicalGroup(0, [ip1], 3) #first graph node is tagged after box corners i.e. 5
    gmsh.model.setPhysicalName(0, 3, "starting_point")
    gmsh.model.geo.addPhysicalGroup(0, [ip2, ip4], 4)
    gmsh.model.setPhysicalName(0, 4, "end_points")


    gmsh.model.geo.synchronize() # Sync CAD representation

    #embed the graph in the mesh
    gmsh.model.mesh.embed(1, [gl1, gl2, gl3], 2, boxsurface) # (dim object_1, object_1, dim object_2, object_2)

    if type == "quad":
        gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.mesh.generate(2)


    #Here we sort the element_qualities
    if type == 'quad':
        eltype = 3
    if type == 'tri':
        eltype = 2
    eltypes, _ = gmsh.model.mesh.getElementsByType(eltype)
    qualitites = gmsh.model.mesh.getElementQualities(eltypes, "gamma")
    element_qualities = np.array([eltypes, qualitites])
    sorted_qualities = np.sort(element_qualities[1, :], axis = 0)
    _, idx = np.where(element_qualities == sorted_qualities[0])
    print(element_qualities[0, idx])
    print("eta = ", sorted_qualities[0])

    if view:
        gmsh.fltk.initialize()
        gmsh.fltk.run()

    if save==True:
        gmsh.write(name)

    gmsh.finalize()

create_mesh(0.25, "tri", "tri_5.msh", save = False)
