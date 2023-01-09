import gmsh
import sys
import numpy as np

def square_mesh(structured, cls, name, view = False):

    gmsh.initialize(["",  "-clscale", str(cls)])
    gmsh.model.add('box')

    model = gmsh.model
    occ = model.occ

    lc = 0.10
    #Add points for square box
    bp1 = occ.addPoint(0, 0, 0, lc)
    bp2 = occ.addPoint(1, 0, 0, lc)
    bp3 = occ.addPoint(1, 1, 0, lc)
    bp4 = occ.addPoint(0, 1, 0, lc)


    # add lines to gmsh
    l1 = occ.addLine(bp1, bp2)
    l2 = occ.addLine(bp2, bp3)
    l3 = occ.addLine(bp3, bp4)
    l4 = occ.addLine(bp4, bp1)


    boxloop = occ.addCurveLoop([bp1, bp2, bp3, bp4])
    boxsurface = occ.addPlaneSurface([boxloop])

    occ.synchronize()

    gmsh.model.addPhysicalGroup(2, [boxsurface], 1)
    gmsh.model.setPhysicalName(2, 1, "My surface") #(dim, tag, label)

    gmsh.model.addPhysicalGroup(1, [bp1, bp2, bp3, bp4], 1)
    gmsh.model.setPhysicalName(1, 1, "boundary")
    gmsh.model.addPhysicalGroup(0, [bp1, bp2, bp3, bp4], 2)
    gmsh.model.setPhysicalName(0, 2, "boundary")

    occ.synchronize() # Sync CAD representation

    eltypes, _ = gmsh.model.mesh.getElementsByType(2)
    qualitites = gmsh.model.mesh.getElementQualities(eltypes, "gamma")

    element_qualities = np.array([eltypes, qualitites])
    print(np.sort(element_qualities[1, :], axis = 0))

    if structured:
        gmsh.model.mesh.setTransfiniteSurface(boxsurface)

    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.mesh.generate(2)

    if view:
        gmsh.fltk.initialize()
        gmsh.fltk.run()

    gmsh.write(name)
    gmsh.finalize()


square_mesh(structured = True, cls = 1.0, name = "mesh_1.msh", view = True)
