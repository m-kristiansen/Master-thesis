import gmsh
import sys
import numpy as np

def square_mesh(cls, structured, name, view = False):

    gmsh.initialize(["", "-clscale", str(cls)])
    gmsh.model.add('box')

    model = gmsh.model
    occ = model.occ


    lc = .15
    #Add edge points
    bp2 = occ.addPoint(1, -1, 0, lc)
    bp3 = occ.addPoint(1, 1, 0, lc)
    bp4 = occ.addPoint(-1, 1, 0, lc)

    #Add points for interface
    ml = occ.addPoint(-1, 0, 0, lc)
    cp = occ.addPoint(0, 0, 0, lc)
    mt = occ.addPoint(0, 1, 0, lc)
    mr = occ.addPoint(1, 0, 0, lc)
    mb = occ.addPoint(0, -1, 0, lc)

    # outer lines
    l3 = occ.addLine(mb, bp2)
    l4 = occ.addLine(bp2, mr)
    l5 = occ.addLine(mr, bp3)
    l6 = occ.addLine(bp3, mt)
    l7 = occ.addLine(mt, bp4)
    l8 = occ.addLine(bp4, ml)

    #center lines
    sl1 = occ.addLine(cp, mb)
    sl2 = occ.addLine(cp, ml)
    sl3 = occ.addLine(cp, mr)
    sl4 = occ.addLine(cp, mt)

    occ.synchronize()

    if structured:
        topright = occ.addCurveLoop([sl3, l5, l6, -sl4])
        topright_surf = occ.addPlaneSurface([topright])
        occ.synchronize()
        model.addPhysicalGroup(2, [topright_surf], 1)
        gmsh.model.setPhysicalName(2, 1, "topright_surf")

        bottomright = occ.addCurveLoop([l3, l4, -sl3, sl1])
        bottomright_surf = occ.addPlaneSurface([bottomright])
        occ.synchronize()
        model.addPhysicalGroup(2, [bottomright_surf], 2)
        gmsh.model.setPhysicalName(2, 2, "bottomright_surf")


        topleft = occ.addCurveLoop([-sl2, sl4, l7, l8])
        topleft_surf = occ.addPlaneSurface([topleft])
        occ.synchronize()
        model.addPhysicalGroup(2, [topleft_surf], 4)
        gmsh.model.setPhysicalName(2, 4, "topleft_surf")

        model.mesh.setTransfiniteSurface(topright_surf, "Left")
        model.mesh.setTransfiniteSurface(bottomright_surf, "Left")
        model.mesh.setTransfiniteSurface(topleft_surf, "Right")

    else:
        curveloop = occ.addCurveLoop([-sl2, sl1, l3, l4, l5, l6, l7, l8])
        surface = occ.addPlaneSurface([curveloop])
        occ.synchronize()
        model.addPhysicalGroup(2, [surface], 1)
        gmsh.model.setPhysicalName(2, 1, "surface")


    model.addPhysicalGroup(1, [-sl2, sl1, l3, l4, l5, l6, l7, l8], 1)
    gmsh.model.setPhysicalName(1, 1, "boundary")

    model.addPhysicalGroup(0, [ml, mb, bp2, mr, bp3, mt, bp4], 1)
    gmsh.model.setPhysicalName(0, 1, "boundary")

    occ.synchronize()

    gmsh.model.mesh.generate(2)

    if view:
        gmsh.fltk.initialize()
        gmsh.fltk.run()

    gmsh.write(name)
    gmsh.finalize()


square_mesh(1, structured = True, name = "nonconvex.vtk", view = True)
