import gmsh
import sys
import numpy as np

def square_mesh(cls, type, name, view = False):

    gmsh.initialize(["", "-clscale", str(cls)])
    gmsh.model.add('box')

    model = gmsh.model
    occ = model.occ


    #Add points for square box
    bp1 = occ.addPoint(-1, -1, 0)
    bp2 = occ.addPoint(1, -1, 0)
    bp3 = occ.addPoint(1, 1, 0)
    bp4 = occ.addPoint(-1, 1, 0)

    #Add points for interface
    ml = occ.addPoint(-1, 0, 0)
    cp = occ.addPoint(0.0, 0.0, 0)
    mt = occ.addPoint(0, 1, 0)
    mr = occ.addPoint(1, 0, 0)
    mb = occ.addPoint(0, -1, 0)

    # outer lines
    l1 = occ.addLine(bp1, mb)
    l2 = occ.addLine(mb, bp2)
    l3 = occ.addLine(bp2, mr)
    l4 = occ.addLine(mr, bp3)
    l5 = occ.addLine(bp3, mt)
    l6 = occ.addLine(mt, bp4)
    l7 = occ.addLine(bp4, ml)
    l8 = occ.addLine(ml, bp1)

    #center lines
    sl1 = occ.addLine(cp, mb)
    sl2 = occ.addLine(cp, ml)
    sl3 = occ.addLine(cp, mr)
    sl4 = occ.addLine(cp, mt)

    occ.synchronize()

    bottomleft = occ.addCurveLoop([l1, -sl1, sl2, l8])
    bottomleft_surf = occ.addPlaneSurface([bottomleft])
    occ.synchronize()
    model.addPhysicalGroup(2, [bottomleft_surf], 1)
    gmsh.model.setPhysicalName(2, 1, "bottomleft_surf")

    bottomright = occ.addCurveLoop([l2, l3, -sl3, sl1])
    bottomright_surf = occ.addPlaneSurface([bottomright])
    occ.synchronize()
    model.addPhysicalGroup(2, [bottomright_surf], 2)
    gmsh.model.setPhysicalName(2, 2, "bottomright_surf")

    topright = occ.addCurveLoop([sl3, l4, l5, -sl4])
    topright_surf = occ.addPlaneSurface([topright])
    occ.synchronize()
    model.addPhysicalGroup(2, [topright_surf], 3)
    gmsh.model.setPhysicalName(2, 3, "topright_surf")

    topleft = occ.addCurveLoop([-sl2, sl4, l6, l7])
    topleft_surf = occ.addPlaneSurface([topleft])
    occ.synchronize()
    model.addPhysicalGroup(2, [topleft_surf], 4)
    gmsh.model.setPhysicalName(2, 4, "topleft_surf")

    model.addPhysicalGroup(1, [l1, l2, l3, l4, l5, l6, l7, l8], 1)
    gmsh.model.setPhysicalName(1, 1, "boundary")

    model.addPhysicalGroup(0, [bp1, bp2, bp3, bp4, ml, mt, mr, mb], 1)
    gmsh.model.setPhysicalName(0, 1, "boundary")

    model.addPhysicalGroup(1, [-sl2, sl1], 2)
    gmsh.model.setPhysicalName(1, 2, "interface")

    model.addPhysicalGroup(1, [-sl2], 3)
    gmsh.model.setPhysicalName(1, 3, "interface_x")

    model.addPhysicalGroup(0, [ml, cp], 5)
    gmsh.model.setPhysicalName(1, 5, "interface_x")

    model.addPhysicalGroup(1, [sl1], 4)
    gmsh.model.setPhysicalName(1, 4, "interface_y")

    model.addPhysicalGroup(0, [cp, mb], 6)
    gmsh.model.setPhysicalName(1, 6, "interface_y")

    model.addPhysicalGroup(0, [ml, mb, cp], 4)
    gmsh.model.setPhysicalName(0, 4, "interface")

    model.addPhysicalGroup(0, [ml], 7)
    gmsh.model.setPhysicalName(0, 7, "iface_left")

    model.addPhysicalGroup(0, [mb], 8)
    gmsh.model.setPhysicalName(0, 8, "iface_right")

    model.addPhysicalGroup(0, [cp], 9)
    gmsh.model.setPhysicalName(0, 9, "center_point")

    occ.synchronize()

    model.mesh.setTransfiniteSurface(bottomleft_surf, "Left")
    model.mesh.setTransfiniteSurface(bottomright_surf, "Left")
    model.mesh.setTransfiniteSurface(topleft_surf, "Right")
    model.mesh.setTransfiniteSurface(topright_surf, "Left")

    if type=="quad":
        gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.mesh.generate(2)

    if view:
        gmsh.fltk.initialize()
        gmsh.fltk.run()

    #gmsh.write(name)
    gmsh.finalize()


square_mesh(0.125, type = "tri", name = "quad_5.msh", view = True)
