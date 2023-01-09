using GridapGmsh: gmsh, GmshDiscreteModel
"""Pretty much this
         ____________(1, 0.5)
         |__Omega1___|
         |__Omega2___|
(-1, -0.5)
"""
function split_square_mesh(clscale::Real, cell_type::Symbol=:tri; structured::Bool=false, view::Bool=false)
    @assert clscale > 0
    @assert cell_type ∈ (:quad, :tri)

    !isdir(".msh_cache") && mkdir(".msh_cache")

    gmsh.initialize(["", "-clscale", string(clscale)])

    model = gmsh.model
    occ = model.occ

    make_line = (p, q) -> occ.addLine(p, q)

    points, normals = [], Dict{String, Vector{}}()
    push!(points, [occ.addPoint(0, 0, 0),
                occ.addPoint(0, 0.5, 0),
                occ.addPoint(1, 0.5, 0),
                occ.addPoint(1, 0, 0),
                occ.addPoint(1, -0.5, 0),
                occ.addPoint(0, -0.5, 0)]...)

    npts = length(points)
    lines = [occ.addLine(points[p], points[1 + p%npts]) for p ∈ eachindex(points)]
    # Add the interface
    append!(lines, occ.addLine(points[1], points[4]))

    top_lines = [lines[1], lines[2], lines[3], -lines[end]]
    top_loop = occ.addCurveLoop(top_lines)
    top_surf = occ.addPlaneSurface([top_loop])

    bottom_lines = [lines[4], lines[5], lines[6], lines[end]]
    bottom_loop = occ.addCurveLoop(bottom_lines)
    bottom_surf = occ.addPlaneSurface([bottom_loop])

    occ.synchronize()

    # Mark the boundary points of the interface top_loop
    names = ("iface_left", "ul", "ur", "iface_right", "lr", "ll")
    for (tag, (point, name)) ∈ enumerate(zip(points, names))
        iface_tag = model.addPhysicalGroup(0, [point], tag)
        gmsh.model.setPhysicalName(0, iface_tag, name)
    end

    top_surf_group = model.addPhysicalGroup(2, [top_surf], 1)
    gmsh.model.setPhysicalName(2, top_surf_group, "top_surface")

    bottom_surf_group = model.addPhysicalGroup(2, [bottom_surf], 2)
    gmsh.model.setPhysicalName(2, bottom_surf_group, "bottom_surface")

    names = ("top_left", "top", "top_right", "bottom_right", "bottom", "bottom_left", "interface")
    for (tag, (line, name)) ∈ enumerate(zip(lines, names))
        line_group = model.addPhysicalGroup(1, [line], tag)
        gmsh.model.setPhysicalName(1, line_group, name)
    end

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    structured && model.mesh.setTransfiniteSurface(top_surf)
    structured && model.mesh.setTransfiniteSurface(bottom_surf)

    cell_type == :quad && gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.mesh.generate(2)


    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    name = "quad_1.msh"
    #gmsh.write(name)

    gmsh.finalize()

    (name, normals)
end

split_square_mesh(0.0625, :quad, view = true, structured = true)
