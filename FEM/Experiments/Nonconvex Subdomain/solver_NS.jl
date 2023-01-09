using Gridap
using GridapGmsh
using Plots

# Here, we define the manufactured functions, we find the correspondng source terms for our manufactured solutions
# only the source terms f, f̂ are used to find our solution (u, û).

γ = 1

r(x) = sqrt(x[1]^2+x[2]^2)
θ(x) = x[2] < 0.0 ? -acos(x[1]/r(x)) : acos(x[1]/r(x))

dr(x) = (2/3)*r(x)^(-1/3)*sin((2*θ(x)+pi)/3)
dθ(x) = (2/3)*r(x)^(-1/3)*cos((2*θ(x)+pi)/3)

u(x) = x[1]< 0.0 && x[2] < 0.0 ? 0.0 : sin((2*θ(x)+π)/3)*r(x)^(2/3)
∇u(x) = x[1]< 0.0 && x[2] < 0.0 ? VectorValue(0.0,0.0) : VectorValue(dr(x), dθ(x))

Δu(x) = 0.0

û(x) = 0.0
∇û(x) = VectorValue(0.0, 0.0)
Δxû(x) = 0.0
Δyû(x) = 0.0

f(x) = -Δu(x)

M(x) = TensorValue(cos(θ(x)), sin(θ(x)), -sin(θ(x)), cos(θ(x))) #Transformation matrix, polar -> cartesian
ν₁(x) = M(x)⋅VectorValue(0, -1) # normal vector from Ω₁ over Γx in polar
ν₂(x) = M(x)⋅VectorValue(1, 0) # normal vector from Ω₁ over Γy in polar, note that we remove the minus because of our definition of Θ.

g₁(x) = -∇u(x)⋅ν₁(x)-γ*(u(x)-û(x))
g₂(x) = -∇u(x)⋅ν₂(x)-γ*(u(x)-û(x))

f̂₁(x) = -Δxû(x)+γ*(û(x)-u(x))
f̂₂(x) = -Δyû(x)+γ*(û(x)-u(x))


function solver(mesh, cell_type::Symbol=:tri, order::Int64=1; write::Bool=false)
    @assert cell_type ∈ (:quad, :tri)

    if cell_type == :quad
        path = pwd()*"/Quadrilateral mesh/quad_"
    end
    if cell_type == :tri
        path = pwd()*"/Triangular mesh/tri_"
    end

    model = GmshDiscreteModel(path*string(mesh)*".msh")

    order = order #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..
    #Define local basis functions
    reffe_u = ReferenceFE(lagrangian, Float64, order)
    reffe_û = ReferenceFE(lagrangian, Float64, order)


    #where should our solution live
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model, tags = ["interface"])
    Γx = BoundaryTriangulation(model, tags = ["interface_x"])
    Γy = BoundaryTriangulation(model, tags = ["interface_y"])
    ∂Ω = BoundaryTriangulation(model, tags = ["boundary"])


    V = TestFESpace(Ω,reffe_u, conformity=:H1, dirichlet_tags = ["boundary"] ) #we must add a constraint for neumann BC's
    V̂ = TestFESpace(Γ,reffe_û,dirichlet_tags = ["iface_left", "iface_right"], conformity=:H1)

    Y = MultiFieldFESpace([V,V̂])


    # Define global trial space
    U = TrialFESpace(V, u) #include u to enforce dirichlet tradtionally
    Û = TrialFESpace(V̂, û)
    X = MultiFieldFESpace([U,Û])

    #integration measure
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΓx = Measure(Γx,degree)
    dΓy = Measure(Γy,degree)


    a((u,û),(v,v̂)) = ∫(∇(u)⋅∇(v))dΩ + ∫(∇(û)⋅∇(v̂))dΓ + γ*∫((v-v̂)*(u-û))dΓ

    l((v,v̂)) = ∫(f*v)dΩ + ∫(f̂₁*v̂)dΓx + ∫(f̂₂*v̂)dΓy - ∫(g₁*v)dΓx-∫(g₂*v)dΓy

    # Build affine FE operator
    op = AffineFEOperator(a,l,X,Y)


    # Solve
    ls = LUSolver()
    solver = LinearFESolver(ls)
    wh = solve(solver, op)
    uh, ûh = wh

    eu = u - uh

    eul2 = sqrt(sum( ∫( eu*eu )dΩ ))
    euh1 = sqrt(sum( ∫( eu*eu + ∇(eu)⋅∇(eu) )dΩ ))
    println("error_u: ", euh1)

    eû = û - ûh
    eûl2 = sqrt(sum( ∫( eû*eû  )dΓ ))
    eûh1 = sqrt(sum( ∫( eû*eû + ∇(eû)⋅∇(eû) )dΓ ))
    println("error_û: ", eûh1)

    H1_error = sqrt(euh1^2+eûh1^2)
    L2_error = sqrt(eul2^2+eûl2^2)

    println("H1_tot: ", H1_error, "L2_tot: ", L2_error)

    if write
        writevtk(Ω,"../../vtu_files/error_u",cellfields=["eu" => eu])
        writevtk(Γ,"../../vtu_files/error_û",cellfields=["euhat" => eû])

        writevtk(Ω,"../../vtu_files/u",cellfields=["u" => u])
        writevtk(Ω, "../../vtu_files/omega", cellfields=["uh" => uh])
        writevtk(Γ, "../../vtu_files/gamma", cellfields=["uhat" => ûh])
    end


    return H1_error, L2_error
end


function conv_test(filename)
    #Create zero arrays to store errors
    P1 = zeros((2, 5))
    P2 = zeros((2, 5))
    Q1 = zeros((2, 5))
    Q2 = zeros((2, 5))

    for i in [1,2,3,4,5]
        #get errors for elementypes and interpolation order
        P1[:, i] .= solver(i, :tri, 1)
        P2[:, i] .= solver(i, :tri, 2)
        Q1[:, i] .= solver(i, :quad, 1)
        Q2[:, i] .= solver(i, :quad, 2)
    end

    h_tri = [0.3536, 0.1768, 0.0942, 0.0484, 0.0248] # these are extracted through inspection.
    h_quad = [0.25, 0.125, 0.0625, 0.0334, 0.0172]


    function slope(hs,errors)
      x = log10.(hs)
      y = log10.(errors)
      linreg = hcat(fill!(similar(x), 1), x) \ y
      round(linreg[2], digits = 2)
    end

    label1=["P1_H1 s="*string(slope(h_tri, P1[1, :])),
            "P1_L2 s="*string(slope(h_tri, P1[2, :])),
            "P2_H1 s="*string(slope(h_tri, P2[1, :])),
            "P2_L2 s="*string(slope(h_tri, P2[2, :]))]

    label2 =["Q1_H1 s="*string(slope(h_quad, Q1[1, :])),
            "Q1_L2 s="*string(slope(h_quad, Q1[2, :])),
            "Q2_H1 s="*string(slope(h_quad, Q2[1, :])),
            "Q2_L2 s="*string(slope(h_quad, Q2[2, :]))]

    plot(h_tri, [P1[1, :], P1[2, :], P2[1, :], P2[2, :]],
    xaxis=:log, yaxis=:log, label=reshape(label1, 1, length(label1)), shape=:utriangle, xlabel="h",ylabel="error norm")

    plot!(h_quad, [Q1[1, :], Q1[2, :], Q2[1, :], Q2[2, :]],
    xaxis=:log, yaxis=:log, label=reshape(label2, 1, length(label2)), shape=:rect, xlabel="h",ylabel="error norm")

    savefig(pwd()*"/Figures/"*filename)
end

#solver(1, :tri, write=true)
conv_test("conv_test_nonconvex.pdf")
