using Gridap
using GridapGmsh
using Plots

# Here, we define the manufactured functions, we find the correspondng source terms for our manufactured solutions
# only the source terms f, f̂ are used to find our solution (u, û) while the gradients will be used later to calculate the

γ = 1


u(x) = sin(abs(x[2]))
∇u(x) = VectorValue(0, sign(x[2])*cos(x[2]))
Δu(x) = -sin(abs(x[2]))

û(x) = 0
∇û(x) = VectorValue(0, 0)
Δû(x) = 0

∇u₁(x) = VectorValue(0, cos(x[2]))
∇u₂(x) = VectorValue(0, -cos(x[2]))

f(x) = -Δu(x)
g(x) = -∇u₁(x)⋅VectorValue(0, -1)-∇u₂(x)⋅VectorValue(0,1)-γ*(u(x)-û(x))
f̂(x) = -Δû(x)+γ*(û(x)-u(x))



import Gridap: ∇
∇(::typeof(u)) = ∇u #get grad u from here

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
    Γ = BoundaryTriangulation(model, tags = ["interface", "iface_left", "iface_right"])

    dirichlet_tags = ["top", "bottom", "bottom_left", "bottom_right",
    "top_left", "top_right", "iface_left", "iface_right", "ll", "lr", "ul", "ur"]

    V = TestFESpace(Ω,reffe_u, conformity=:H1, dirichlet_tags = dirichlet_tags) #we must add a constraint for neumann BC's
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


    a((u,û),(v,v̂)) = ∫(∇(u)⋅∇(v))dΩ + ∫(∇(û)⋅∇(v̂))dΓ + γ*∫((v-v̂)*(u-û))dΓ

    l((v,v̂)) = ∫(f*v)dΩ+ ∫(f̂*v̂)dΓ - ∫(g*v)dΓ

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

    h_tri = [0.1767, 0.0914, 0.048, 0.0246, 0.0124] # these are extracted through inspection.
    h_quad = [0.125, 0.0625, 0.0333, 0.0172, 0.00877]


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

solver(1, :tri, write=true)
#conv_test("conv_test_sinabs.pdf")
