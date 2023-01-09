using Gridap
using GridapGmsh
using Plots

# Here, we define the manufactured functions, we find the correspondng source terms for our manufactured solutions
# only the source terms f, f̂ are used to find our solution (u, û) while the gradients will be used later to calculate the

γ = 1

#
# u(x) = sin(x[1])+cos(x[2])
# ∇u(x) = VectorValue(cos(x[1]), -sin(x[2]))
# Δu(x) = -sin(x[1])-cos(x[2])
#
# û(x) = cos(x[1])+sin(x[2])
# ∇û(x) = VectorValue(-sin(x[1]), cos(x[2]))
# Δû(x) = -cos(x[1])-sin(x[2])

u(x) = x[1]^3-x[2]^3
∇u(x) = VectorValue(3*x[1]^2, -3*x[2]^2)
Δu(x) = 6*x[1]-6*x[2]

û(x) = x[1]^3+x[2]^3
∇û(x) = VectorValue(3*x[1]^2, 3*x[2]^2)

f(x) = -Δu(x)
g(x) = -γ*(u(x)-û(x))

∂xx(x) = 6*x[1]
∂yy(x) = 6*x[2]

f̂₁(x) = -VectorValue(∂xx(x)/sqrt(2), ∂yy(x)/sqrt(2))⋅VectorValue(1/sqrt(2),1/sqrt(2))+γ*(û(x)-u(x))
f̂₂(x) = -∂xx(x)+γ*(û(x)-u(x))
f̂₃(x) = -∂yy(x)+γ*(û(x)-u(x))



#import Gridap: ∇
#∇(::typeof(u)) = ∇u #get grad u from here

function solver(mesh, order ;write::Bool=false)

    path = pwd()*"/Triangular mesh/tri_"

    model = GmshDiscreteModel(path*string(mesh)*".msh")

    order = order #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

    #Define local basis functions
    reffe_u = ReferenceFE(lagrangian, Float64, order)
    reffe_û = ReferenceFE(lagrangian, Float64, order)


    #where should our solution live
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model, tags = ["Graph"])
    Γxy = BoundaryTriangulation(model, tags = ["Graph_xy"])
    Γx = BoundaryTriangulation(model, tags = ["Graph_x"])
    Γy = BoundaryTriangulation(model, tags = ["Graph_y"])

    V = TestFESpace(Ω,reffe_u, conformity=:H1, dirichlet_tags = "boundary") #we must add a constraint for neumann BC's
    V̂ = TestFESpace(Γ,reffe_û,dirichlet_tags = ["starting_point", "end_points"], conformity=:H1)

    Y = MultiFieldFESpace([V,V̂])


    # Define global trial space
    U = TrialFESpace(V, u) #include u to enforce dirichlet tradtionally
    Û = TrialFESpace(V̂, û)
    X = MultiFieldFESpace([U,Û])

    #integration measure
    degree = order+1
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΓxy = Measure(Γxy,degree)
    dΓx = Measure(Γx,degree)
    dΓy = Measure(Γy,degree)


    a((u,û),(v,v̂)) = ∫(∇(u)⋅∇(v))dΩ + ∫(∇(û)⋅∇(v̂))dΓ + γ*∫((v-v̂)*(u-û))dΓ

    l((v,v̂)) = ∫(f*v)dΩ + ∫(f̂₁*v̂)dΓxy + ∫(f̂₂*v̂)dΓx + ∫(f̂₃*v̂)dΓy  - ∫(g*v)dΓ

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

    P1 = zeros((2, 5))
    P2 = zeros((2, 5))
    for i in [1,2,3,4,5]
        P1[:, i] .= solver(i, 1)
        P2[:, i] .= solver(i, 2)
    end

    h_tri = [0.3504, 0.1572, 0.0773, 0.0437, 0.0214] #these are extracted through inspection

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

    plot(h_tri, [P1[1, :], P1[2, :], P2[1, :], P2[2, :]],
        xaxis=:log, yaxis=:log, label=reshape(label1, 1, length(label1)), shape=:utriangle, xlabel="h",ylabel="error norm")


    savefig(pwd()*"/Figures/"*filename)
end


conv_test("conv_test_BC.pdf")
#solver(1, 1, write = true)
