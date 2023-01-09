using Gridap
using GridapGmsh
using Plots


function solver(γ)

    u(x) = cos(x[1])+sin(x[2])
    ∇u(x) = VectorValue(-sin(x[1]),cos(x[2]))
    Δu(x) = -cos(x[1])-sin(x[2])

    û(x) = sin(x[1])
    ∇û(x) = VectorValue(cos(x[1]),0)
    Δû(x) = -sin(x[1])

    # u(x) = x[1]-x[2]
    # ∇u(x) = VectorValue(1,-1)
    # Δu(x) = 0
    #
    # û(x) = x[1]^2
    # ∇û(x) = VectorValue(2*x[1],0)
    # Δû(x) = 2

    f(x) = -Δu(x)
    g(x) = -γ*(u(x)-û(x))
    f̂(x) = -Δû(x)+γ*(û(x)-u(x))

    model = GmshDiscreteModel("/Users/martinkristiansen/Desktop/Master thesis/FEM/Experiments/Straight Line/Triangular mesh/tri_1.msh")
    #model = simplexify(model)

    order = 2 #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

    #Define local basis functions
    reffe_u = ReferenceFE(lagrangian, Float64, order)
    reffe_û = ReferenceFE(lagrangian, Float64, order)

    boundary = ["top", "bottom", "bottom_left", "bottom_right",
    "top_left", "top_right", "iface_left", "iface_right", "ll", "lr", "ul", "ur"]

    graph = ["interface", "iface_left", "iface_right"]

    #where should our solution live
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model, tags = graph)
    ∂Ω = BoundaryTriangulation(model, tags = boundary)

    V = TestFESpace(Ω,reffe_u, conformity=:H1, dirichlet_tags = boundary) #we must add a constraint for neumann BC's
    V̂ = TestFESpace(Γ,reffe_û,dirichlet_tags = ["iface_left", "iface_right"], conformity=:H1)

    Y = MultiFieldFESpace([V,V̂])

    # Define global trial space
    U = TrialFESpace(V, u) #include u to enforce dirichlet tradtionally
    Û = TrialFESpace(V̂, û) # 1 dirichlet on start of graph
    X = MultiFieldFESpace([U,Û])

    #integration measure
    degree = order+1
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    d∂Ω = Measure(∂Ω, degree)
    ν = get_normal_vector(∂Ω)

    #h(x) =  ∇u(x) #Neumann on bounding box

    a((u,û),(v,v̂)) = ∫(∇(u)⋅∇(v))dΩ + ∫(∇(û)⋅∇(v̂))dΓ + γ*∫((v-v̂)*(u-û))dΓ

    l((v,v̂)) = ∫(f*v)dΩ + ∫(f̂*v̂)dΓ - ∫(g*v)dΓ

    # Build affine FE operator
    op = AffineFEOperator(a,l,X,Y)

    # Solve
    ls = LUSolver()
    solver = LinearFESolver(ls)
    wh = solve(solver, op)
    uh, ûh = wh

    eu = u - uh
    eΩ = sqrt(sum( ∫( eu*eu + ∇(eu)⋅∇(eu) )dΩ ))
    println("error_u: ", eΩ)
    eû = û - ûh
    eΓ = sqrt(sum( ∫( eû*eû + ∇(eû)⋅∇(eû) )dΓ ))
    println("error_û: ", eΓ)
    etot = sqrt(eΩ^2+eΓ^2)
    println("error_tot: ", etot)


    return etot, eΩ, eΓ
end


etot = Float64[]
eΩ = Float64[]
eΓ = Float64[]
γ = LinRange(1, 1000, 99)

for i in γ
    a, b, c = solver(i)
    push!(etot, a)
    push!(eΩ, b)
    push!(eΓ, c)

end

l = @layout [a{1.0w}
[grid(1,2)]]

plot(γ, [etot, eΩ, eΓ], legend = false, title = ["$i" for j in 1:1, i in ["e_tot", "||u-uₕ||", "||û-ûₕ||"]],
                  shape=:auto, xlabel="γ",ylabel="error norm", layout = l)
savefig("gamma_increase_2.pdf")
#savefig("gamma_error_2.pdf")
