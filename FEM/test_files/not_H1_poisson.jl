using Gridap
using GridapGmsh

# Here, we define the manufactured functions, we find the correspondng source terms for our manufactured solutions
# only the source terms f, f̂ are used to find our solution (u, û) while the gradients will be used later to calculate the
u(x) = sin(abs(x[2]-0.5))-2*abs(x[2]-0.5)
Δu(x) = -sin(abs(x[2]-0.5))
f(x) = -Δu(x)


function solver(mesh ;write::Bool=false)

    model = GmshDiscreteModel("/Users/martinkristiansen/Desktop/Master thesis/Meshes/mesh_"*string(mesh)*".msh")
    #model = simplexify(model)

    order = 1 #Choose order of local interpolation, order 1 => piecewise linear, order 2 => piecewise quadratic etc..

    #Define local basis functions
    reffe_u = ReferenceFE(lagrangian, Float64, order)

    #where should our solution live
    Ω = Triangulation(model)

    V = TestFESpace(Ω,reffe_u, dirichlet_tags = ["boundary"], conformity=:H1)

    # Define global trial space
    U = TrialFESpace(V, u) #include u to enforce dirichlet tradtionally

    #integration measure
    degree = order+1
    dΩ = Measure(Ω,degree)

    a(u,v) = ∫(∇(u)⋅∇(v))dΩ

    l(v) = ∫(f*v)dΩ  #+∫((h⋅ν)*v)d∂Ω

    # Build affine FE operator
    op = AffineFEOperator(a,l,U,V)

    # Solve
    ls = LUSolver()
    solver = LinearFESolver(ls)
    uh = solve(solver, op)
    eu = u - uh
    etot = sqrt(sum( ∫( eu*eu + ∇(eu)⋅∇(eu) )dΩ ))
    println(etot)

    if write
        writevtk(Ω,"/Users/martinkristiansen/Desktop/Master thesis/Coupled-system/vtu_files/error_u",cellfields=["eu" => eu])
        writevtk(Ω, "/Users/martinkristiansen/Desktop/Master thesis/Coupled-system/vtu_files/omega", cellfields=["uh" => uh])
        writevtk(Ω,"/Users/martinkristiansen/Desktop/Master thesis/Coupled-system/vtu_files/u",cellfields=["u" => u])
    end

    return etot
end


using Plots
function conv_test(figname)

    e = Float64[]
    for i in [1,2,3,4,5]
        etot = solver(i, write = true)
        push!(e, etot)
        println(e)
    end

    h = [0.1767, 0.0914, 0.048, 0.0246, 0.0124] #these are extracted through inspection

    function slope(hs,errors)
      x = log10.(hs)
      y = log10.(errors)
      linreg = hcat(fill!(similar(x), 1), x) \ y
      linreg[2]
    end

    println("slope: ",slope(h, e))

    plot(h, [e, h], xaxis=:log, yaxis=:log, label=["e_tot" "h^0.5"],
                      shape=:auto, xlabel="h",ylabel="error norm")

    savefig(figname)
end


conv_test("h1_possion.pdf")
#solver(1, write = true)
