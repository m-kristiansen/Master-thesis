using Gridap
using GridapGmsh

# Here, we define the manufactured functions, we find the correspondng source terms for our manufactured solutions
# only the source terms f, f̂ are used to find our solution (u, û) while the gradients will be used later to calculate the
γ = 100


model = GmshDiscreteModel("/Users/martinkristiansen/Desktop/Master thesis/Meshes/straight_line.msh")
#model = simplexify(model)

order = 2 #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

#Define local basis functions
reffe_u = ReferenceFE(lagrangian, Float64, order)
reffe_û = ReferenceFE(lagrangian, Float64, order)

#where should our solution live
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model, tags = ["Graph"])
∂Ω = BoundaryTriangulation(model, tags = ["boundary"])

V = TestFESpace(Ω,reffe_u, conformity=:H1) #we must add a constraint for neumann BC's
V̂ = TestFESpace(Γ,reffe_û, dirichlet_tags = ["starting_point", "end_point"], conformity=:H1)

Y = MultiFieldFESpace([V,V̂])


# Define global trial space
U = TrialFESpace(V) #include u to enforce dirichlet tradtionally
Û = TrialFESpace(V̂, [1, 0]) # 1 dirichlet on start of graph
X = MultiFieldFESpace([U,Û])

#integration measure
degree = order+1
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
d∂Ω = Measure(∂Ω, degree)
ν = get_normal_vector(∂Ω)

h(x) = 0 #Neumann on bounding box
f(x) = 0
f̂(x) = 0

a((u,û),(v,v̂)) = ∫(∇(u)⋅∇(v))dΩ + γ*∫((v-v̂)*(u-û))dΓ + ∫(∇(û)⋅∇(v̂))dΓ

l((v,v̂)) = ∫(f*v)dΩ + ∫(f̂*v̂)dΓ +∫(h*v)d∂Ω

# Build affine FE operator
op = AffineFEOperator(a,l,X,Y)

# Solve
ls = LUSolver()
solver = LinearFESolver(ls)
wh = solve(solver, op)
uh, ûh = wh


writevtk(Ω, "../vtu_files/omega", cellfields=["uh" => uh])
writevtk(Γ, "../vtu_files/gamma", cellfields=["uhat" => ûh])
