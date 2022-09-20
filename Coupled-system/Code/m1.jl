using Gridap
using GridapGmsh

# Here, we define the manufactured functions, we find the correspondng source terms for our manufactured solutions
# only the source terms f, f̂ are used to find our solution (u, û) while the gradients will be used later to calculate the
γ = 0

u(x) = x[1]^2+x[2]^2
∇u(x) = VectorValue(2*x[1], 2*x[2])
Δu(x) = 4
f(x) = -Δu(x)+γ*(u(x)-û(x))

û(x) = x[1]+x[2]
∇û(x) = VectorValue(2*x[1], 2*x[2])
Δû(x) = 4
f̂(x) = -Δû(x)+γ*(û(x)-u(x))


import Gridap: ∇
∇(::typeof(u)) = ∇u #get grad u from here

model = GmshDiscreteModel("straight_line.msh")
#model = simplexify(model)

order = 2 #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

#Define local basis functions
reffe_u = ReferenceFE(lagrangian, Float64, order)
reffe_û = ReferenceFE(lagrangian, Float64, order)

#where should our solution live
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model, tags = ["Graph"])
∂Ω = BoundaryTriangulation(model, tags = ["boundary"])

V = TestFESpace(Ω,reffe_u, dirichlet_tags = ["boundary", "starting_point"], conformity=:H1)
V̂ = TestFESpace(Γ,reffe_û,dirichlet_tags = ["starting_point"], conformity=:H1)

Y = MultiFieldFESpace([V,V̂])


# Define global trial space
U = TrialFESpace(V, u) #include u to enforce dirichlet tradtionally
Û = TrialFESpace(V̂, û) # 1 dirichlet on start of graph
X = MultiFieldFESpace([U,Û])

#integration measure
degree = order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
d∂Ω = Measure(∂Ω, degree)
ν = get_normal_vector(∂Ω)

h(x) = ∇u(x) #Neumann on bounding box


a((u,û),(v,v̂)) = ∫(∇(u)⋅∇(v))dΩ + γ*∫((v-v̂)*(u-û))dΓ + ∫(∇(û)⋅∇(v̂))dΓ

l((v,v̂)) = ∫(f*v)dΩ + ∫(f̂*v̂)dΓ #+∫((h⋅ν)*v)d∂Ω

# Build affine FE operator
op = AffineFEOperator(a,l,X,Y)

# Solve
ls = LUSolver()
solver = LinearFESolver(ls)
wh = solve(solver, op)
uh, ûh = wh
eu = u - uh
println("error_u: ", sqrt(sum( ∫( eu*eu )dΩ )))
eû = û - ûh
println("error_û: ", sqrt(sum( ∫( eû*eû )dΓ )))

writevtk(Ω,"error_u",cellfields=["eu" => eu])
writevtk(Γ,"error_û",cellfields=["euhat" => eû])
writevtk(Ω, "omega", cellfields=["uh" => uh])
writevtk(Γ, "gamma", cellfields=["uhat" => ûh])

writevtk(∂Ω,"normal",cellfields=["ν" => ν])
