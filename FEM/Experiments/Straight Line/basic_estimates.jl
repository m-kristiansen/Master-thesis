using Gridap
using GridapGmsh

# Here, we define the manufactured functions, we find the correspondng source terms for our manufactured solutions
# only the source terms f, f̂ are used to find our solution (u, û) while the gradients will be used later to calculate the

# γ = 1
#
# u(x) = x[1]-x[2]
# ∇u(x) = VectorValue(1,-1)
# Δu(x) = 0
#
# û(x) = x[1]
# ∇û(x) = VectorValue(1,0)
# Δû(x) = 0
#
# f(x) = -Δu(x)
# g(x) = -γ*(u(x)-û(x))
# f̂(x) = -Δû(x)+γ*(û(x)-u(x))


# γ = 1
#
# u(x) = x[1]-x[2]
# ∇u(x) = VectorValue(1,-1)
# Δu(x) = 0
#
# û(x) = x[1]^2
# ∇û(x) = VectorValue(2*x[1],0)
# Δû(x) = 2
#
# f(x) = -Δu(x)
# g(x) = -γ*(u(x)-û(x))
# f̂(x) = -Δû(x)+γ*(û(x)-u(x))

γ = 0

u(x) = cos(x[1])+sin(x[2])
∇u(x) = VectorValue(-sin(x[1]),cos(x[2]))
Δu(x) = -cos(x[1])-sin(x[2])

û(x) = sin(x[1])
∇û(x) = VectorValue(cos(x[1]),0)
Δû(x) = -sin(x[1])

f(x) = -Δu(x)
g(x) = -γ*(u(x)-û(x))
f̂(x) = -Δû(x)+γ*(û(x)-u(x))


import Gridap: ∇
∇(::typeof(u)) = ∇u #get grad u from here


model = GmshDiscreteModel("/Users/martinkristiansen/Desktop/Master thesis/Meshes/mesh_1.msh")
#model = simplexify(model)

order = 1 #Choose order of local interpolation, order 0 => piecewise linear, order 1 => hat function etc..

#Define local basis functions
reffe_u = ReferenceFE(lagrangian, Float64, order)
reffe_û = ReferenceFE(lagrangian, Float64, order)

dirichlet_tags = ["top", "bottom", "bottom_left", "bottom_right",
"top_left", "top_right", "iface_left", "iface_right", "ll", "lr", "ul", "ur"]


#where should our solution live
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model, tags = ["interface", "iface_left", "iface_right"])
∂Ω = BoundaryTriangulation(model, tags = dirichlet_tags)


V = TestFESpace(Ω,reffe_u, conformity=:H1, dirichlet_tags = dirichlet_tags ) #we must add a constraint for neumann BC's
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
L2 = sqrt(sum( ∫( eu*eu )dΩ ))
println("error_L2: ", L2)

euh1 = sqrt(sum( ∫( eu*eu + ∇(eu)⋅∇(eu) )dΩ ))
println("error_u: ", euh1)
eû = û - ûh
eûh1 = sqrt(sum( ∫( eû*eû + ∇(eû)⋅∇(eû) )dΓ ))
println("error_û: ", eûh1)
etot = sqrt(euh1^2+eûh1^2)
println("error_tot: ", etot)


writevtk(Ω,"../vtu_files/error_u",cellfields=["eu" => eu])
writevtk(Γ,"../vtu_files/error_û",cellfields=["euhat" => eû])

writevtk(Ω, "../vtu_files/f₁", cellfields = ["f" => f] )
writevtk(Γ, "../vtu_files/f₂", cellfields = ["g" => g] )


writevtk(Ω, "../vtu_files/omega", cellfields=["uh" => uh])
writevtk(Γ, "../vtu_files/gamma", cellfields=["uhat" => ûh])

writevtk(Γ,"../vtu_files/normal",cellfields=["ν" => ν])
