# =========================================================================== #
# Compliant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("glouton_const2.jl")
include("glouton_descent.jl")
include("GRASP.jl")

# =========================================================================== #

# Loading a SPP instance
#println("\nLoading...")
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
# @show C
# #@show A

# grd, zg = greedy2(C, A)
# println("greedy = ",grd, ", z= ", zg)

# dea, v=deepest_d(grd)
# println("deapest = ", dea, ", value = ", v)

g = grasp(A, C, 0.7)
println(g, z(g, C))


# function utility_total(A, C)
#     #retourne la liste des utilités 
#         n = length(C) # Nombre de variables
        
#         # Initialiser une liste des utilités
#         U = zeros(n)
        
#         # Calculer l'utilité pour chaque variable (colonne de A)
#         for i in 1:n
#             # L'utilité est C[i] divisé par la somme des contraintes de la colonne i
#             U[i] = C[i] /sum(A[: , i]) 
#         end
        
#     return(U)
# end

#println(utility_total(A, C))

# =========================================================================== #

# Solving a SPP instance with GLPK
# println("\nSolving...")
# solverSelected = GLPK.Optimizer
# spp = setSPP(C, A)

# set_optimizer(spp, solverSelected)
# optimize!(spp)

# # Displaying the results
# println("z = ", objective_value(spp))
# print("x = "); println(value.(spp[:x]))

# # Solving a SPP instance with GLPK
# println("\nSolving...")
# solverSelected = GLPK.Optimizer
# spp = setSPP(C, A)

# set_optimizer(spp, solverSelected)
# optimize!(spp)

# # Displaying the results
# println("z = ", objective_value(spp))
# print("x = "); println(value.(spp[:x]))

# Collecting the names of instances to solve
# println("\nCollecting...")
# target = "Data"
# fnames = getfname(target)

# println("\nThat's all folks !")
