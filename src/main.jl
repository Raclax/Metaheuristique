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

# =========================================================================== #

# Loading a SPP instance
println("\nLoading...")
fname = "Data/pb_500rnd0100.dat"
C, A = loadSPP(fname)
@show C
#@show A

grd = greedy(C, A)
println("greedy = ",grd, ", z= ", z(grd, C))

dea, v=deepest_d(grd)
println("deapest = ", dea, ", value = ", v)


# # Solving a SPP instance with GLPK
# println("\nSolving...")
# solverSelected = GLPK.Optimizer
# spp = setSPP(C, A)

# set_optimizer(spp, solverSelected)
# optimize!(spp)

# # Displaying the results
# println("z = ", objective_value(spp))
# print("x = "); println(value.(spp[:x]))

# =========================================================================== #

# Collecting the names of instances to solve
# println("\nCollecting...")
# target = "Data"
# fnames = getfname(target)

# println("\nThat's all folks !")
