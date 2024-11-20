# =========================================================================== #
# Compliant julia 1.x

# Using the following packages
#using JuMP, GLPK

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")

# =========================================================================== #

# Loading a SPP instance
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)

# =========================================================================== #

#Solving a SPP instance with GLPK
#println("\nSolving...")
solverSelected = GLPK.Optimizer
spp = setSPP(C, A)

set_optimizer(spp, solverSelected)
optimize!(spp)

# Displaying the results
println("z = ", objective_value(spp))
# print("x = "); println(value.(spp[:x]))


#Collecting the names of instances to solve
# println("\nCollecting...")
# target = "../Data"
# fnames = getfname(target)

#println("\nThat's all folks !")
