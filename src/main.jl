# =========================================================================== #
# Compliant julia 1.x

# Using the following packages
using GLPK

include("loadSPP.jl")
include("setSPP.jl")


# =========================================================================== #

#Solving a SPP instance with GLPK
function solveGLPK(fname)
    C, A = loadSPP(fname)


    solverSelected = GLPK.Optimizer
    spp = setSPP(C, A)

    set_optimizer(spp, solverSelected)
    optimize!(spp)

    # Displaying the results
    println("z = ", objective_value(spp))
    print("x = "); println(value.(spp[:x]))
end 
