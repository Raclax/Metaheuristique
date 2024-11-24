# =========================================================================== #
# Compliant julia 1.x

# Using the following packages
using GLPK
using LinearAlgebra

include("loadSPP.jl")


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

function setSPP(C, A)
    m, n = size(A)
    spp = Model()
    @variable(spp, x[1:n], Bin)
    @objective(spp, Max, dot(C, x))
    @constraint(spp , cte[i=1:m], sum(A[i,j] * x[j] for j=1:n) <= 1)
    return spp
  end