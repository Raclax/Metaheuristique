include("loadSPP.jl")
include("glouton_const2.jl")
include("glouton_descent.jl")
using InvertedIndices
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)



function grasp(A, C, alpha)
    S=zeros(Int, n)
    candidat_set = [i for i in 1:length(C)]
    utility_values = utility(A,C, length(C))
    umin = minimum(utility_values)
    ulim = umin + alpha*(maximum(utility_values) - umin)
    # cand=[]
    # for i=1:n
    #     push!(cand, i)
    # end
    while ulim != 0
        # println("Taille du set =", size(candidat_set))
        # print(ulim, "argmax = ", maximum(utility_values) - umin," ", alpha)
        
        RCL = [e for e in 1:length(utility_values) if utility_values[e] >= ulim]
        # println("RCL ", RCL)
        # Select an element e from the RCL at random
        e = rand(RCL)
        #candidat_set=filter!(a->a!=e, candidat_set)
        S[e] = 1
        utility_values[e]=0
        # for i in eachindex(candidat_set)
        #     for j in 1:m
        #         if(A[j,candidat_set[i]]+A[j,e]>1)
        #             push!(col2Remov,candidat_set[i])
        #         end
        #     end
        # end
        for i in 1:m
            if A[i, e] == 1
                for j in 1:n
                    if A[i, j] == 1 && j != e
                        utility_values[j]=0
                    end
                end
            end
        end
                        
        println("avant")
        println(candidat_set)
        umin = minimum(utility_values)
        ulim = umin + alpha*(maximum(utility_values) - umin)
        println("après")

    # while(!isempty(cand))
    #     utility_values = utility(A,C, length(cand))
    #     evali=copy(utility_values[cand])
    #     min=minimum(evali)
    #     max=maximum(evali)
    #     limit=min+alpha*(max-min)
    #     RCL=[]
    #     for i in eachindex(evali)
    #         if(evali[i]>=limit)
    #             push!(RCL,i)
    #         end
    #     end
    #     e=rand(collect(RCL))
    #     e2=cand[e]
    #     cand=filter!(a->a!=e, cand)
    #     del=[]
    #     for i in eachindex(cand)
    #         for j in 1:m
    #             if(A[j,cand[i]]+A[j,e2]>1)
    #                 push!(del,cand[i])
    #             end
    #         end
    #     end
    #     cand=setdiff(cand,del)
    #     x0[e2]=1
    # end
    # println("résultat pour GRASP ", z(C,x0))
    # return x0

    
        # Update the candidate set C and reevaluate the utility
        # C, A, utility_values = conflict_update(C, A, e, utility_values)
    end

    return S
end
