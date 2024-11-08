include("loadSPP.jl")
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)

C = Int.(C)
A = Int.(A)


function utility(A, C, n)

    U = zeros(n)
    for i in 1:n
        U[i] = C[i] /sum(A[: , i]) 
    end
    return U
end


function greedy2(C, A)
    m, n = size(A)

    u = sortperm(utility(A, C, n), rev=true)

    x = zeros(Int, n)

    for un in u
        #println("un = ", un)
        feasible = true
        for i in 1:m
            if A[i, un] == 1
                for j in 1  :n
                    if x[j] == 1 && A[i, j] == 1
                        feasible = false
                        break
                    end
                end
            end
            if !feasible
                break
            end
        end
        if feasible
            x[un] = 1
        end
    end
    return x, dot(x, C)
end




