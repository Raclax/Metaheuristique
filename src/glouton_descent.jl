include("loadSPP.jl")
#include("glouton_const2.jl")
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)


function utility(A, C, n)

    U = zeros(n)
    for i in 1:n
        U[i] = C[i] /sum(A[: , i]) 
    end
    return U
end


function greedy2(C, A)
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



function z(x, values)
    return dot(x, values)
end

function valide(sol)
    colonnes = findall(x -> x == 1, sol)
    cidx = [colonne[1] for colonne in colonnes]

    if any(sum(A[:, cidx], dims=2) .> 1)
        return false
    end
    return true
end


function deepest_d(solution)
    best = solution
    best_valeur = z(best, C)
    # en 2-1
    new=true
    while new
        #println("Début-recherche dans un voisinage avec mouvement 2-1")
        new=false

        x, zx = echange_xx(2, 1, best)
        if zx > best_valeur
            best = x
            best_valeur = zx

            #println("nouveau meilleur = ", best_valeur)
            new=true
        end

    #println("Fin--recherche dans un voisinage avce mouvement 2-1")
    end

    println("______________")
           

    # en 1-1
    new=true

    while new
        #println("Début-recherche dans un voisinage avec mouvement 1-1")
        new=false

        x, zx = echange_xx(1, 1, best)
        if zx > best_valeur
            best = x
            best_valeur = zx

            #println("nouveau meilleur = ", best_valeur)
            new=true
        end

        
    #println("Fin--recherche dans un voisinage avce mouvement 1-1")           
    end

    println("______________")
        

    # en 0-1
    new=true

    while new
        #println("Début-recherche dans un voisinage avec mouvement 0-1")
        new=false
        voisins, idx= echange_xx(0, 1, best)
        # println("Il y en a ",length(voisins), " voisins")

        x, zx = echange_xx(0, 1, best)
        if zx > best_valeur
            best = x
            best_valeur = zx

            #println("nouveau meilleur = ", best_valeur)
            new=true
        end
    #println("Fin--recherche dans un voisinage avce mouvement 0-1")
    end

    println("______________")
        

    
    return best, best_valeur
end 



function echange_xx(k, p, solution)
    # l = Vector{Vector{Int}}() 
    # idx = Vector{Tuple{Vector{Int}, Vector{Int}}}()
    ones = findall(x -> x == 1, solution)
    zeros = findall(x -> x == 0, solution)
    valactuelle = z(solution, C)

    if k == 0
        for j in zeros
            voisin = copy(solution)
            voisin[j] = 1  

            if valide(voisin) 
                voisin_value = z(voisin, C)
                if voisin_value > valactuelle
                    solution = voisin
                    valactuelle = voisin_value
                end
            end
        end
        
    elseif k == 1

        for i in ones
            for j in zeros
                voisin = copy(solution)
                voisin[i] = 0  
                voisin[j] = 1 

                if valide(voisin)
                    voisin_value = z(voisin, C)
                    if voisin_value > valactuelle
                        solution = voisin
                        valactuelle = voisin_value
                    end
                end
            end      
        end

    else

        for i in 1:length(ones)-1
            for k_idx in i+1:length(ones)
                for j in zeros
                    voisin = copy(solution)
                    voisin[ones[i]] = 0 
                    voisin[ones[k_idx]] = 0  
                    voisin[j] = 1  

                    if valide(voisin)
                        voisin_value = z(voisin, C)
                        if voisin_value >= valactuelle
                            solution = voisin
                            valactuelle = voisin_value
                        end
                    end
                end
            end
        end

    end
    return solution, valactuelle
end
