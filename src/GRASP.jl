include("loadSPP.jl")
include("glouton_const2.jl")
include("glouton_descent.jl")
using InvertedIndices
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)


function grasp(A, C, alpha, repeat=10)
    best = zeros(Int, n)
    best_valeur = 0
    for i in 1:repeat
        S = grconst(A, C, alpha)
        #v=z(S,C)
        S, v = deepest_d(S)
        if v > best_valeur
            best = S
            best_valeur = v
        end
    end
    return best, best_valeur
end


function grconst(A, C, alpha)
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
        
        for i in 1:m
            if A[i, e] == 1
                for j in 1:n
                    if A[i, j] == 1 && j != e
                        utility_values[j]=0
                    end
                end
            end
        end
                        
        # println("avant")
        # println(candidat_set)
        umin = minimum(utility_values)
        ulim = umin + alpha*(maximum(utility_values) - umin)
        #println("après")

    end

    return S
end


function z(x, values)
    # println("Dimensions de x: ", x)
    # println("Dimensions de values: ", size(values))
    return dot(x, values)
end

function valide(sol, A)
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
    ones = findall(x -> x == 1, solution)
    zeros = findall(x -> x == 0, solution)
    # en 2-1
    new=true
    while new
        #println("Début-recherche dans un voisinage avec mouvement 2-1")
        new=false

        x, zx = echange_xx(2, 1, best, zeros, ones)
        if zx > best_valeur
            best = x
            best_valeur = zx
            ones = findall(x -> x == 1, best)
            zeros = findall(x -> x == 0, best)

            #println("nouveau meilleur = ", best_valeur)
            new=true
        end

    #println("Fin--recherche dans un voisinage avce mouvement 2-1")
    end

    #println("______________")
           

    # en 1-1
    new=true

    while new
        #println("Début-recherche dans un voisinage avec mouvement 1-1")
        new=false

        x, zx = echange_xx(1, 1, best, zeros, ones)
        if zx > best_valeur
            best = x
            best_valeur = zx
            ones = findall(x -> x == 1, best)
            zeros = findall(x -> x == 0, best)

            #println("nouveau meilleur = ", best_valeur)
            new=true
        end

        
    #println("Fin--recherche dans un voisinage avce mouvement 1-1")           
    end

    #println("______________")
        

    # en 0-1
    new=true

    while new
        #println("Début-recherche dans un voisinage avec mouvement 0-1")
        new=false
        voisins, idx= echange_xx(0, 1, best, zeros, ones)
        # println("Il y en a ",length(voisins), " voisins")

        x, zx = echange_xx(0, 1, best)
        if zx > best_valeur
            best = x
            best_valeur = zx
            ones = findall(x -> x == 1, best)
            zeros = findall(x -> x == 0, best)

            #println("nouveau meilleur = ", best_valeur)
            new=true
        end
    #println("Fin--recherche dans un voisinage avce mouvement 0-1")
    end

    #println("______________")
        

    
    return best, best_valeur
end 

function echange_xx(k, p, solution, zeros, ones)
    # l = Vector{Vector{Int}}() 
    # idx = Vector{Tuple{Vector{Int}, Vector{Int}}}()
    valactuelle = z(solution, C)

    if k == 0
        for j in zeros
            voisin = copy(solution)
            voisin[j] = 1 

            voisin_value = z(voisin, C)
            if voisin_value > valactuelle
                if valide(voisin, A)
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

            
                voisin_value = z(voisin, C)
                if voisin_value > valactuelle
                    if valide(voisin, A)
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

                    voisin_value = z(voisin, C)
                    if voisin_value > valactuelle
                        if valide(voisin, A)
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