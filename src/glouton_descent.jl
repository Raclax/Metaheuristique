include("loadSPP.jl")
include("glouton_const2.jl")
fname = "Data/pb_500rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)

function z(x, values)
    return dot(x, values)
end

function valide(sol, A)
    colonnes = findall(x -> x == 1, sol)
    cidx = [colonne[1] for colonne in colonnes]

    for row in 1:size(A,1)

        if sum(A[row, cidx])>1
            return false
        end

    end

    return true
end

function valides(sol, A)

    for i in 1:m

        if sum(A[i, j] * sol[j] for j in 1:n)>1
            return false
        end

    end

    return true 
end

function deepest_d(solution)
    best = solution
    best_valeur = z(best, C)

    
    # en 2-1
    new=true
    while new
        println("Début-recherche dans un voisinage avec mouvement 2-1")
        new=false

        voisins, idx, t = echange_xx(2, 1, best)
        println("Il y en a ",length(voisins), ", et ", t, " testés")

        for voisin in voisins

            if z(voisin, C) > best_valeur
                println("ancien meilleur = ", best_valeur)
                best = voisin
                best_valeur = z(voisin, C)

                println("nouveau meilleur = ", best_valeur)
                new=true
            end

        end

    println("Fin--recherche dans un voisinage avce mouvement 2-1")
    end

    println("______________")
           

    # en 1-1
    new=true

    while new
        println("Début-recherche dans un voisinage avec mouvement 1-1")
        new=false

        voisins, idx, t= echange_xx(1, 1, best)
        println("Il y en a ",length(voisins), ", et ", t, " testés")

        for voisin in voisins

            if z(voisin, C) > best_valeur
                println("ancien meilleur = ", best_valeur)
                best = voisin
                best_valeur = z(voisin, C)
                println("nouveau meilleur = ", best_valeur)
                new=true
            end

        end

    println("Fin--recherche dans un voisinage avce mouvement 1-1")           
    end

    println("______________")
        

    # en 0-1
    new=true

    while new
        println("Début-recherche dans un voisinage avec mouvement 0-1")
        new=false

        voisins, idx, t= echange_xx(0, 1, best)
        println("Il y en a ",length(voisins), ", et ", t, " testés")

        for voisin in voisins

            if z(voisin, C) > best_valeur
                println("ancien meilleur = ", best_valeur)
                best = voisin
                best_valeur = z(voisin, C)

                println("nouveau meilleur = ", best_valeur)
                new=true
            end

        end

    println("Fin--recherche dans un voisinage avce mouvement 0-1")
    end

    println("______________")
        

    
    return best, best_valeur
end 



function echange_xx(k, p, solution)
    l = [] 
    idx = []
    ones = findall(x -> x == 1, solution)
    zeros = findall(x -> x == 0, solution)
    trys = 0

    if k == 0

        for j in zeros
            voisin = copy(solution)
            voisin[j] = 1  
            trys+=1

            if valide(voisin, A) && z(voisin, C) >= z(solution, C)
                push!(l, voisin)
                push!(idx, ([], [j]))
            end
        end
        
    elseif k == 1

        for i in ones
            for j in zeros
                voisin = copy(solution)
                voisin[i] = 0  
                voisin[j] = 1 
                trys+=1

                if valide(voisin, A) && z(voisin, C) >= z(solution, C)
                    push!(l, voisin)
                    push!(idx, ([i], [j]))
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
                    trys+=1

                    if valide(voisin, A) && z(voisin, C) >= z(solution, C)
                        push!(l, voisin)
                        push!(idx, ([ones[i], ones[k_idx]], [j]))
                    end
                end
            end
        end

    end
    return l, idx, trys
end
