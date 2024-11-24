include("loadSPP.jl")
using LinearAlgebra

# fname = "../Data/pb_100rnd0300.dat"
# C, A = loadSPP(fname)
# m, n = size(A)

function utility(A, C, n)
    U = zeros(n)
    for i in 1:n
        U[i] = C[i] /sum(A[: , i]) 
    end
    return U
end


# Construction avec algorithme glouton
function greedy(fname) 
    C, A = loadSPP(fname)
    m, n = size(A)
    utilites = sortperm(utility(A, C, n), rev=true) # trier les utilités
    x = zeros(Int, n) # solution

    for best_u in utilites 
        feasible = true
        for i in 1:m
            if A[i, best_u] == 1 # si 1 dans la colonne
                for j in 1:n
                    if x[j] == 1 && A[i, j] == 1 # si 1 dans la ligne + dans la solution, alors colonne non validée
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
            x[best_u] = 1 # si non, colonne validée
        end
    end
    return x, dot(x, C)
end


# Calcul de la valeur de la fonction
function z(x, values)
    return dot(x, values)
end

# Fonction de validation d'une solution
function valide(sol, A)
    colonnes = findall(x -> x == 1, sol)
    cidx = [colonne[1] for colonne in colonnes]

    if any(sum(A[:, cidx], dims=2) .> 1)
        return false
    end
    return true
end


# Fonction de recherche locale en decente profonde
function resoudreSPP(fname)
    C, A = loadSPP(fname)
    m, n = size(A)
    best = greedy(fname)[1]
    best_valeur = z(best, C)
    ones = findall(x -> x == 1, best)
    zeros = findall(x -> x == 0, best)


    # en 2-1
    new=true # variable pour savoir si on a trouvé une meilleure solution
    while new
        new=false

        x, zx = echange_xx(2, 1, best, zeros, ones, A, C)
        if zx > best_valeur
            best = x
            best_valeur = zx
            ones = findall(x -> x == 1, best)
            zeros = findall(x -> x == 0, best)

            new=true # si on en a trouvé une, alors on recommence ce kp
        end
    end

    # en 1-1
    new=true

    while new
        new=false

        x, zx = echange_xx(1, 1, best, zeros, ones, A, C)
        if zx > best_valeur
            best = x
            best_valeur = zx
            ones = findall(x -> x == 1, best)
            zeros = findall(x -> x == 0, best)

            new=true
        end
    end        

    # en 0-1
    new=true

    while new
        new=false

        x, zx = echange_xx(0, 1, best, zeros, ones, A, C)
        if zx > best_valeur
            best = x
            best_valeur = zx
            ones = findall(x -> x == 1, best)
            zeros = findall(x -> x == 0, best)

            new=true
        end
    end
    
    return best, best_valeur
end 


#Fonction unique pour 0-1, 1-1 et 2-1
function echange_xx(k, p, solution, zeros, ones, A, C)

    valactuelle = z(solution, C)

    if k == 0 # 0-1
        for j in zeros
            voisin = copy(solution)
            voisin[j] = 1  # Change un 0 en 1

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
                voisin[j] = 1 # Change un 0 en 1 et un 1 en 0

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
                    voisin[j] = 1  # Change un 0 en 1 et deux 1 en 0

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



function experimentationSPP()
    println("didactic", resoudreSPP("../Data/didactic.dat"))
    println("pb_100rnd0100", resoudreSPP("../Data/pb_100rnd0100.dat"))
    println("pb_100rnd0300", resoudreSPP("../Data/pb_100rnd0300.dat"))
    println("pb_200rnd0100", resoudreSPP("../Data/pb_200rnd0100.dat"))
    println("pb_200rnd0500", resoudreSPP("../Data/pb_200rnd0500.dat"))
    println("pb_500rnd0100", resoudreSPP("../Data/pb_500rnd0100.dat"))
    println("pb_500rnd0100", resoudreSPP("../Data/pb_500rnd0100.dat"))
    println("pb_500rnd1700", resoudreSPP("../Data/pb_500rnd1700.dat"))
    println("pb_1000rnd0100", resoudreSPP("../Data/pb_1000rnd0100.dat"))
    println("pb_1000rnd0200", resoudreSPP("../Data/pb_1000rnd0200.dat"))
    println("pb_2000rnd0100", resoudreSPP("../Data/pb_2000rnd0100.dat"))
end