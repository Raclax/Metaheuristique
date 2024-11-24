include("loadSPP.jl")
include("EI1.jl")
using LinearAlgebra
using InvertedIndices
using Distributions
using StatsBase


function grasp(fname, alpha=0.7, repeat=12)
    C, A = loadSPP(fname)
    m, n = size(A)
    best = zeros(Int, n) # solution
    best_valeur = 0 # valeur de la solution

    for i in 1:repeat
        S = grconst(A, C, alpha)
        v = z(S, C)
        S, v = deepest_d(S, A, C) # amélioration, peut être évitée pour améliorer le CPUt
        
        if v > best_valeur
            best = S
            best_valeur = v # devient la meilleure valeur si meilleure que la précédente
        end
    end
    return best, best_valeur
end

# Construction randomisée aléatoire
function grconst(A, C, alpha)
    m, n = size(A)
    S=zeros(Int, n) # solution
    utility_values = utility(A,C, length(C)) # utilités
    umin = minimum(utility_values)
    ulim = umin + alpha*(maximum(utility_values) - umin)

    while ulim > 0  #tant que toutes les utilités ne sont pas nulles
        RCL = [e for e in 1:length(utility_values) if utility_values[e] >= ulim] # liste des candidats restreints

        if isempty(RCL)
            break
        end
        e = rand(RCL)
        S[e] = 1 # on prend une colonne aléatoirement dans la RCL
        utility_values[e]=0 # et on met son utilité à 0 pour pas qu'elle ne soit reprise
        
        for i in 1:m
            if A[i, e] == 1
                for j in 1:n
                    if A[i, j] == 1 && j != e
                        utility_values[j]=0 # on met à 0 les utilités des colonnes qui ne peuvent plus être prises également
                    end
                end
            end
        end
                        
        umin = minimum(utility_values)
        ulim = umin + alpha*(maximum(utility_values) - umin) #on update les bornes

    end
    return S
end


function z(x, values)
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


function deepest_d(solution, A, C)
    best = solution
    best_valeur = z(best, C)
    ones = findall(x -> x == 1, solution)
    zeros = findall(x -> x == 0, solution)

    # en 2-1
    new=true
    while new
        new=false

        x, zx = echange_xx(2, 1, best, zeros, ones, A, C)
        if zx > best_valeur
            best = x
            best_valeur = zx
            ones = findall(x -> x == 1, best)
            zeros = findall(x -> x == 0, best)
            new=true
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


function echange_xx(k, p, solution, zeros, ones, A, C)

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


### REACTIVE GRASP


function ReactiveGrasp(fname, m=10, maxIteration=85, Nalpha=15)  
    C, A = loadSPP(fname)

    pk = fill(1 / m, m) 
    alpha_values = [(i - 1) / m for i in 1:m]
    n_alpha = zeros(Int, m)

    # Initialisation pour avoir un average non null
    z_avg = zeros(Float64, m)
    solu = []

    for i in 1:m
        c, z_avg[i] = grasp(fname, i / m)
        push!(solu, c)
    end

    z_best = maximum(z_avg)
    s_best = solu[argmax(z_avg)]
    z_worst = minimum(z_avg)
    best_alpha = alpha_values[argmax(z_avg)]

    # Construire la solution en utilisant l'algorithme de construction
    solution, valactuelle = grasp(fname, alpha_values[1])

    for iter in 1:maxIteration
        # Choisir une valeur de alpha en fonction des probabilités pk
        # println(Weights(pk))
        alpha_idx = sample(1:m, Weights(pk))
        alpha = alpha_values[alpha_idx]
    
        solution, valactuelle = grasp(fname, alpha)

        # Mettre à jour z_best, z_worst, moyenne
        if valactuelle > z_best
            z_best = valactuelle
            s_best = solution
            best_alpha = alpha  # Mettre à jour best_alpha avec la nouvelle meilleure valeur
        end

        z_worst = min(z_worst, valactuelle)
        
        n_alpha[alpha_idx] += 1

        # moyenne d'un alpha = somme des valeurs de graps obtenues aplha / nombre d'essaies du alpha,
        # on a  la somme = ancienne moyenne * nombre d'essaies réalisé avant la nouvelle itération + nouveau resultat,
        # il suffit ensuite juste de diviser par le nombre actuel d'itérations.
        z_avg[alpha_idx] = ((n_alpha[alpha_idx] - 1) * z_avg[alpha_idx] + valactuelle) / n_alpha[alpha_idx]
        
        if iter % Nalpha == 0 
            # on calcul periodiquement selon Nalpha
            # Recalculer les probabilités pk toutes les N_alpha itérations
            qk = [(z_avg[k] - z_worst) / (z_best - z_worst) for k in 1:m]
            pk = [qk[k] / sum(qk) for k in 1:m]
        end

    end
    
    return s_best, z_best, best_alpha # Retourne la meilleure solution trouvée


end

function experimentationSPP()
    println("didactic", genetic_algorithm("../Data/didactic.dat"))
    println("pb_100rnd0100", genetic_algorithm("../Data/pb_100rnd0100.dat"))
    println("pb_100rnd0300", genetic_algorithm("../Data/pb_100rnd0300.dat"))
    println("pb_200rnd0100", genetic_algorithm("../Data/pb_200rnd0100.dat"))
    println("pb_200rnd0500", genetic_algorithm("../Data/pb_200rnd0500.dat"))
    println("pb_200rnd1600", genetic_algorithm("../Data/pb_200rnd1600.dat"))
    println("pb_500rnd0100", genetic_algorithm("../Data/pb_500rnd0100.dat"))
    println("pb_500rnd1700", genetic_algorithm("../Data/pb_500rnd1700.dat"))
    println("pb_1000rnd0100", genetic_algorithm("../Data/pb_1000rnd0100.dat"))
    println("pb_1000rnd0200", genetic_algorithm("../Data/pb_1000rnd0200.dat"))
    println("pb_2000rnd0100", genetic_algorithm("../Data/pb_2000rnd0100.dat"))
end