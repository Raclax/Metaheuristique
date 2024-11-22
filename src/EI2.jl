include("loadSPP.jl")
include("EI1.jl")
using LinearAlgebra
using InvertedIndices
using Distributions
using StatsBase

using Statistics
using Dates

fname = "../Data/pb_200rnd0500.dat"
C, A = loadSPP(fname)
m, n = size(A)


function grasp(A, C, alpha, repeat=12)
    best = zeros(Int, n) # solution
    best_valeur = 0 # valeur de la solution

    for i in 1:repeat
        S = grconst(A, C, alpha)
        v = z(S, C)
        #S, v = deepest_d(S) # construction et amélioration d'une solution
        
        if v > best_valeur
            best = S
            best_valeur = v # devient la meilleure valeur si meilleure que la précédente
        end
    end
    return best, best_valeur
end

# Construction randomisée aléatoire
function grconst(A, C, alpha)
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


function deepest_d(solution)
    best = solution
    best_valeur = z(best, C)
    ones = findall(x -> x == 1, solution)
    zeros = findall(x -> x == 0, solution)

    # en 2-1
    new=true
    while new
        new=false

        x, zx = echange_xx(2, 1, best, zeros, ones)
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

        x, zx = echange_xx(1, 1, best, zeros, ones)
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

        x, zx = echange_xx(0, 1, best, zeros, ones)
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


function echange_xx(k, p, solution, zeros, ones)

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


function ReactiveGrasp(A, C, m, maxIteration, Nalpha)  
    
    pk = fill(1 / m, m) 
    alpha_values = [(i - 1) / m for i in 1:m]
    n_alpha = zeros(Int, m)

    # Initialisation pour avoir un average non null
    z_avg = zeros(Float64, m)
    solu = []

    for i in 1:m
        c, z_avg[i] = grasp(A, C, i / m)
        push!(solu, c)
    end

    z_best = maximum(z_avg)
    s_best = solu[argmax(z_avg)]
    z_worst = minimum(z_avg)
    best_alpha = alpha_values[argmax(z_avg)]

    # Construire la solution en utilisant l'algorithme de construction
    solution, valactuelle = grasp(A, C, alpha_values[1])

    for iter in 1:maxIteration
        # Choisir une valeur de alpha en fonction des probabilités pk
        alpha_idx = sample(1:m, Weights(pk))
        alpha = alpha_values[alpha_idx]
    
        solution, valactuelle = grasp(A, C, alpha)

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

#s,v, a = ReactiveGrasp(A, C, 10,85,15)

#grasp(A, C, 5, 10)


function run_reactivegrasp_multiple_times(A, C, population_size, iterations, Nalpha)
    valactuelle_values = []
    #running_times = []
    alpha_values = []

    for _ in 1:10
        #start_time = now()
        _, valactuelle, alpha = ReactiveGrasp(A, C, population_size, iterations, Nalpha)
        #end_time = now()
        
        push!(valactuelle_values, valactuelle)
        #push!(running_times, end_time - start_time)
        push!(alpha_values, alpha)
    end

    max_valactuelle = maximum(valactuelle_values)
    min_valactuelle = minimum(valactuelle_values)
    mean_valactuelle = mean(valactuelle_values)
    #mean_running_time = mean(running_times)
    mean_alpha = mean(alpha_values)


    return max_valactuelle, min_valactuelle, mean_valactuelle, mean_alpha
end

run_reactivegrasp_multiple_times(A, C, 10, 85, 15)