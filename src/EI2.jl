include("loadSPP.jl")
include("EI1.jl")
using LinearAlgebra
using InvertedIndices
using Distributions
using StatsBase
using Random
using Base.Threads

fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)


function grasp(A, C, alpha, repeat=10)
    best = zeros(Int, n)
    best_valeur = 0

    for i in 1:repeat
        S = grconst(A, C, alpha)
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

    while ulim > 0

        
        RCL = [e for e in 1:length(utility_values) if utility_values[e] >= ulim]

        if isempty(RCL)
            break
        end
        e = rand(RCL)
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
                        
        umin = minimum(utility_values)
        ulim = umin + alpha*(maximum(utility_values) - umin)

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

    z_avg = zeros(Float64, m)
    for i in 1:m
        z_avg[i] = grasp(A, C, i / m)[2]
    end

    z_best = maximum(z_avg)
    z_worst = minimum(z_avg)

    solution, valactuelle = grasp(A, C, alpha_values[1])

    for iter in 1:maxIteration
        alpha_idx = sample(1:m, Weights(pk))
        alpha = alpha_values[alpha_idx]
    
        solution, valactuelle = grasp(A, C, alpha)

        n_alpha[alpha_idx] += 1
        z_avg[alpha_idx] = ((n_alpha[alpha_idx] - 1) * z_avg[alpha_idx] + valactuelle) / n_alpha[alpha_idx]
    
        z_best = max(z_best, valactuelle)
        z_worst = min(z_worst, valactuelle)
        for i in 1:m
            pk[i] = exp(-z_avg[i] / (z_best - z_worst + 1e-6))
        end
        pk /= sum(pk) 
    end
    
    return solution, valactuelle

end

s,v = ReactiveGrasp(A, C, 10,100,10)