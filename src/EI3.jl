include("loadSPP.jl")
include("EI2.jl")
using LinearAlgebra
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)


function initialize_population(population_size, chromosome_length)
    return [rand(0:1, chromosome_length) for _ in 1:population_size]
end

# ----------------------------------------------------------------------------------------------------

function evaluate_population(population)
    return [fitness(chromosome) for chromosome in population]
end

# ----------------------------------------------------------------------------------------------------

function fitness(chromosome)
    return dot(chromosome, C) 
end

# ----------------------------------------------------------------------------------------------------

function selection(population)
    cand2 = population[rand(1:length(population))]
    cand1 = population[rand(1:length(population))]

    if fitness(cand2) > fitness(cand1)
        cand1 = cand2
    end

    filter!(elem -> elem[1] != cand1, population)

    return cand1  
end 
# ----------------------------------------------------------------------------------------------------

function crossover(p1, p2, crossover_rate)
    enfantA=copy(p1)
    enfantB=copy(p2)
    
    if rand() < crossover_rate
        crossover_point = rand(1:length(p1))
        enfantA = vcat(p1[1:crossover_point], p2[crossover_point+1:end])
        enfantB = vcat(p2[1:crossover_point], p1[crossover_point+1:end])
    end
    
    return enfantA, enfantB
end

# ----------------------------------------------------------------------------------------------------

function mutation(chromosome, mutation_rate)
    for i in eachindex(chromosome)
        if rand() < mutation_rate
            if chromosome[i] == 0
                chromosome[i] = 1
            else
                chromosome[i] = 0
            end
        end
    end
    return repare(chromosome)
end

# ----------------------------------------------------------------------------------------------------

function repare(chromosome)

    utilite = map(i -> C[i] / sum(A[:, i]), 1:n)
    while !valide(chromosome)
        conflict_indices = findall(x -> x > 1, A * chromosome)
        
        for idx in conflict_indices
            conflicting_elements = findall(x -> x == 1, A[idx, :])
            min_value_index = argmin(utilite[conflicting_elements])
            chromosome[conflicting_elements[min_value_index]] = 0
        end
    end
    return chromosome
end

# ----------------------------------------------------------------------------------------------------

# Valide et répare
function survivantEnfant(enfMut)
    if valide(enfMut)
        return enfMut
    else
        for i in findall(x -> x == 1, enfMut)
            enfMut[i] = 0
            if valide(enfMut)
                return enfMut
            end
            enfMut[i] = 1
        end
    end
end

# ----------------------------------------------------------------------------------------------------
# Recupere la nouvelle generation comme population de base
function changeGeneration(newGen::Vector{Any}, popSize::Int64)
    if isempty(newGen)
        error("New generation is empty")
    end
    if length(newGen) < popSize
        error("New generation has fewer elements than the population size")
    end
    return newGen[1:popSize]
end
# ----------------------------------------------------------------------------------------------------

function IdentifieMeilleur(population)
    best_value = 0
    best_solution = nothing
    for chromosome in population
        value = fitness(chromosome)
        if value > best_value && valide(chromosome)
            best_value = value
            best_solution = chromosome
        end
    end
    return best_value, best_solution
end

# ----------------------------------------------------------------------------------------------------

function valide(sol)
    colonnes = findall(x -> x == 1, sol)
    cidx = [colonne[1] for colonne in colonnes]

    if any(sum(A[:, cidx], dims=2) .> 1)
        return false
    end
    return true
end 

# ----------------------------------------------------------------------------------------------------

population_size = 10
generations = 100
mutation_rate = 0.01
crossover_rate = 0.7




function genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)

    population = []

    
    for i in 1:population_size
        alpha = rand(0.1 : 0.7)
        s = grconst(A, C, alpha)
        push!(population, s)
    end

    fitness_values = evaluate_population(population)    

    for _ in 1:generations
        #newGen=[]
        
        p1 = selection(population)
        p2 = selection(population)
        
        enfantA, enfantB = crossover(p1, p2, crossover_rate)
        
        mutA = survivantEnfant(mutation(enfantA, mutation_rate))
        #push!(newGen, mutA)
        mutB = survivantEnfant(mutation(enfantB, mutation_rate))
        #push!(newGen, mutB)
        
        sorted_indices = sortperm(fitness_values, rev=true)
        population = population[sorted_indices[1:population_size]]
        fitness_values = fitness_values[sorted_indices[1:population_size]]
        
        
        push!(population, mutA)
        push!(population, mutB)
        push!(fitness_values, fitness(mutA))
        push!(fitness_values, fitness(mutB))
        # population = changeGeneration(newGen, population_size)
    end
    
    best_value, best_solution = IdentifieMeilleur(population)
    
    return best_solution, best_value
end

best_solution, best_value = genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)