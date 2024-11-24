using LinearAlgebra

include("loadSPP.jl")
include("EI2.jl")
fname = "../Data/pb_200rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)


function evaluate_population(population)
    return [fitness(chromosome) for chromosome in population]
end

# ----------------------------------------------------------------------------------------------------

function fitness(chromosome)
    #équivalant à z(x, values)
    return dot(chromosome, C) 
end

# ----------------------------------------------------------------------------------------------------

function selection(population)

    #sélection par tournois
    cand2 = population[rand(1:length(population))]
    cand1 = population[rand(1:length(population))]

    if fitness(cand2) > fitness(cand1)
        cand1 = cand2
    end

    return cand1  
end 
# ----------------------------------------------------------------------------------------------------

function crossover(p1, p2, crossover_rate)
    enfantA=copy(p1)
    enfantB=copy(p2)
    
    if rand() < crossover_rate #crossover_rate est élevé : il y a souvent des croisements mais pas systématiquement
        crossover_point = rand(1:length(p1)) # on choisit un point de coupe aléatoire
        enfantA = vcat(p1[1:crossover_point], p2[crossover_point+1:end])
        enfantB = vcat(p2[1:crossover_point], p1[crossover_point+1:end])
        
    end
    
    return enfantA, enfantB
end

# ----------------------------------------------------------------------------------------------------

function mutation(chromosome, mutation_rate)
    for i in eachindex(chromosome)
        if rand() < mutation_rate #mutation_rate pas : il y a pas souvent de mutations 
            chromosome[i] = 1 - chromosome[i]
        end
    end
    return repare(chromosome) #la solution est systématiquement réparée
end

# ----------------------------------------------------------------------------------------------------

function repare(chromosome)

    if valide(chromosome)
        return chromosome
    end
    utilite = map(i -> C[i] / sum(A[:, i]), 1:n) #utilité de chaque colonne
    while !valide(chromosome)
        conflict_indices = findall(x -> x > 1, A * chromosome) #indices des colonnes conflictuelles
        
        for idx in conflict_indices
            conflicting_elements = findall(x -> x == 1, A[idx, :]) #indices des éléments conflictuels
            min_value_index = argmin(utilite[conflicting_elements]) #on prend l'élément de plus petite utilité
            chromosome[conflicting_elements[min_value_index]] = 0 #et on le retire
        end
    end
    return chromosome
end

# ----------------------------------------------------------------------------------------------------

function IdentifieMeilleur(population, fitness_values, sorted_indices)
    best_index = sorted_indices[1]
    best_solution = population[best_index]
    best_value = fitness_values[best_index]
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

function genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)
    
    # Initialisation de la population, avec la construction de grasp car les solutions sont meilleures + les solutions sont valides 
    population = [grconst(A, C, rand(0.1:0.7)) for _ in 1:population_size]
    fitness_values = evaluate_population(population)    

    for _ in 1:generations
        
        p1 = selection(population)
        p2 = selection(population)
        
        enfantA, enfantB = crossover(p1, p2, crossover_rate)
        
        mutA = mutation(enfantA, mutation_rate)
        mutB = mutation(enfantB, mutation_rate)

        push!(population, mutA)
        push!(population, mutB)
        push!(fitness_values, fitness(mutA))
        push!(fitness_values, fitness(mutB))

        sorted_indices = sortperm(fitness_values, rev=true) #on trie les indices des solutions par ordre décroissant de fitness
        population = population[sorted_indices[1:population_size]] #on garde les meilleures solutions
        fitness_values = fitness_values[sorted_indices[1:population_size]] #on garde les meilleures valeurs de fitness

    end
    
    best_solution = population[1]
    best_value = fitness_values[1]


    return best_solution, best_value
    #return best_value

end

# ----------------------------------------------------------------------------------------------------

population_size = 10
generations = 10
mutation_rate = 0.01
crossover_rate = 0.9


best_solution, best_val = genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)
# function run_multiple_times(A, C, population_size, generations, mutation_rate, crossover_rate, num_runs)
#     results = Float64[]
#     times = Float64[]

#     for _ in 1:num_runs
#         elapsed_time = @elapsed best_value = genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)
#         push!(results, best_value)
#         push!(times, elapsed_time)
#     end

#     min_value = minimum(results)
#     max_value = maximum(results)
#     mean_value = mean(results)
#     avg_time = mean(times)

#     return min_value, max_value, mean_value, avg_time
# end

#run_multiple_times(A, C, population_size, generations, mutation_rate, crossover_rate, 10)