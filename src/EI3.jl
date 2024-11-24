using LinearAlgebra

include("loadSPP.jl")
include("EI2.jl")
#fname = "../Data/pb_500rnd0100.dat"
# C, A = loadSPP(fname)
# m, n = size(A)


function evaluate_population(population, C)
    return [fitness(chromosome, C) for chromosome in population]
end

# ----------------------------------------------------------------------------------------------------

function fitness(chromosome, C)    #équivalant à z(x, values)
    return dot(chromosome, C) 
end

# ----------------------------------------------------------------------------------------------------

function selection(population, C)

    #sélection par tournois
    cand2 = population[rand(1:length(population))]
    cand1 = population[rand(1:length(population))]

    if fitness(cand2, C) > fitness(cand1, C)
        cand1 = cand2
    end

    return cand1  
end 
# ----------------------------------------------------------------------------------------------------

function crossover(p1, p2, crossover_rate, C, A)
    enfantA=copy(p1)
    enfantB=copy(p2)
    
    if rand() < crossover_rate #crossover_rate est élevé : il y a souvent des croisements mais pas systématiquement
        crossover_point = rand(1:length(p1)) # on choisit un point de coupe aléatoire
        enfantA = vcat(p1[1:crossover_point], p2[crossover_point+1:end])
        enfantB = vcat(p2[1:crossover_point], p1[crossover_point+1:end])
        
    end
    enfantA = reparation(C,A,enfantA)
    enfantB = reparation(C,A,enfantB)
    return enfantA, enfantB
end

# ----------------------------------------------------------------------------------------------------

function mutation(enf, mutation_rate,C, A)
    
    chromosome = copy(enf)
    for i in eachindex(chromosome)
        if rand() < mutation_rate #mutation_rate pas : il y a pas souvent de mutations 
            chromosome[i] = 1 - chromosome[i]
        end
    end
    return reparation(C,A,chromosome)
end

# ----------------------------------------------------------------------------------------------------


function reparation(C,A,sol)
    m, n = size(A)
    utilite = map(i -> C[i] / sum(A[:, i]), 1:n)
    u_idx=sortperm(utilite,rev=false)
    for i=1:m
        if(!valideLigne(A,sol,i))
            for j in eachindex(sol)
                if(A[i,u_idx[j]]==1 && sol[u_idx[j]]==1)
                    sol[u_idx[j]]=0
                    if(valideLigne(A,sol,i))
                        break
                    end
                end
            end
        end
    end
    return sol
end

# ----------------------------------------------------------------------------------------------------

function valideLigne(A,sol,i)
    m, n = size(A)
    sum=0
    for j=1:n
        if(sol[j]==1 && A[i,j]==1)
            sum=1+sum
        end
        if(sum>1)
            return false
        end
    end
    return true
end

# ----------------------------------------------------------------------------------------------------

function genetic_algorithm(fname, population_size = 200, generations = 50, mutation_rate = 0.9, crossover_rate = 0.01)
    
    C, A = loadSPP(fname)
    m, n = size(A)

    # Initialisation de la population, avec la construction de grasp car les solutions sont meilleures + les solutions sont valides 
    population = [grconst(A, C, rand(0.1:0.7)) for _ in 1:population_size]
    fitness_values = evaluate_population(population, C)    

    for _ in 1:generations
        
        p1 = selection(population, C)
        p2 = selection(population, C)
        
        enfantA, enfantB = crossover(p1, p2, crossover_rate, C, A)
        
        mutA = mutation(enfantA, mutation_rate, C, A)
        mutB = mutation(enfantB, mutation_rate, C, A)

        push!(population, mutA)
        push!(population, mutB)
        push!(fitness_values, fitness(mutA, C))
        push!(fitness_values, fitness(mutB, C))

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

# population_size = 200
# generations = 50
# mutation_rate = 0.01
# crossover_rate = 0.9


# best_solution, best_val = genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)
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

# run_multiple_times(A, C, population_size, generations, mutation_rate, crossover_rate, 10)


function experimentationSPP()
    println("didactic", genetic_algorithm("../Data/didactic.dat"))
    println("pb_100rnd0100", genetic_algorithm("../Data/pb_100rnd0100.dat"))
    println("pb_100rnd0300", genetic_algorithm("../Data/pb_100rnd0300.dat"))
    println("pb_200rnd0100", genetic_algorithm("../Data/pb_200rnd0100.dat"))
    println("pb_200rnd0500", genetic_algorithm("../Data/pb_200rnd0500.dat"))
    println("pb_500rnd0100", genetic_algorithm("../Data/pb_500rnd0100.dat"))
    println("pb_500rnd0100", genetic_algorithm("../Data/pb_500rnd0100.dat"))
    println("pb_500rnd1700", genetic_algorithm("../Data/pb_500rnd1700.dat"))
    println("pb_1000rnd0100", genetic_algorithm("../Data/pb_1000rnd0100.dat"))
    println("pb_1000rnd0200", genetic_algorithm("../Data/pb_1000rnd0200.dat"))
    println("pb_2000rnd0100", genetic_algorithm("../Data/pb_2000rnd0100.dat"))
end