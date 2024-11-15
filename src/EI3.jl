include("loadSPP.jl")
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
m, n = size(A)


function initialize_population(population_size, chromosome_length)
    return [rand(Bool, chromosome_length) for _ in 1:population_size]
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

function selection(population, fitness_values)
    cand2 = population[rand(1:length(population))]
    cand1 = population[rand(1:length(population))]

    if fitness_values[cand2] > fitness_values[cand1]
        cand1 = cand2
    end

    filter!(elem -> elem[1] != cand1, population)

    return cand1  

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

function mutation(population, mutation_rate)
    for chromosome in population
        for i in eachindex(chromosome)
            if rand() < mutation_rate
                chromosome[i] = !chromosome[i]
            end
        end
    end
    return population
end

# ----------------------------------------------------------------------------------------------------

# Selectionne un individu survivant entre deux individus
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
function changeGeneration(newGen, popSize)

    pop = Vector(undef, popSize)
    NbRealisable = 0
    maxFitness = 0
    minFitness = 100

    for i=1:popSize
        pop[i] = pop!(newGen)
        ind, fitness, realisable = pop[i]
        if realisable
            NbRealisable +=1
            minFitness = min(fitness,minFitness)
        end
        maxFitness = max(fitness,maxFitness)
    end
    println("Nbre Realisable = ", NbRealisable, " minFitnessRealisable = ", minFitness, " maxFitness = ", maxFitness)
    return pop
end

# ----------------------------------------------------------------------------------------------------

function IdentifieMeilleur(population)
    best_value = -Inf
    best_solution = nothing
    for chromosome in population
        value = fitness(chromosome)
        if value > best_value
            if valide(chromosome)
                best_value = value
                best_solution = chromosome
            end
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

population_size = 100
generations = 1000
mutation_rate = 0.01
crossover_rate = 0.7

best_solution, best_value = genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)



function genetic_algorithm(A, C, population_size, generations, mutation_rate, crossover_rate)
    population = initialize_population(population_size, size(A, 2))

    # best_solution = nothing
    # best_value = -Inf

    for gen in 1:generations
        newGen=[]

        fitness_values = evaluate_population(population)

        p1 = selection(population, fitness_values)
        p2 = selection(population, fitness_values)

        enfantA, enfantB = crossover(p1, p2, crossover_rate)

        mutA = survivantEnfant(mutation(enfantA, mutation_rate))
        push!(newGen, mutA)
        mutB = survivantEnfant(mutation(enfantB, mutation_rate))
        push!(newGen, mutA)

        population = changeGeneration(newGen, population_size)
    end

    best_value, best_solution = IdentifieMeilleur(population)

    return best_solution, best_value
end