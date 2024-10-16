include("loadSPP.jl")
fname = "Data/pb_500rnd0100.dat"
C, A = loadSPP(fname)

C = Int.(C)
A = Int.(A)

function greedy(C, A)::Vector{Int}
    # Initialiser l'ensemble solution S à vide
    S = zeros(Int,length(C))
    ivar = [i for i in 1:length(C)]
    println("Construction gloutonne evaluée d'une solution admissible")
    # Boucle tant que C n'est pas vide
    while !isempty(C)
        println("ivar : ", ivar)
        println("C : ", C)
        # Sélectionner l'indice de l'élément e avec la plus grande utilité dans C
        ie = utility(A, C) 
        println("jselect : ", ie)
        # Incorporer l'élément e dans la solution S
        S[ivar[ie]]=1
        @show S
        # Mettre à jour l'ensemble des candidats en supprimant les éléments en conflit
        A, C, ivar = conflict(A, C, ie, ivar)
        println("_________")
    end

    # Retourner l'ensemble solution S
    return S
end

function utility(A, C)
#retourne l'indice de l'élement de l'indice dans C de l'élément ayant la plus grande utilité  
    n = length(C) # Nombre de variables
    
    # Initialiser une liste des utilités
    U = zeros(n)

    # Calculer l'utilité pour chaque variable (colonne de A)
    for i in 1:n
        # L'utilité est C[i] divisé par la somme des contraintes de la colonne i
        U[i] = C[i] /sum(A[: , i]) 
    end
    
    println("U : ", U)
    return argmax(U)
end

function conflict(A, C, ie, ivar)
    # Mettre à jour l'ensemble des candidats en supprimant les éléments en conflit : 
    # Dans A : la colonne correspondant à l'element d'indice ie, ainsi que toutes les colonnes ayant des 1 sur la meme ligne que ie en a
    # Dans C : l'element ie ainsi que tous les elements dont les colonnes ont ete supprime
    
    # Identifier la colonne à supprimer
    col_to_remove = A[:, ie]
    # Trouver les lignes où la colonne a des 1
    rows_with_ones = findall(col_to_remove .== 1)

    # Identifier les colonnes à supprimer
    columns_to_remove = []
    for row in rows_with_ones
        columns_with_ones = findall(A[row, :] .== 1)
        for col in columns_with_ones
            push!(columns_to_remove, col[1])  # Ajouter les colonnes à supprimer
        end
    end
    
    # Supprimer la colonne correspondante à ie et toutes les colonnes à supprimer
    A = A[:, setdiff(1:size(A, 2), [ie; collect(columns_to_remove)])]
        

    remaining_indices = setdiff(1:length(C), collect(columns_to_remove))
    C = C[remaining_indices]
    ivar = ivar[remaining_indices]  # Corriger ivar après suppression
    
    return A, C, ivar
end

