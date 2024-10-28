include("loadSPP.jl")
fname = "../Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)

C = Int.(C)
A = Int.(A)


function utility(A, C, n)
#retourne l'indice de l'élement de l'indice dans C de l'élément ayant la plus grande utilité  
    
    # Initialiser une liste des utilités
    U = zeros(n)

    # Calculer l'utilité pour chaque variable (colonne de A)
    for i in 1:n
        # L'utilité est C[i] divisé par la somme des contraintes de la colonne i
        U[i] = C[i] /sum(A[: , i]) 
    end
    
    # println("U : ", U)
    return U
end


function greedy2(C, A)
    m, n = size(A)
    #println("m = ", m, " n = ", n)

    # Trier les indices en fonction des ratios
    u = sortperm(utility(A, C, n), rev=true)

    # Initialiser le vecteur x avec des zéros
    x = zeros(Int, n)

    # Ajouter des 1 aux bons indices
    for un in u
        #println("un = ", nounde)
        feasible = true
        for i in 1:m
            if A[i, un] == 1
                for j in 1  :n
                    if x[j] == 1 && A[i, j] == 1
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
            x[un] = 1
        end
    end
    #affiche x avec ses indices
    
    #for (index, value) in enumerate(x)
    #    println("Index: $index, Valeur: $value")
    #end
    z = sum(C[i] * x[i] for i in 1:n)
    #println(res)
    return x, z
end


# function conflict(A, C, ie, ivar)
#     # Mettre à jour l'ensemble des candidats en supprimant les éléments en conflit : 
#     # Dans A : la colonne correspondant à l'element d'indice ie, ainsi que toutes les colonnes ayant des 1 sur la meme ligne que ie en a
#     # Dans C : l'element ie ainsi que tous les elements dont les colonnes ont ete supprime
    
#     # Identifier la colonne à supprimer
#     col_to_remove = A[:, ie]
#     # Trouver les lignes où la colonne a des 1
#     rows_with_ones = findall(col_to_remove .== 1)

#     # Identifier les colonnes à supprimer
#     columns_to_remove = []
#     for row in rows_with_ones
#         columns_with_ones = findall(A[row, :] .== 1)
#         for col in columns_with_ones
#             push!(columns_to_remove, col[1])  # Ajouter les colonnes à supprimer
#         end
#     end
    
#     # Supprimer la colonne correspondante à ie et toutes les colonnes à supprimer
#     A = A[:, setdiff(1:size(A, 2), [ie; collect(columns_to_remove)])]
        

#     remaining_indices = setdiff(1:length(C), collect(columns_to_remove))
#     C = C[remaining_indices]
#     ivar = ivar[remaining_indices]  # Corriger ivar après suppression
    
#     return A, C, ivar
# end


