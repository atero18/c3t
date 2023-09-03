rm(list = ls())
source("profiling/core_profiling.R")

prefixeSauvegarde = "contiguity"

pb7x7 <- c3t_grid_simulation(7L, 7L, distance = "euclidean",
                               m = 0.0, M = Inf,
                               calculToutesValeurs = TRUE)

# Étude du temps de recherche des classes contigues
# -- 1 élément = 1 cluster, contiguité donnée par matrice de contiguïté
# L'élément donnant la contiguité envoyé est la matrice de contiguïté

partition_vec = 1:n
profil1 <- profvis({clusters_contiguity_list(partition_vec, pb7x7$contiguity)})

# -- 1 élément = 1 cluster, contiguité donnée par matrice contiguityMat

partition_vec = 1:n
contiguite = contiguityMat(g$contiguity)
profil2 <- profvis({clusters_contiguity_list(partition_vec, contiguite)})
profil2

# 22/06/2023 : Suppression du namespace checkmate, l'appel du namespace a un cout non nul.
# suppression de la vérification constante des éléments lors de l'appel à get_elem_symSparseMat. La vérification est
# faite en amont lors de l'appel de [i,j]. Réduction d'envriron 90%

# la conversion graphe de contiguïté -> contiguityMat prend du temps. Pour la recherche de contiguité entre classes, si un
# graphe est envoyé alors il est transformé en matrice creuse booléenne.

# 27/06/2023 : Rendre le calcul de position dans la matrice vectoriel a grandement amélioré le temps d'accès à une matrice
# symétrique. De même utiliser rep plutôt que expand.grid est beaucoup plus rapide.

# -- Comparaison avec appel matrice classique et de contiguite

m <- mark(matriceClassique = clusters_contiguity_list(partition_vec, g$contiguity),
          contiguityMat = clusters_contiguity_list(partition_vec, contiguite),
          iterations = 1, check = FALSE)

plot(m)
