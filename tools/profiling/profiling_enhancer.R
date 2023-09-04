source("tools/profiling/setup.R") # nolint

profvis <- function(...) profiling(...,  prefix = "enhancer")

# -- Grille 7x7

pb7x7 <- c3t_grid_simulation(7L, 7L, distance = "euclidean",
                             m = 0.0, M = Inf,
                             calculToutesValeurs = TRUE)

set.seed(123L)
resPerso <- AHR_single(pb7x7, linkage = "single")


profil1 <- profvis({res7x7 <- enhancer_pb(pb7x7, resPerso[[5L]]$partition)},
                   name = "7x7")

# 29/06/2023 : L'utilisation du graphe de contiguité et de igraph plutôt
# que de la matrice de contiguité dans la recherche des points transférables
# réduit de 1/3 le temps de calcul

# -- Grille 20x20

pb20x20 <- c3t_grid_simulation(20L, 20L, distance = "euclidean",
                               m = 0.0, M = 3000.0,
                               calculToutesValeurs = TRUE)

set.seed(123L)
resPerso <- AHR_single(pb20x20, linkage = "single")

profil2 <- profvis({res20x20 <- enhancer_pb(pb20x20,
                                            resPerso[[105L]]$partition)},
                   name = "20x20")
