source("tools/profiling/setup_profiling.R") # nolint

# -- Grille 7x7

pb7x7 <- c3t_grid_simulation(7L, 7L, distance = "euclidean",
                             m = 0.0, M = Inf,
                             calculToutesValeurs = TRUE)

set.seed(123L)
resPerso <- AHR_single(pb7x7, linkage = "saut_min")
grilleAmeliorant <- enhancer_grid_creation()
grid_enhancer_pb(pb7x7, resPerso[[5L]]$partition,
                 grilleAmeliorant, enregistrerTemps = TRUE)


profil1 <- profvis({res7x7 <- enhancer_pb(pb7x7, resPerso[[5L]]$partition)})
profil1

# 29/06/2023 : L'utilisation du graphe de contiguité et de igraph plutôt
# que de la matrice de contiguité dans la recherche des points transférables
# réduit de 1/3 le temps de calcul

# -- Grille 20x20

pb20x20 <- c3t_grid_simulation(20L, 20L, distance = "euclidean",
                               m = 0.0, M = 3000.0,
                               calculToutesValeurs = TRUE)

set.seed(123L)
resPerso <- AHR_single(pb20x20, linkage = "saut_min")

profil2 <- profvis({res20x20 <- enhancer_pb(pb20x20, resPerso[[5L]]$partition)})
profil2
