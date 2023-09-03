rm(list = ls())
source("profiling/core_profiling.R")

prefixeSauvegarde = "silhouette"

# Grille 7x7

pb7x7 <- c3t_grid_simulation(7L, 7L, distance = "euclidean",
                               m = 0.0, M = Inf,
                               calculToutesValeurs = TRUE)

# Cas avec partition triviale
profvis(silhouette(pb7x7, 1:49))

# Cas avec partition plus complexe
partition = AHR_single(pb7x7)[[10]]$partition
profvis(silhouette(pb7x7, partition))

# Ajout du renvoi des détails
partition = AHR_single(pb7x7)[[10]]$partition
profvis(silhouette(pb7x7, partition, details = TRUE))


# 27/06/2023 : La première méthode de calcul de b prenait trop de temps
# cas une nouvelle sous-matrice de distance était calculée pour chaque élément
# . Maintenant le calcul est fait cluster par cluster

# 28/06/2023 : Pour le calcul de b si le détail n'est pas demandé
# on se restreint au calcul des éléments qui se trouvent dans un cluster
# de taille > 1
# Ajout d'un retour direct lorsque la partition est la partition triviale et
# qu'aucun détail n'est demandé

# Grille 20x20
pb20x20 <- c3t_grid_simulation(20L, 20L, distance = "euclidean",
                               m = 0.0, M = Inf,
                               calculToutesValeurs = TRUE)


# Cas avec partition triviale
profvis(silhouette(pb20x20, 1:400), rerun = TRUE)

# Cas avec partition plus complexe
partition = AHR_single(pb20x20)[[10]]$partition
profvis(silhouette(pb20x20, partition))


# Comparaison avec l'existant

pb30x50 <- c3t_grid_simulation(30L, 50L, distance = "euclidean",
                                 m = 0.0, M = Inf,
                                 calculToutesValeurs = TRUE)

set.seed(123L)
partition = sample(1:20, replace = TRUE, size = 30*50)

m <- mark(c3t = c3t::silhouette(pb30x50, partition),
          cluster = cluster::silhouette(partition, dmatrix = matrice_dist30x50),
          check = FALSE)
m
