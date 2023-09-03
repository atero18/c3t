rm(list = ls())
source("profiling/core_profiling.R")

prefixeSauvegarde = "RAH"

# Comparaison stockage vecteur et matrice

data(grid_queen_vide0_metropole0_x7_y7_indivMoy100_quant3_qual0)
g7x7 = grid_queen_vide0_metropole0_x7_y7_indivMoy100_quant3_qual0
matrice_dist7x7 = as.matrix(dist(g7x7$context, method = "euclidean"))

pbVecteur <- c3t_grid_simulation(7L, 7L, distance = "euclidean",
                                   m = 0.0, M = Inf,
                                    calculToutesValeurs = TRUE,
                                   storageMode = "vector")

pbMatrice <- c3t_grid_simulation(7L, 7L, distance = "euclidean",
                                   m = 0.0, M = Inf,
                                   calculToutesValeurs = TRUE,
                                   storageMode = "matrix")


set.seed(123L)
m <- mark(vecteur = AHR_single(pbVecteur, linkage = "saut_min"),
          matrice = AHR_single(pbMatrice, linkage = "saut_min"),
          iterations = 4, check = FALSE)


# Comparaison mode calcul distances


m <- mark(mode2 = AHR_single(pb20x20, linkage = "saut_min", modeEvaluationDistInter = 2),
          mode3 = AHR_single(pb20x20, linkage = "saut_min", modeEvaluationDistInter = 3),
          mode4 = AHR_single(pb20x20, linkage = "saut_min", modeEvaluationDistInter = 4),
          mode6 = AHR_single(pb20x20, linkage = "saut_min", modeEvaluationDistInter = 6))

m <- mark(mode2 = AHR_single(pb7x7, linkage = "saut_min", modeEvaluationDistInter = 2),
          mode3 = AHR_single(pb7x7, linkage = "saut_min", modeEvaluationDistInter = 3),
          mode4 = AHR_single(pb7x7, linkage = "saut_min", modeEvaluationDistInter = 4),
          mode6 = AHR_single(pb7x7, linkage = "saut_min", modeEvaluationDistInter = 6), min_iterations = 3)
m

# 11/07/2023 : À ce jour le meilleur mode semble être :
# - le troisième pour la complexité temporelle (très proche du second cependant)
# - le 4 pour l'allocation mémoire, le 2 en 2nde position
# - le pire (et de loin) étant le mode 4 pour la complexité temporelle
# le mode retenu sera le 2.

# Unique itération

# Étude sans calcul de distance en interne et sans contrainte de taille
# -- Grille 7x7
pb7x7 <- c3t_grid_simulation(7L, 7L, distance = "euclidean",
                               m = 0.0, M = Inf,
                               calculToutesValeurs = TRUE)
set.seed(123L)
profil1 <- profvis::profvis({res = AHR_single(pb7x7, linkage = "saut_min")})
#sauvegarde_profiling(profil1, prefixeSauvegarde, "profil1")
profil1

# -- Grille 20x20
pb20x20 <- c3t_grid_simulation(20L, 20L, distance = "euclidean",
                               m = 0.0, M = Inf,
                               calculToutesValeurs = TRUE)
set.seed(123L)
profil2 <- profvis::profvis({res = AHR_single(pb20x20, linkage = "saut_min")})
profil2

# 22/06/2023 : La recherche initiale des classes contigues prenait beaucoup trop
# de temps (600 ms sur 1000). Éviter un appel de `which` à chaque fois en stockant
# la liste des éléments de chaque cluster a permis une grosse réduction
# Le calcul des distances était relancé à chaque itération, correction.
# Si la distance d'ensemble est saut min ou saut max, optimisation du calcul
# en faisant une mise à jour à chaque fin d'itération à partir des données déjà
# calculées
# Le calcul de masques logiques prend du temps. Optimisation du calcul de masques
# en prenant mieux en compte la symétrie des distances (réduction de 45%)
# Le calcul des distances inter-classes non déterminées à chaque début d'itération
# fait moins appel a un calcul de masque logique
# Utilisation moindre des masques en modifiant .get_elem_symSparseMat pour que la fonction
# gère une demande d'accès à un ensemble de données plutôt qu'un accès un à un

# -- Grille 30x50
pb30x50 <- c3t_grid_simulation(30L, 50L, distance = "euclidean",
                                 m = 0.0, M = Inf,
                                 calculToutesValeurs = TRUE)
set.seed(123L)
profil3 <- profvis::profvis({res = AHR_single(pb30x50, linkage = "saut_min", Cpp = TRUE)})
profil3

# 23/06/2023 : L'accès aux matrices de distances creuses prenant du temps, stockage
# des données de distances dans un vecteur, en considérant la symétrie et la nullité au
# niveau de la diagonale. Beaucoup plus rapide (70% environ) mais pour de grandes dimensions
# la place prise en mémoire est importante ((n-1)n/2 éléments à stocker).

# 28/06/2023 : Vérifier qu'une partition n'a pas besoin d'être simplifiée accélère le calcul
# L'utilisation de la fonction `match` plutôt que as.factor pour la simplification des partitions
# est beaucoup plus efficace, rendant le calcul quasi instantaté (gain de 2 secondes)

# 29/06 : L'utilisation d'un masque sur la recherche des indices i > j pour symMat est plus efficace
# que l'utilisation de pmax et pmin.

# 30/06/2023 : Tentative d'utiliser une récupération en une étape des distances inter éléments pour
# le calcul des distances inter-clusters. Le résultat est semblable voire un peu moins bon en allocation mémoire.
# Ce qui prend le plus de temps dans le calcul des distances inter-clusters c'est l'accès aux données

# -- Grille 50x70
pb50x70 <- c3t_grid_simulation(50L, 70L, distance = "euclidean",
                                 m = 0.0, M = Inf,
                                 calculToutesValeurs = TRUE)
set.seed(123L)
profil4 <- profvis::profvis({res = AHR_single(pb50x70, linkage = "saut_min")})
profil4

# 28/06/2023 : Amélioration du calcul de la contiguité inter-classe, important lorsque
# l'on ne démarre pas le calcul du bas de l'arbre.

profvis::profvis({res = AHR_pb(pb7x7, linkage = "saut_min", nbTries = 1)})



# Comparaison calcul complet des distances au début ou non
#
pb20x20 <- c3t_grid_simulation(20L, 20L, distance = "euclidean",
                                 m = 0.0, M = Inf,
                                 calculToutesValeurs = FALSE)

pb20x20_complete <- c3t_grid_simulation(20L, 20L, distance = "euclidean",
                                          m = 0.0, M = Inf,
                                          calculToutesValeurs = TRUE)

m <- mark(nonComplet = as_tibble(AHR_single(pb20x20$copy(), linkage = "saut_min", calculDistComplet = FALSE)),
          complet = as_tibble(AHR_single(pb20x20_complete, linkage = "saut_min", calculDistComplet = TRUE)), min_iterations = 4)

# Ajout de contraintes sur les fusions (relatives aux tailles)
profvis::profvis({res = AHR_single(pb4x4, linkage = "saut_min", fusionConstraint = "minAny")})

# Comparaison globale
listePb = list("49" = pb7x7,
               "400" = pb20x20,
               "1500" = pb30x50,
               "3500" = pb50x70)

n = c(49, 400, 1500, 3500)
grilleProfiling = data.frame(n = rep(n, each = 3), Cpp = rep(c(T, T, F), length(n)), calculDistComplet = rep(c(T, F, F), length(n)))
comparaisonRAH = press(.grid = grilleProfiling,
            {
              mark(AHR_single(listePb[[as.character(n)]], linkage = "saut_min", Cpp = Cpp, calculDistComplet = calculDistComplet), filter_gc = FALSE, check = FALSE)
            }
)

comparaisonRAH = select(comparaisonRAH, -memory)

#save(comparaisonRAH, file = "profiling/saves/.....Rdata")

# t0Comp = Sys.time()
# t0Comp
# RAHgrandeGrilleComp = AHR_single(pb175x200, linkage = "saut_min", Cpp = TRUE, calculDistComplet = TRUE)
# t1Comp = Sys.time()
# t1Comp
# tempsTotalComp = t1Comp - t0Comp
# tempsTotalComp
#
# t0NonComp = Sys.time()
# t0NonComp
# RAHgrandeGrilleNonComp = AHR_single(pb175x200, linkage = "saut_min", Cpp = TRUE, calculDistComplet = FALSE)
# t1NonComp = Sys.time()
# t1NonComp
# tempsTotalNonComp = t1NonComp - t0NonComp
# tempsTotalNonComp
