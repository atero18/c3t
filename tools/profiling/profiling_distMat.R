source("tools/profiling/setup_profiling.R") # nolint

prefixeSauvegarde <- "DistMat"


# Comparaison d'accès en fonction du type de stockage
# et par rapport aux matrices classiques
temps_acces <- function(M)
{
  n_int <- nrow(M)

  data <- data.frame(acces_unitaire_total = NA)

  set.seed(123L)

  # Accès unitaires
  nbAccesUnitaires <- ceiling(n_int / 5L)
  vecteurI <- sample(seq_len(n_int), nbAccesUnitaires, replace = TRUE)
  vecteurJ <- sample(seq_len(n_int), nbAccesUnitaires, replace = TRUE)
  tempsCumule <- 0.0

  for (k in seq_len(nbAccesUnitaires))
  {
    i <- vecteurI[k]
    j <- vecteurJ[k]

    t <- Sys.time()
    M[i, j] == 1.0
    tempsCumule <- tempsCumule + Sys.time() - t
  }
  data$acces_unitaire_total <- tempsCumule
  #data$acces_unitaire_moyen = t / nbAccesUnitaires # nolint

  # Accès lignes / colonnes
  nbAccesLC <- ceiling(n_int / 10L)
  vecteurIndices <- sample(seq_len(n_int), nbAccesLC, replace = TRUE)
  tempsCumule <- 0.0
  for (k in seq_len(nbAccesLC))
  {
    indice <- vecteurIndices[k]

    t <- Sys.time()
    M[indice, ] == 1.0
    tempsCumule <- tempsCumule + Sys.time() - t
  }
  data$acces_ligne_total <- tempsCumule
  #data$acces_ligne_moyen = t / nbAccesLC # nolint

  tempsCumule <- 0.0
  for (k in seq_len(nbAccesLC))
  {
    indice <- vecteurIndices[k]

    t <- Sys.time()
    M[indice, ] == 1.0
    tempsCumule <- tempsCumule + Sys.time() - t
  }
  data$acces_colonne_total <- tempsCumule
  #data$acces_colonne_moyen = t / nbAccesLC # nolint


  # Accès sous-matrice
  nbSousMatrices <- ceiling(n_int / 5L)
  tailleI <- rbinom(n_int, nbSousMatrices, 1L / 2L)
  tailleJ <- rbinom(n_int, nbSousMatrices, 1L / 2L)
  tempsCumule <- 0.0
  for (k in seq_len(nbSousMatrices))
  {
    vecteurI <- sample(seq_len(n_int), size = tailleI[k])
    vecteurJ <- sample(seq_len(n_int), size = tailleJ[k])

    t <- Sys.time()
    M[vecteurI, vecteurJ] == 1L
    tempsCumule <- tempsCumule + Sys.time() - t

  }
  t <- Sys.time() - t
  data$temps_acces_sousmat_total <- tempsCumule

  # Accès matrice complète
  t <- Sys.time()
  M[] == 1.0
  t <- Sys.time() - t
  data$temps_acces_complete <- t

  data[1L, ] <- as.numeric(data[1L, ])
  data
}

n <- 3000L

set.seed(123L)
contexte <- gen_context(n, nbQuantitatives_int = 1L)

setClass("matriceS4", slots = list(mat = "ANY"))
setMethod("nrow", signature = "matriceS4", function(x) nrow(x@mat))

setMethod("[", signature(x = "matriceS4", i = "numeric",
                         j = "missing", drop = "ANY"),
          function(x, i, j, drop) x[i, seq_len(nrow(x))])

setMethod("[", signature(x = "matriceS4", i = "missing",
                         j = "numeric", drop = "ANY"),
          function(x, i, j, drop) (x[seq_len(nrow(x)), j]))

setMethod("[", signature(x = "matriceS4", i = "missing",
                         j = "missing", drop = "ANY"),
          function(x, i, j, drop) (x[seq_len(nrow(x)), seq_len(nrow(x))]))


setMethod("[", signature(x = "matriceS4", i = "numeric",
                         j = "numeric", drop = "ANY"),
          function(x, i, j, drop) (x@mat[i, j]))


base <- as.matrix(dist(contexte))
matriceS4 <- new("matriceS4", mat = base)
dMvecteur <- constructor_DistMat(base, modeStockage = "vector")
dMmatrice <- constructor_DistMat(base, modeStockage = "matrix")

profvis({tempsBase <- temps_acces(base)}, simplify = FALSE)
tempsMatS4 <- temps_acces(matriceS4)
tempsdMvecteur <- temps_acces(dMvecteur)
profil1 <- profvis({tempsdMmatrice <- temps_acces(dMmatrice)}, simplify = FALSE)
profil1



comparaisonTemps <- rbind(tempsBase,
                          tempsMatS4,
                          tempsdMvecteur,
                          tempsdMmatrice)

for (k in colnames(comparaisonTemps))
{
  comparaisonTemps[, k] <-
    as.numeric(comparaisonTemps[, k]) / as.numeric(comparaisonTemps[1L, k])
}


comparaisonTemps$stockage <- c("matriceR",  "matS4",
                               "distMat_vec", "distMat_mat")
