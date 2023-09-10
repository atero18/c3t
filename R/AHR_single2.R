#' @include partition_class.R
#' @include AHCTree.R
#' @include pbCon.R


AHR_single2 <- function(pb, linkage = "saut_max", partitionInit = NULL, # nolint: cyclocomp_linter
                        minNbClusters = 1L,
                        arretMin = FALSE, modeEvaluationDistInter = 2L,
                        calculDistComplet = TRUE,
                        fusionConstraint = FALSE,
                        fusionConstraintMode = "deterministic")
{
  linkage <- checkLinkage(linkage)
  assertFusionConstraint(fusionConstraint)
  fusionConstraint <- corresponding_fusion_const(fusionConstraint)

  n <- pb$n()

  if (is.null(partitionInit))
    partitionInit <- seq_len(n)

  # Gestion de la contrainte et du mode de fusion
  penaliteTaille <- fusionConstraint == "penaliteTaille"

  if (fusionConstraint %in% c(FALSE, "no", "non", "penaliteTaille"))
    fusionConstraint <- FALSE

  else
  {
    typeContrainteFusion <- codeGroupesLiaison(fusionConstraint)

    if (typeContrainteFusion == codeGroupesLiaison("default"))
      fusionConstraint <- FALSE

    else
      fusionConstraint <- TRUE
  }

  if (fusionConstraint)
    fusionConstraintDeter <- fusionConstraintMode == "deterministic"

  # Transformation de la partition initiale pour obtenir une séquence
  # d'ids pour les clusters
  partition <- standardize_partition(partitionInit)
  nbClusters <- max(partition)
  nomClusters <- as.character(seq_len(nbClusters))
  listeClusters <- lapply(seq_len(nbClusters),
                          function(c) which(partition == c))
  names(listeClusters) <- nomClusters
  nbElementsClusters <- lengths(listeClusters)
  nbSingletons <- sum(nbElementsClusters == 1L)
  remainingClusters <- rep(TRUE, nbClusters)

  clustersSizes <- .clusters_sizes(partition, pb$sizes)

  # Vérification de la contrainte minimale
  nbClustersMin <- pb$nbTooSmallClusters(partition)
  contrainteMinVerifiee <- nbClustersMin == 0L

  # Définition d'une taille "correcte" pour une région
  tailleCorrecte <- 0.0
  if (pb$hasMinConstraint())
    tailleCorrecte <- pb$m

  if (pb$hasMaxConstraint())
    tailleCorrecte <- mean(c(tailleCorrecte, pb$M))

  # Recherche des fusions possibles
  contiguiteClasses <-
    pb$listeContiguiteClasses(partition, avecTailleFusion = TRUE)

  # `Inf` : fusions that could exist in the future
  linkageDistances <- matrix(Inf, nrow = nbClusters, ncol = nbClusters)

  if (pb$hasMaxConstraint())
  {
    indexes <- as.matrix(expand.grid(seq_len(nbClusters),
                                     seq_len(nbClusters)))

    indexes <- indexes[indexes[, 1L] < indexes[, 2L], , drop = FALSE]
    fusionSizes <- clustersSizes[indexes[, 1L]] + clustersSizes[indexes[, 2L]]
    maskTooHeavy <- fusionSizes > pb$M
    # `NaN` : fusions that can't be done now or later
    if (any(maskTooHeavy))
    {
      linkageDistances[indexes[maskTooHeavy, 1L:2L, drop = FALSE]] <- NaN
      linkageDistances[indexes[maskTooHeavy, 2L:1L, drop = FALSE]] <- NaN
      contiguiteClasses <-
        contiguiteClasses[contiguiteClasses[, "fusionSize"] <= pb$M, ,
                          drop = FALSE]
    }
  }

  # Clusters that are not on the same connected component will never
  # be able to merge.
  if (!pb$isConnected())
  {
    connectedComponents <- split(seq_len(n), pb$connected_components())
    for (i in seq_along(connectedComponents))
    {
      for (j in seq_along(connectedComponents)[-i])
      {
        linkageDistances[connectedComponents[[i]],
                         connectedComponents[[j]]] <- NaN
      }
    }
  }

  diag(linkageDistances) <- NaN

  # Fusions between contiguous clusters and with respect of the maximal
  # size constraint are allowed (`NA` = authorized fusions with
  # distances to be determined)
  indexes <- contiguiteClasses[, 1L:2L, drop = FALSE]


  if (is_medoid_linkage(linkage))
  {
    medoidLinkage <- TRUE
    medoids <- medoids_partition(pb, partition)
  }
  else
  {
    medoidLinkage <- FALSE
    medoids <- NULL
  }



  # Calcul of the different linkage distances
  distances <- .calcul_distances_inter(pb, partition,
                                       indexes, linkage,
                                       medoides = medoids)

  linkageDistances[indexes] <- distances
  linkageDistances[indexes[, 2L:1L, drop = FALSE]] <- distances
  rm(distances)


  # Définition d'un modèle de partition pour l'enregistrement dans l'arbre
  contraintes <-         c(TRUE, FALSE,  TRUE)
  names(contraintes) <- c("connectivity", "min", "max")

  if (contrainteMinVerifiee)
    scoreMinSizeConst <- 0.0
  else
  {
    scoreMinSizeConst <-
      score_constraints_min(pb$m, clustersSizes = clustersSizes)
  }

  partitionRef <-
    constructor_Partition(pb$m, partition = partition, method = "RAH",
                          contiguity = TRUE,
                          minConstraint = contrainteMinVerifiee,
                          maxConstraint = TRUE,
                          scoreMinSizeConst = scoreMinSizeConst,
                          scoreMaxSizeConst = 0.0)


  # Classe ArbreCAH stockant les partitions
  arbreCAH <- constructor_AHCTree(pb, partitionRef)
  rm(partitionRef)

  nbClusters <- nbClusters - 1L

  nearestNeighbor <- nearest_neighbor(linkageDistances, inner = TRUE)

  satisfaisable <- !all(is.nan(nearestNeighbor))


  while (!(arretMin && contrainteMinVerifiee) &&
        satisfaisable &&
        nbClusters >= minNbClusters)
  {
    # Calcul des distances inter-classes non encore calculées /
    # qui ont été modifiées par une fusion


    clustersMin <-
      which(linkageDistances ==
              min(linkageDistances[remainingClusters, remainingClusters],
                  na.rm = TRUE), arr.ind = TRUE)

    if (nrow(clustersMin) == 0L)
    {
      satisfaisable <- FALSE
      break
    }
    #minDist <- min_distance(linkageDistances, nearestNeighbor)

    cluster1 <- clustersMin[1L, 1L]
    cluster2 <- clustersMin[1L, 2L]
    #cluster1 <- minDist[1L]
    #cluster2 <- minDist[2L]

    remainingClusters[cluster2] <- FALSE

    if (nbSingletons > 0L)
    {
      cluster1EstSingleton <- nbElementsClusters[cluster1] == 1L
      cluster2EstSingleton <- nbElementsClusters[cluster2] == 1L

      nbSingletons <- nbSingletons - cluster1EstSingleton - cluster2EstSingleton
      if (fusionConstraint && nbSingletons == 0L)
        typeContrainteFusion <- max(codeGroupesLiaison("minMin"),
                                    typeContrainteFusion)
    }
    if (nbClustersMin > 0L)
    {
      cluster1EstMin <- clustersSizes[cluster1] < pb$m
      cluster2EstMin <- clustersSizes[cluster2] < pb$m

      # Si au moins un cluster est min, comme un disparaît
      # on perd un cluster min
      if (cluster1EstMin || cluster2EstMin)
        nbClustersMin <- nbClustersMin - 1L

      # Si les deux sont min et que la fusion engendre un cluster assez gros
      # on perd un autre min
      if (cluster1EstMin &&
          cluster2EstMin &&
          clustersSizes[cluster1] + clustersSizes[cluster2] >= pb$m)
        nbClustersMin <- nbClustersMin - 1L

      if (nbClustersMin == 0L)
      {
        fusionConstraint <- FALSE
        contrainteMinVerifiee <- TRUE
      }
    }

    # Partition update
    oldPartition <- partition
    partition[partition == cluster2] <- cluster1
    listeClusters[[cluster1]] <- c(listeClusters[[cluster1]],
                                   listeClusters[[cluster2]])
    listeClusters[[cluster2]] <- numeric(0L)


    # Medoid update
    if (medoidLinkage)
    {
      elementsCluster1 <- listeClusters[[cluster1]]
      medoids[cluster1] <- elementsCluster1[medoid(pb[elementsCluster1,
                                                      elementsCluster1,
                                                      drop = FALSE])]
      medoids[cluster2] <- NaN
    }

    partitionArbre <- arbreCAH[[nbClusters]]
    partitionArbre$partition <- standardize_partition(partition)
    partitionArbre$contraintes["min"] <- contrainteMinVerifiee

    if (contrainteMinVerifiee)
      partitionArbre$scoreSizeConsts[c("min", "total")] <- 0.0
    else
    {
      partitionArbre$scoreSizeConsts[c("min", "total")] <-
        score_constraints_min(pb, clustersSizes = clustersSizes)
    }


    # Ajout des nouvelles contiguités à cluster1
    newContiguitiesC1 <- update_contiguities(linkageDistances,
                                             cluster1,
                                             cluster2)

    if (length(newContiguitiesC1) > 0L)
    {
      linkageDistances[cluster1, newContiguitiesC1] <-
        linkageDistances[newContiguitiesC1, cluster1] <- NA_real_
    }

    # Ajout des nouvelles contiguïtés à cluster2 (peut être utilisé pour
    # la mise à jour des distances)
    newContiguitiesC2 <- update_contiguities(linkageDistances,
                                             cluster2,
                                             cluster1)

    if (length(newContiguitiesC2) > 0L)
    {
      linkageDistances[cluster2, newContiguitiesC2] <-
        linkageDistances[newContiguitiesC2, cluster2] <- NA_real_
    }

    # Update sizes
    clustersSizes <- update_sizes(clustersSizes, cluster1, cluster2)

    # Update authorized fusions, with respect to maximum size constraint
    if (pb$hasMaxConstraint())
    {
      linkageDistances <- .remove_too_heavy_fusions(linkageDistances,
                                                    clustersSizes, cluster1,
                                                    cluster2, pb$M)
    }



    # Update linkage distances
    linkageDistances <-
      .update_linkage(linkageDistances,
                      cluster1, cluster2,
                      linkage,
                      usefulClusters =
                        seq_along(nomClusters)[remainingClusters], # nolint: indentation_linter
                      elemDistances = pb,
                      partitionBefore = oldPartition,
                      nbElementsClustersBefore = nbElementsClusters,
                      availableClusters = remainingClusters,
                      medoides = medoids)


    nbElementsClusters[cluster1] <- nbElementsClusters[cluster1] +
                                    nbElementsClusters[cluster2]
    nbElementsClusters[cluster2] <- NaN
    nbClusters <- nbClusters - 1L

    # S'il n'y a plus de fusions possibles
    if (all(is.nan(linkageDistances) | is.infinite(linkageDistances)))
    {
      satisfaisable <- FALSE
    }

  }

  return(arbreCAH[!(names(arbreCAH) %in% 1L:nbClusters)])
}
