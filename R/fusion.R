#' Update the partition after fusion.
#'
#' This function updates the partition vector after merging two clusters.
#'
#' @param partition A numeric vector representing the initial partition.
#' @param cluster1 The integer representing the first cluster to merge.
#' @param cluster2 The integer representing the second cluster to merge.
#'
#' @returns A numeric vector representing the updated partition.
#' @noRd
#' @importFrom checkmate assertNumber
#' @examples
#' partition <- c(1, 2, 3, 1, 3)
#' update_partition(partition, 1, 3)
update_partition <- function(partition, cluster1, cluster2)
{
  # Argument verification
  assertPartition(partition)
  assertNumber(cluster1, finite = FALSE)
  assertNumber(cluster2, finite = FALSE)

  if (!(cluster1 %in% partition) || !(cluster2 %in% partition))
  {
    stop("`cluster1` and `cluster2` must be elements of `partition`")
  }

  .update_partition(partition, cluster1, cluster2)
}

.update_partition <- function(partition, cluster1, cluster2)
{
  partition[partition == cluster2] <- cluster1

  return(partition)
}


#' Calculate size after fusions in contiguity list.
#'
#' This function calculates the size of clusters after merging two clusters in
#' a contiguity list.
#'
#' @param liste_contiguite A matrix indicating on each row a contiguity..
#' @param taillesClusters A numeric vector representing the sizes
#' of the clusters.
#'
#' @returns A numeric vector or a dataframe representing the updated
#' sizes after fusions.
#' @noRd
size_after_fusion <- function(liste_contiguite,
                              taillesClusters)
{
  if ("fusionSize" %in% colnames(liste_contiguite))
  {
    taillesInconnues <- is.na(liste_contiguite[, "fusionSize"])
    if (any(taillesInconnues))
    {
      liste_contiguite[taillesInconnues, "fusionSize"] <-
       taillesClusters[as.character(liste_contiguite[taillesInconnues, 1L])] +
       taillesClusters[as.character(liste_contiguite[taillesInconnues, 2L])]

      return(liste_contiguite)
    }
  }
  else
  {
    return(taillesClusters[as.character(liste_contiguite[, 1L])] +
           taillesClusters[as.character(liste_contiguite[, 2L])])
  }
}

update_sizes <- function(oldClusterSizes, cluster1, cluster2)
{
  oldClusterSizes[cluster1] <- oldClusterSizes[cluster1] +
                               oldClusterSizes[cluster2]

  oldClusterSizes[cluster2] <- NaN

  oldClusterSizes
}

# Contiguity must have been updated before
.remove_too_heavy_fusions <- function(oldDistances,
                                      newClusterSizes,
                                      cluster1, cluster2,
                                      M = Inf)
{
  if (is.infinite(M))
    return(oldDistances)

  usefulClusters <- which(!is.nan(newClusterSizes[-c(cluster1, cluster2)]))

  if (length(usefulClusters) == 0L)
    return(oldDistances)

  maskTooHeavy <-
    newClusterSizes[cluster1] + newClusterSizes[usefulClusters] > M

  if (any(maskTooHeavy))
  {
    oldDistances[cluster1, usefulClusters[maskTooHeavy]] <- NaN

    if (!inherits(oldDistances, c("DistMat", "AbstractSymMat")))
      oldDistances[usefulClusters[maskTooHeavy], cluster1] <- NaN
  }

  oldDistances
}

#' Update contiguity list after cluster fusion.
#'
#' This function updates the contiguity list after merging two clusters.
#'
#' @param liste_contiguite A matrix representing the contiguity list.
#' @param cluster1 The cluster identifier for the first cluster to merge.
#' @param cluster2 The cluster identifier for the second cluster to merge.
#' @param clustersSizes A numeric vector representing the sizes
#' of the clusters.
#'
#' @returns A matrix representing the updated contiguity list.
#' @noRd
update_contiguity_list <- function(liste_contiguite,
                                   cluster1,
                                   cluster2,
                                   clustersSizes = NULL)
{
  # Determination of neighborhoods
  masque1A <- liste_contiguite[, 1L] == cluster1
  masque1B <- liste_contiguite[, 2L] == cluster1
  masque1 <- masque1A | masque1B
  voisinsCluster1 <- c(liste_contiguite[masque1A, 2L],
                       liste_contiguite[masque1B, 1L])

  # Reset distances with neighbors of cluster 1
  if ("distance" %in% colnames(liste_contiguite))
    liste_contiguite[masque1, "distance"] <- NA

  # Update sizes after fusion with cluster 1
  if ("fusionSize" %in% colnames(liste_contiguite))
  {
    if (!is.null(clustersSizes))
    {
      liste_contiguite[masque1, "fusionSize"] <-
      liste_contiguite[masque1, "fusionSize"] +
      clustersSizes[as.character(cluster2)]
    }
    else
      liste_contiguite[masque1, "fusionSize"] <- NA
  }


  masque2A <- liste_contiguite[, 1L] == cluster2
  voisinsCluster2A <- liste_contiguite[masque2A, 2L]
  posVoisinsCluster2A <- which(masque2A)
  masque2B <- liste_contiguite[, 2L] == cluster2
  voisinsCluster2B <- liste_contiguite[masque2B, 1L]
  posVoisinsCluster2B <- which(masque2B)

  posVoisinsCluster2 <- c(posVoisinsCluster2A, posVoisinsCluster2B)

  voisinsCluster2 <- c(voisinsCluster2A, voisinsCluster2B)
  nouvellesConnexions <- setdiff(voisinsCluster2, c(voisinsCluster1, cluster1))


  posSuppression <- posVoisinsCluster2
  if (length(nouvellesConnexions) > 0L)
  {
    masqueChangement <- voisinsCluster2 %in% nouvellesConnexions
    posChangement <- posVoisinsCluster2[masqueChangement]
    posSuppression <- posVoisinsCluster2[!masqueChangement]
    liste_contiguite[posChangement, 1L] <- pmin(cluster1, nouvellesConnexions)
    liste_contiguite[posChangement, 2L] <- pmax(cluster1, nouvellesConnexions)

    # Reset distances
    if ("distance" %in% colnames(liste_contiguite))
      liste_contiguite[posChangement, "distance"] <- NA_real_

    # Update sizes after fusion
    if ("fusionSize" %in% colnames(liste_contiguite))
    {
      if (!is.null(clustersSizes))
      {
        liste_contiguite[posChangement, "fusionSize"] <-
        liste_contiguite[posChangement, "fusionSize"] +
        clustersSizes[as.character(cluster1)]
      }
      else
        liste_contiguite[posChangement, "fusionSize"] <- NA_real_
    }

  }

  # Remove contiguities to cluster 2
  if (length(posSuppression) > 0L)
    liste_contiguite <-
    liste_contiguite[!(seq_len(nrow(liste_contiguite)) %in% posSuppression),
                     , drop = FALSE]

  liste_contiguite
}


#' Respect contiguity list size constraint.
#'
#' This function filters the contiguity list to respect the size constraint.
#'
#' @param pb An object representing the partition problem.
#' @param liste_contiguite A matrix representing the contiguity list.
#' @param taillesClusters A numeric vector representing the sizes
#' of the clusters.
#'
#' @returns A matrix representing the filtered contiguity list.
#' @noRd
size_respect_contiguities <- function(pb, liste_contiguite,
                                      taillesClusters)

{
  if (!pb$hasMaxConstraint() || nrow(liste_contiguite) == 0L)
    return(liste_contiguite)

  if ("fusionSize" %in% colnames(liste_contiguite))
  {
    liste_contiguite <-
      size_after_fusion(liste_contiguite, taillesClusters)
    masqueConservation <- liste_contiguite[, "fusionSize"] <= pb$M

  }
  else
  {
    taillesClustersApresFusion <-
      size_after_fusion(liste_contiguite, taillesClusters)

    masqueConservation <- taillesClustersApresFusion <= pb$M
  }


  liste_contiguite[masqueConservation, , drop = FALSE]
}

codeGroupesLiaison <- function(code) ifelse(is.numeric(code), code,
                                            switch(tolower(code),
                                                   singletonmin = 1L,
                                                   singletonany = 2L,
                                                   minmin = 3L,
                                                   minany = 4L,
                                                   5L))

#' Transform contiguity list.
#'
#' This function performs various transformations on the contiguity list
#' based on specific criteria.
#'
#' @param pb An object of type `pbCon`.
#' @param liste_contiguite A matrix representing the contiguity list.
#' @param taillesClusters A numeric vector representing the
#' sizes of the clusters.
#' @param nbElementsClusters A numeric vector representing the
#' number of elements in each cluster.
#' @param priorisation A logical value indicating whether
#' to perform prioritization.
#' @param alea A logical value indicating whether to apply
#' randomization.
#' @param reduction A character string specifying the
#' type of reduction to apply.
#' @param reductionPartielle A logical value indicating
#' whether to allow partial reduction.
#' @param retour A character string specifying the type of result to return.
#' @param byName A logical value indicating whether to use cluster names.
#' @param orderByDistance A logical value indicating whether
#' to order by distance.
#'
#' @returns A matrix representing the transformed contiguity list.
#' @noRd
transform_contiguity_list <- function(pb, liste_contiguite, # nolint: cyclocomp_linter
                                      taillesClusters,
                                      nbElementsClusters,
                                      priorisation = FALSE,
                                      alea = FALSE, reduction = "non",
                                      reductionPartielle = TRUE,
                                      retour = "liste",
                                      byName = TRUE,
                                      orderByDistance = TRUE)
{
  nbContiguites <- nrow(liste_contiguite)
  posElements <- seq_len(nrow(liste_contiguite))

  if (orderByDistance && !("distance" %in% colnames(liste_contiguite)))
  {
    warning("It is requested to order the contiguity list by distance, but distances are not indicated") # nolint: line_length_linter
    orderByDistance <- FALSE
  }


  if (orderByDistance)
    alea <- FALSE

  # Case where only ordering is requested
  if (!priorisation && reduction == "non")
  {
    if (alea)
      posElements <- sample(seq_len(nbContiguites))
    else if (orderByDistance)
      posElements <- order(liste_contiguite[, "distance"])

    if (retour == "liste")
      return(liste_contiguite[posElements, , drop = FALSE])
    else
      return(posElements)
  }

  if (is.null(names(nbElementsClusters)) || is.null(names(taillesClusters)))
    byName <- FALSE

  if (byName)
  {
    col1 <- as.character(liste_contiguite[, 1L])
    col2 <- as.character(liste_contiguite[, 2L])
  }
  else
  {
    col1 <- liste_contiguite[, 1L]
    col2 <- liste_contiguite[, 2L]
  }


  # Categorization of fusions
  if (pb$hasMinConstraint())
  {
    singletonMin <- (nbElementsClusters[col1] == 1L &
                       taillesClusters[col2] < pb$m) |
      (nbElementsClusters[col2] == 1L &
         taillesClusters[col1] < pb$m)

    singletonAny <- (nbElementsClusters[col1] == 1L |
                       nbElementsClusters[col2] == 1L)

    minMin <- taillesClusters[col1] < pb$m &
      taillesClusters[col2] < pb$m &
      !singletonMin

    minAny <- (taillesClusters[col1] < pb$m |
                 taillesClusters[col2] < pb$m) &
      !singletonMin
  }
  else
  {
    singletonMin <- minAny <- minMin <- rep(FALSE, nrow(liste_contiguite))
    singletonAny <- (nbElementsClusters[col1] == 1L |
                       nbElementsClusters[col2] == 1L)


  }

  restant <- !singletonAny & !minAny


  dataType_mat <- rbind(singletonMin,
                        singletonAny,
                        minMin,
                        minAny,
                        restant)

  # Application of reduction (partial)
  if (!(reduction %in% c("non", codeGroupesLiaison("default"))))
  {
    nbLignesData <- nrow(dataType_mat)
    puissanceReduction <- codeGroupesLiaison(reduction)
    reductionTerminee <- FALSE
    while (!reductionTerminee)
    {
      if (puissanceReduction >= nbLignesData)
        reductionTerminee <- TRUE

      else if (reductionPartielle && !any(dataType_mat[puissanceReduction, ]))
        puissanceReduction <- puissanceReduction + 1L

      else
      {
        reductionTerminee <- TRUE
        dataType_mat[(puissanceReduction + 1L):nbLignesData, ] <- FALSE
      }
    }
  }

  if (!priorisation)
  {
    if (puissanceReduction == 1L)
      posElements <- which(dataType_mat[1L, ])

    else if (puissanceReduction < nbLignesData)
      posElements <- which(dataType_mat[puissanceReduction, ])

    if (alea)
      posElements <- sample(posElements)
    else if (orderByDistance)
      posElements <- order(liste_contiguite[, "distance"])

    if (retour == "liste")
      return(liste_contiguite[posElements, , drop = FALSE])

    else
      return(posElements)
  }

  posElements <- apply(dataType_mat, 1L, which)

  if (alea)
    posElements <- lapply(posElements, sample)

  else if (orderByDistance)
    posElements <- lapply(posElements, sort)

  posElements <- do.call("c", posElements)

  if (retour == "liste")
    return(liste_contiguite[posElements, , drop = FALSE])

  else
    return(posElements)
}

update_contiguities <- function(oldDistances, cluster1, cluster2)
{
  usefulClusters <- seq_len(nrow(oldDistances))[-c(cluster1, cluster2)]

  if (length(usefulClusters) == 0L)
    return(integer(0L))

  usefulClusters <-
    usefulClusters[!is.nan(oldDistances[cluster1, usefulClusters])]

  if (length(usefulClusters) == 0L)
    return(integer(0L))

  usefulClusters[is.infinite(oldDistances[cluster1, usefulClusters]) &
                   is.finite(oldDistances[cluster2, usefulClusters])]
}
