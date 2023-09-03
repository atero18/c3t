#' Generate initial partition.
#'
#' This function generates an initial partition considering
#' connectivity and maximum size constraints.
#'
#' @param pb A `pbCon` object.
#' @param mode A character string specifying the mode of generating
#' the initial partition. Possible values are
#' "unitary" (default) or "random".
#' @param nbFusions An integer specifying the number of fusions
#' to be performed in "random" mode.
#'
#' @returns A numeric vector representing the initial partition.
#' @keywords internal
gen_initial_partition <- function(pb, mode = "unitary", nbFusions = 1L)
{
  n <- nrow(pb)

  partition <- seq_len(n)
  clustersSizes <- clusters_sizes(partition, pb$sizes)

  # If "unitary" mode, return the initial partition without any fusion
  if (mode == "unitary")
    return(partition)

  # If "random" mode, perform random fusions until the specified
  # number of fusions (nbFusions) is reached
  else if (mode == "random")
  {

    # Generate the list of contiguous classes
    clusters_contiguity_list <- clusters_contiguity_list(partition,
                                                         pb$contiguity)

    # Ensure that the contiguous list respects the size constraint
    clusters_contiguity_list <-
      size_respect_contiguities(pb,
                                clusters_contiguity_list,
                                clustersSizes)
    k <- 1L

    while (k <= nbFusions && nrow(clusters_contiguity_list) > 0L)
    {
      # Randomly select a pair of clusters to merge
      l <- sample(seq_len(nrow(clusters_contiguity_list)), size = 1L)

      cluster1 <- clusters_contiguity_list[l, 1L]
      nomCluster1 <- as.character(cluster1)
      cluster2 <- clusters_contiguity_list[l, 2L]
      nomCluster2 <- as.character(cluster2)

      # Get the sizes of the selected clusters.
      pos1 <- names(clustersSizes) == nomCluster1
      tailleCluster1 <- clustersSizes[pos1]
      pos2 <- names(clustersSizes) == nomCluster2
      tailleCluster2 <- clustersSizes[pos2]


      # Update the partition and sizes after the fusion
      partition <- update_partition(partition,
                                    cluster1,
                                    cluster2)

      clustersSizes[pos1] <- tailleCluster1 + tailleCluster2
      clustersSizes <- clustersSizes[!pos2]
      clusters_contiguity_list <-
        update_contiguity_list(clusters_contiguity_list,
                               cluster1, cluster2)

      # Ensure the contiguous list respects the size constraint
      # after the fusion
      clusters_contiguity_list <-
        size_respect_contiguities(pb, clusters_contiguity_list,
                                  clustersSizes)


      k <- k + 1L
    }

    if (k < nbFusions)
    {
      warning("Impossible to complete the requested number of fusions")
    }

    return(standardize_partition(partition))
  }
  else
  {
    stop("Unknown `mode`. Must be unitary or random")
  }
}

#' @describeIn gen_initial_partition generates multiple
#' initial partitions
#' @param modes A character vector specifying the modes of generating
#' the initial partitions.
#' @keywords internal
#' @importFrom purrr partial
#' @importFrom checkmate assertSubset
gen_initial_partitions <- function(pb, modes, nbFusions = 1L)
{
  assertSubset(modes, c("unitary", "random"), empty.ok = FALSE)
  posUnitaire <- which(modes == "unitary")
  if (length(posUnitaire) > 1L)
    modes <- modes[-posUnitaire[-1L]]

  res <-
    vapply(modes, partial(gen_initial_partition,
                          pb = pb,
                          nbFusions = nbFusions),
           integer(nrow(pb)))

  res <- apply(unique(res, MARGIN = 2L), MARGIN = 2L,
               identity, simplify = FALSE)

  names(res) <- unicite_labels(modes)

  res
}
