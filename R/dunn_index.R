#' Dunn Index
#'
#' Calculates the Dunn index for a given partition based on
#' the inter-element distance matrix. The Dunn index assesses the compactness
#' and separation of clusters in a partition.
#'
#' @param distances The inter-element distance matrix.
#' @param partition The partition of interest.
#' @param linkage The method for computing inter-cluster distances.
#'                Defaults to "single" (minimum distance between clusters).
#' @param valueOnly If TRUE, returns only the Dunn index value (D)
#' instead of the full result.
#'
#' @returns If `valueOnly` is TRUE, returns the Dunn index value (D) as a
#' floating-point number. Otherwise, returns a list containing the
#' following components:
#'           - `diameters`: A vector containing the diameters of each cluster.
#'           - `linkageDistances`: A matrix containing the inter-cluster
#'           distances.
#'           - `clustersMoinsSepares`: The IDs of the two clusters with
#'           the smallest separation.
#'           - `separation`: The separation between the two clusters with
#'           the smallest separation.
#'           - `leastCompactCluster`: The ID of the least compact cluster.
#'           - `compactness`: The diameter of the least compact cluster.
#'           - `D`: The Dunn index value, calculated as the separation divided
#'           by the compactness.
#'           - `linkage`: The method used for computing inter-cluster distances.
#'
#' @importFrom purrr partial
#' @name dunn
#' @references Dunn, Joseph C. (1973) A fuzzy relative of the ISODATA
#' process and its use in detecting compact well-separated clusters,
#' Journal of Cybernetics
#'
#' @keywords internal
dunn_index <- function(distances, partition, linkage = "single",
                       valueOnly = FALSE, ...)
{

  # Convert the partition to a list of clusters
  listeClusters <- partition_to_list(partition)
  nbClusters <- length(listeClusters)
  nomsClusters <- names(listeClusters)
  clustersIDs <- as.integer(nomsClusters)

  # Compactness
  # -- Calculation of the diameter for each cluster
  calculate_diameter <-
    function(indices) diameter_cluster(distances = distances,
                                       indices = indices)
  diameters <- c3t_sapply(listeClusters, calculate_diameter)

  # -- Finding the least compact cluster
  leastCompactCluster <- which.max(diameters)
  compactness <- unname(diameters[leastCompactCluster])
  leastCompactCluster <- clustersIDs[leastCompactCluster]

  # Separation
  # -- Creating a grid to store all cluster fusion combinations
  grid <- cbind(rep(clustersIDs, each  = nbClusters),
                rep(clustersIDs, times = nbClusters))

  grid <- grid[grid[, 1L] < grid[, 2L], , drop = FALSE]

  # -- Calculating inter-cluster distances
  linkageDistances_vec <- .calcul_distances_inter(distances,
                                                partition,
                                                grid, linkage)

  linkageDistances_mat <- matrix(Inf, nrow = nbClusters, ncol = nbClusters)
  linkageDistances_mat[lower.tri(linkageDistances_mat)] <-
    linkageDistances_mat[upper.tri(linkageDistances_mat)] <-
    linkageDistances_vec

  rownames(linkageDistances_mat) <-
    colnames(linkageDistances_mat) <- clustersIDs

  # -- Finding the two clusters with the smallest separation
  clustersMoinsSepares <- which.min(linkageDistances_vec)[1L]
  separation <- unname(linkageDistances_vec[clustersMoinsSepares])
  clustersMoinsSepares <- clustersIDs[grid[clustersMoinsSepares, ]]

  if (compactness == 0.0)
    D <- Inf

  else
    D <- separation / compactness

  if (valueOnly)
    return(D)

  list(diameters = diameters, linkageDistances = linkageDistances_mat,
       clustersMoinsSepares = clustersMoinsSepares,
       separation = separation, leastCompactCluster = leastCompactCluster,
       compactness = compactness,
       D = D,
       linkage = linkage)
}

#' Update Dunn index
#'
#' @inheritParams update_criterion_params
#' @keywords internal
update_dunn_index <- function(dataCriterion, distances, partitionBefore,
                              donor, receiver,
                              givenElement)
{
  newPartition <- partitionBefore
  newPartition[givenElement] <- receiver

  listeClusters <- partition_to_list(newPartition)
  nomsClusters <- names(listeClusters)
  clustersIDs <- as.integer(nomsClusters)

  nomClusterDon <- as.character(donor)
  nomClusterRec <- as.character(receiver)

  # Compactness
  # -- Update the diameters of donor and receiver clusters
  diameterDonor  <- max(distances[listeClusters[[nomClusterDon]],
                                  listeClusters[[nomClusterDon]]])
  dataCriterion$diameters[donor] <- diameterDonor

  diameterReceiver  <- max(distances[listeClusters[[nomClusterRec]],
                                     listeClusters[[nomClusterRec]]])
  dataCriterion$diameters[receiver] <- diameterReceiver

  # -- Update the compactness
  # The diameter of the receiver can only increase (>=)
  if (diameterReceiver  >= dataCriterion$compactness)
  {
    dataCriterion$compactness <- diameterReceiver
    dataCriterion$leastCompactCluster <- receiver
  }
  # The diameter of the donor can only decrease (<=)
  else if (dataCriterion$leastCompactCluster == donor)
  {
    posMax <- which.max(dataCriterion$diameters)
    dataCriterion$compactness <- dataCriterion$diameters[posMax]
    dataCriterion$leastCompactCluster <- clustersIDs[posMax]
  }

  # Separation
  # -- Update inter-cluster distances
  # --- With the donor
  clustersDifferentsDonneur <- clustersIDs[clustersIDs != donor]
  nomClustersDifferentsDonneur <- as.character(clustersDifferentsDonneur)
  distancesADonneur <- .calcul_distances_inter(distances,
                                               newPartition,
                                               cbind(donor,
                                                     clustersDifferentsDonneur),
                                               dataCriterion$linkage)

  posMin <- which.min(distancesADonneur)
  minDistDonneur <- distancesADonneur[posMin]
  changement <- FALSE
  if (minDistDonneur < dataCriterion$separation)
  {
    changement <- TRUE
    dataCriterion$separation <- minDistDonneur
    dataCriterion$clustersMoinsSepares <-
      sort(c(donor, clustersDifferentsDonneur[posMin]))
  }

  dataCriterion$linkageDistances[nomClusterDon, nomClustersDifferentsDonneur] <-
    distancesADonneur
  dataCriterion$linkageDistances[nomClustersDifferentsDonneur, nomClusterDon] <-
    distancesADonneur

  # --- With the receiver
  clustersDifferentsReceveur <- clustersIDs[clustersIDs != receiver]
  nomClustersDifferentsReceveur <- as.character(clustersDifferentsReceveur)
  distancesAReceveur <-
    .calcul_distances_inter(distances,
                            newPartition,
                            cbind(receiver,
                                  clustersDifferentsReceveur),
                            dataCriterion$linkage)

  # -- If there was no change and the previous separation involved the receiver
  # or donor cluster, that distance is obsolete
  posMin <- which.min(distancesAReceveur)[1L]
  minDistReceveur <- distancesAReceveur[posMin]
  if (minDistReceveur < dataCriterion$separation)
  {
    changement <- TRUE
    dataCriterion$separation <- minDistReceveur
    dataCriterion$clustersMoinsSepares <-
      sort(c(receiver, nomClustersDifferentsReceveur[posMin]))
  }
  dataCriterion$linkageDistances[nomClusterRec,
                                 nomClustersDifferentsReceveur] <-
    distancesAReceveur
  dataCriterion$linkageDistances[nomClustersDifferentsReceveur,
                                 nomClusterRec] <-
    distancesAReceveur

  # If no change has been ben made and the previous separation were concerning
  # the donor and / or receiver cluster, this distance doesn't exist anymore
  if (!changement &&
      any(c(donor, receiver) %in%
          dataCriterion$clustersMoinsSepares))
  {
    separation <- min(dataCriterion$linkageDistances)
    clustersMoinsSepares <- which(dataCriterion$linkageDistances == separation,
                                  arr.ind = TRUE)[1L, ]
    dataCriterion$separation <- separation
    dataCriterion$clustersMoinsSepares <-
      sort(clustersIDs[clustersMoinsSepares])
  }

  dataCriterion$D <- ifelse(dataCriterion$compactness != 0.0,
                            dataCriterion$separation / dataCriterion$compactness,
                            Inf)

  return(dataCriterion)
}
