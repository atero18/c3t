#' @noRd
frontier_art_to_transferable <- function(frontier_mat,
                                         articulationPoints_vec)
{
  if (nrow(frontier_mat) == 0L)
    return(integer(0L))

  frontierElems_vec <- unique(c(frontier_mat[, "x1"]))
  setdiff(frontierElems_vec, articulationPoints_vec)
}


#' Generates a list of transferable elements for `pbCon` class
#'
#' Determines the list of transferable elements for each cluster in a partition.
#'
#' @param regionalisation The partition (vector of positive integers)
#' @param pb A `pbCon` class object
#' @param contiguityGraph An `igraph` contiguity graph
#' @param contiguityMat A `contiguity` matrix
#' @returns A list with three elements: 'frontier', 'articulation',
#' 'transferable' containing information about boundary, articulation,
#' and transferable elements.
#' @keywords internal
#' @name elements_transferables
#' @importFrom igraph is_named vertex_attr V
transferable_elements_pb <- function(regionalisation, pb,
                                     contiguityGraph, contiguityMat)
{
  clustersList_vec <- unique(regionalisation)

  if (missing(pb) && missing(contiguityGraph) && missing(contiguityMat))
  {
    stop("At least the graph or the contiguity matrix must be provided")
  }

  else if (missing(pb))
  {
    if (missing(contiguityMat))
      contiguityMat <- graph_to_contiguity_matrix(contiguityGraph)

    if (missing(contiguityGraph))
      contiguityGraph <- contiguity_matrix_to_graph(contiguityMat)
  }

  if (missing(pb))
  {
    interiorFront_mat <-
      interior_boundary(contiguityMat, regionalisation, removeSymmetry = FALSE)
  }

  else
  {
    interiorFront_mat <-
      interior_boundary(pb, regionalisation, removeSymmetry = FALSE)
  }


  if (length(interiorFront_mat) == 0L)
    return(interiorFront_mat)

  if (missing(pb))
  {
    articulationPoints_vec <-
      articulation_pts_clusters(contiguityGraph,
                                regionalisation,
                                clustersList_vec)
  }

  else
  {
    articulationPoints_vec <-
      articulation_pts_clusters(pb,
                                regionalisation,
                                clustersList_vec)
  }


  list(frontier = interiorFront_mat,
       articulation = articulationPoints_vec,
       transferable =
         frontier_art_to_transferable(interiorFront_mat,
                                      articulationPoints_vec))

}

#' Updates the list of transferable elements after partition modification
#'
#' Determines the updated list of transferable elements after a modification
#' in the partition.
#'
#' @describeIn elements_transferables Determines the updated list of
#' transferable elements after a modification in the partition.
#' @param anciennePartition_vec The initial partition
#' (vector of positive integers)
#' @param nouvellePartition_vec The modified partition
#' (vector of positive integers)
#' @param oldTransferablePoints_list The previous list of
#' transferable elements
#' @param donor Identifier of the donor cluster
#' @param receiver Identifier of the receiver cluster
#' @keywords internal
maj_elements_transferables <- function(anciennePartition_vec,
                                       nouvellePartition_vec,
                                       oldTransferablePoints_list,
                                       pb, contiguityGraph, contiguityMat,
                                       donor = NULL,
                                       receiver = NULL)
{

  # Finding the clusters that have been modified
  # If clusterDonneur and clusterReceveur are given by the user, we assume only
  # those clusters have been modified
  if (is.null(donor) || is.null(receiver))
  {
    mask <- anciennePartition_vec != nouvellePartition_vec
    clustersModifies <- unique(c(anciennePartition_vec[mask],
                                 nouvellePartition_vec[mask]))
  }
  else
    clustersModifies <- c(donor, receiver)

  # Calculating the new boundary
  interiorFront_mat <- oldTransferablePoints_list$frontier
  mask <- !(interiorFront_mat[, "cluster1"] %in% clustersModifies) &
            !(interiorFront_mat[, "cluster2"] %in% clustersModifies)

  interiorFront_mat <- interiorFront_mat[mask, , drop = FALSE]

  if (missing(pb))
  {
    newInteriorFront_mat <-
      lapply(clustersModifies, function(c)
        interior_boundary_cluster(contiguityMat, c,
                                  nouvellePartition_vec,
                                  removeSymmetry = FALSE))
  }
  else
  {
    newInteriorFront_mat <-
      lapply(clustersModifies,
             partial(interior_boundary_cluster,
                     contiguity = pb,
                     partition  = nouvellePartition_vec,
                     removeSymmetry = FALSE))
  }



  newInteriorFront_mat <- do.call("rbind", newInteriorFront_mat)
  mask <- !(newInteriorFront_mat[, "cluster2"] %in% clustersModifies)
  interiorFront_mat <- rbind(interiorFront_mat,
                             newInteriorFront_mat,
                             newInteriorFront_mat[mask,
                                                  c(2L, 1L, 4L, 3L),
                                                  drop = FALSE])

  # Calculating articulation points
  if (missing(pb))
  {
    pointsArticulationsMAJ_vec <-
      update_art_pts_partition(contiguityGraph,
                               anciennePartition_vec,
                               nouvellePartition_vec,
                               oldTransferablePoints_list$articulation,
                               clustersModifies)
  }
  else
  {
    pointsArticulationsMAJ_vec <-
      update_art_pts_partition(pb, anciennePartition_vec,
                               nouvellePartition_vec,
                               oldTransferablePoints_list$articulation,
                               clustersModifies)
  }

  list(frontier =
         interiorFront_mat, articulation = pointsArticulationsMAJ_vec,
       transferable =
         frontier_art_to_transferable(interiorFront_mat,
                                      pointsArticulationsMAJ_vec))
}

#' @keywords internal
#' @importFrom checkmate assertFlag
remove_transferables <- function(pb, # nolint: cyclocomp_linter
                                 partition,
                                 transferablePoints,
                                 clustersSizes =
                                   clusters_sizes(partition, pb$sizes),
                                 rmFromSingleton = FALSE,
                                 rmCreationMin = TRUE,
                                 rmToMin = FALSE,
                                 rmFromMin = FALSE,
                                 rmCreationMax = TRUE,
                                 rmFromMax = FALSE,
                                 rmToMax = TRUE,
                                 rmNormToNorm = FALSE,
                                 smallClusters = NULL,
                                 bigClusters = NULL,
                                 nbElemsClusters = table(partition))
{

  transferables <- transferablePoints$transferable

  if (length(transferables) == 0L)
    return(transferablePoints)

  assertFlag(rmFromSingleton)

  singletonsClusters <- which(nbElemsClusters == 1L)
  rmFromSingleton <- rmFromSingleton && length(singletonsClusters) > 0L

  assertFlag(rmNormToNorm)

  if (pb$hasMinConstraint())
  {
    assertFlag(rmCreationMin)
    assertFlag(rmFromMin)
    assertFlag(rmToMin)
  }
  else
  {
    rmCreationMin <- FALSE
    rmFromMin <- FALSE
    rmToMin <- FALSE

  }

  if (pb$hasMaxConstraint())
  {
    assertFlag(rmCreationMax)
    assertFlag(rmFromMax)
    assertFlag(rmToMax)
  }
  else
  {
    rmCreationMax <- FALSE
    rmFromMax <- FALSE
    rmToMax <- FALSE
  }

  if (rmFromMin || rmToMin || rmCreationMin || rmNormToNorm)
  {
    if (is.null(smallClusters))
      smallClusters <- which(clustersSizes < pb$m)

    nbSmallClusters <- length(smallClusters)
    rmFromMin <- rmFromMin && nbSmallClusters > 0L
    rmToMin <- rmToMin && nbSmallClusters > 0L
  }


  if (rmFromMax || rmToMax || rmCreationMax || rmNormToNorm)
  {
    if (is.null(bigClusters))
      bigClusters <- which(clustersSizes > pb$M)

    nbBigClusters <- length(bigClusters)
    rmFromMax <- rmFromMax && nbBigClusters > 0L
    rmToMax <- rmToMax && nbBigClusters > 0L
  }

  # Normal to normal
  if (rmNormToNorm && length(transferables) > 0L)
  {

    maskNormToNorm <-
      transferablePoints$frontier[, "x1"] %in% transferables &
      !transferablePoints$frontier[, "cluster1"] %in% smallClusters &
      !transferablePoints$frontier[, "cluster1"] %in% bigClusters &
      !transferablePoints$frontier[, "cluster2"] %in% smallClusters &
      !transferablePoints$frontier[, "cluster2"] %in% bigClusters

    if (any(maskNormToNorm))
    {
      transferablePoints$frontier <-
        transferablePoints$frontier[!maskNormToNorm, , drop = FALSE]
      transferables <-
        transferables[transferables %in% transferablePoints$frontier[, "x1"]]
    }
  }

  # Singletons
  if (rmFromSingleton && length(transferables) > 0L)
  {
    transferables <-
      transferables[!partition[transferables] %in% which(nbElemsClusters == 1L)]
  }

  # Min clusters
  if (rmFromMin && length(transferables) > 0L)
  {
    maskFromMin <- transferables %in% which(partition %in% smallClusters)

    if (!rmFromSingleton && length(singletonsClusters) > 0L)
      maskFromMin <- maskFromMin & !transferables %in% singletonsClusters

    transferables <- transferables[!maskFromMin]
  }

  if (rmToMin && length(transferables) > 0L)
  {
    maskToMin <-
      transferablePoints$frontier[, "x1"] %in% transferables &
      transferablePoints$frontier[, "cluster2"] %in% smallClusters

    if (any(maskToMin))
    {
      transferablePoints$frontier <-
        transferablePoints$frontier[!maskToMin, , drop = FALSE]

      transferables <-
        transferables[transferables %in% transferablePoints$frontier[, "x1"]]

    }
  }

  if (rmCreationMin && length(transferables) > 0L)
  {
    maskCreationMin <-
      clustersSizes[partition[transferables]] >= pb$m &
      clustersSizes[partition[transferables]] - pb$sizes[transferables] < pb$m

    transferables <- transferables[!maskCreationMin]
  }

  # Max clusters
  if (rmFromMax && length(transferables) > 0L)
  {
    maskFromMax <- transferables %in% which(partition %in% bigClusters)
    transferables <- transferables[!maskFromMax]

  }

  if (rmToMax && length(transferables) > 0L)
  {
    maskToMax <-
      transferablePoints$frontier[, "x1"] %in% transferables &
      transferablePoints$frontier[, "cluster2"] %in% bigClusters

    if (any(maskToMax))
    {
      transferablePoints$frontier <-
        transferablePoints$frontier[!maskToMax, , drop = FALSE]
      transferables <-
        transferables[transferables %in% transferablePoints$frontier[, "x1"]]
    }
  }

  if (rmCreationMax && length(transferables) > 0L)
  {
    maskCreationMax <-
      transferablePoints$frontier[, "x1"] %in% transferables &
      !transferablePoints$frontier[, "cluster2"] %in% bigClusters &
      clustersSizes[transferablePoints$frontier[, "cluster2"]] +
      pb$sizes[transferablePoints$frontier[, "x1"]] > pb$M

    if (any(maskCreationMax))
    {
      transferablePoints$frontier <-
        transferablePoints$frontier[!maskCreationMax, , drop = FALSE]

      transferables <-
        transferables[transferables %in% transferablePoints$frontier[, "x1"]]
    }

  }

  transferablePoints$transferable <- transferables
  transferablePoints
}
