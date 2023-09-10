#' @include pbCon.R
NULL


#' Find articulation points
#'
#' @description find articulation points of a partition.
#' @param partition A partition. For `is_articulation_pt_cluster`,
#' must be a regionalisation.
#' (vector of strictly positive integers)
#' @param contiguity A contiguity matrix or a contiguity graph
#' (package `igraph`). (contiguity matrix or `igraph` graph)
#' @name articulation_points
NULL


# Check is an alement is an articulation point ----------------------------

#' @describeIn articulation_points Indicates whether a point in a
#' set is an articulation point or not. Based on a BFS
#' traversal, optimized as it can stop earlier if `x` is not
#' an articulation point.
#'
#' @param x The point for which to determine whether it is an articulation
#' point or not. Can be an identifier (strictly positive integer) or a name,
#' if `contiguity` is an `igraph` (strictly positive integer or string).
#' @returns `TRUE` if `x` is an articulation point for the connected set,
#' `FALSE` otherwise. (boolean)
#' @noRd
#' @keywords internal
#' @importFrom igraph is_igraph V neighbors components
#' @importFrom igraph is_named delete_vertices induced_subgraph
#' @importFrom igraph vertex_attr<-
is_articulation_pt <- function(x, contiguity) # nolint: cyclocomp_linter
{

  # If the contiguity is not an igraph object
  if (!is_igraph(contiguity))
  {
    if (nrow(contiguity) <= 2L)
      return(FALSE)

    contiguity <- contiguity_matrix_to_graph(contiguity)
  }

  # If contiguity graph is not named, assign names
  if (!is_named(contiguity))
  {
    vertex_attr(contiguity) <- list(name = seq_along(V(contiguity)))
  }


  n <- length(V(contiguity))

  # If number of vertices is less than or equal to 2,
  # there is no articulation point
  if (n <= 2L)
    return(FALSE)

  # Find neighbors of x in the graph
  neighborsX <- names(neighbors(contiguity, x))
  nbNeighborsX <- length(neighborsX)

  # If x has no neighbor
  if (nbNeighborsX == 0L)
    return(FALSE)

  # If x has only one neighbor, it's not an articulation point
  if (nbNeighborsX == 1L)
    return(FALSE)

  # Create a subgraph induced by neighbors of x
  subgraphNeighborsX <- induced_subgraph(contiguity, neighborsX)

  # Identify connected components of the subgraph
  connectedCompNeighborsX <-
    components(subgraphNeighborsX, mode = "weak")$membership
  nbCCNeighborsX <- length(unique(connectedCompNeighborsX))

  # If the subgraph is connected, x is not an articulation point
  if (nbCCNeighborsX == 1L)
    return(FALSE)

  # Create a list of connected components of neighbors of x
  listCCNeighborsX <-
    lapply(seq_len(nbCCNeighborsX),
           function(c) neighborsX[connectedCompNeighborsX == c])

  # On va chercher si tous les voisins de x communiquent dans
  # la région une fois x supprimé

  # Create a graph without x
  graphWithoutX <- delete_vertices(contiguity, x)

  # Choose a connected component at random
  c <- sample(seq_len(nbCCNeighborsX), size = 1L)
  listNeighbors <- c
  sizeListNeighbors <- 1L

  # Choose a neighbor of x from this connedted component at random
  if (length(listCCNeighborsX[[c]]) > 1L)
  {
    stack <- sample(listCCNeighborsX[[c]], size = 1L)
  } else
    stack <- listCCNeighborsX[[c]]


  # Create an empty list for marked vertices
  markedVertices <- integer(0L)

  while (length(stack) > 0L &&
        sizeListNeighbors != nbCCNeighborsX)
  {
    # Take the first vertex from the stack
    p <- stack[1L]

    # Find neighbors of p that have not been marked
    V <- setdiff(names(neighbors(graphWithoutX, p)), markedVertices)

    markedVertices <- c(markedVertices, V)

    if (length(V) != 0L)
    {
      # If unmarked neighbors of x have been found we add all their connected
      # component(s) to the list
      for (vx in intersect(V, neighborsX))
      {
        connectedCompVX <- connectedCompNeighborsX[vx]
        if (!(connectedCompVX %in% listNeighbors))
        {
          listNeighbors <- c(listNeighbors, connectedCompVX)
          sizeListNeighbors <- sizeListNeighbors + 1L
        }
      }

      stack <- c(V, stack[-1L])
    }
    else
      stack <- stack[-1L]
  }

  # If all connected components of neighbors of x are observed,
  # it's not an articulation point
  if (sizeListNeighbors == nbCCNeighborsX)
    return(FALSE)

  else
    return(TRUE)

}

#' @describeIn articulation_points Indicates whether a point
#' in a region resulting from a regionalisation is an articulation point or not.
#' @importFrom igraph is_igraph induced_subgraph
#' @noRd
#' @keywords internal
is_articulation_pt_cluster <- function(x, contiguity, partition)
{

  regionX <- partition[x]
  elementsRegionX <- which(partition == regionX)

  if (length(elementsRegionX) == 1L)
    return(FALSE)

  if (is_igraph(contiguity))
  {
    contiguityRegionX <- induced_subgraph(contiguity, elementsRegionX)
    return(is_articulation_pt(x, contiguityRegionX))
  }

  else
  {
    contiguityMatrixRegionX <- contiguity[elementsRegionX, elementsRegionX]
    posX <- which(elementsRegionX == x)
    return(is_articulation_pt(posX, contiguityMatrixRegionX))
  }
}


# Find all articulation points --------------------------------------------


#' @describeIn articulation_points Returns the list of articulation points
#' in a cluster resulting from a partition.
#' @param cluster The cluster identifier, not empty for `partition`.
#' (strictly positive integer)
#' @importFrom igraph induced_subgraph articulation_points
#' @returns A vector composed of indices of elements that are articulation
#' points. (vector of strictly positive integers)
#' @keywords internal
articulation_pts_cluster <- function(contiguity,
                                     partition,
                                     cluster)
{
  clusterPoints <- which(partition == cluster)
  if (length(clusterPoints) <= 2L)
    return(integer(0L))

  if (is_pbCon(contiguity))
    contiguityCluster <- induced_subgraph(contiguity$contiguity,
                                          clusterPoints)

  else
    contiguityCluster <- induced_subgraph(contiguity, clusterPoints)

  clusterPoints[articulation_points(contiguityCluster)]
}

#' @describeIn articulation_points Returns the list of points that are
#' articulation points for different clusters.
#' @param clusters List of non-empty clusters for `partition` for which to
#' determine their articulation point(s). Default is all clusters.
#' @returns For `articulation_pts_clusters` a vector composed of indices of
#' elements from clusters
#' in `clusters` that are articulation points.
#' (vector of strictly positive integers)
#' @seealso [igraph::articulation_points()]
#' @importFrom purrr partial
#' @export
articulation_pts_clusters <- function(contiguity,
                                      partition,
                                      clusters = clustersIDs(partition))
{
  # For each cluster of `clusters` we look for its articulation points.
  pointsArticulationClusters <-
    lapply(clusters,
           partial(articulation_pts_cluster,
                   partition = partition,
                   contiguity = contiguity))

  do.call("c", pointsArticulationClusters)
}


#' @describeIn articulation_points Updates the set of
#' articulation points for a modified partition.
#' @param oldPartition The initial partition.
#' (vector of strictly positive integers)
#' @param newPartition The modified partition.
#' (vector of strictly positive integers)
#' @param oldArticulationPoints Set of articulation points for the
#' initial partition (vector of strictly positive integers)
#' @param modifiedClusters Set of clusters that have been modified.
#' If not specified, will be calculated. (vector of strictly positive integers)
#' @returns For `update_art_pts_partition` the set of
#' articulation points for the modified partition.
#' (vector of strictly positive integers)
#' @importFrom igraph induced_subgraph articulation_points
#' @keywords internal
#' @noRd
update_art_pts_partition <- function(contiguity, oldPartition,
                                     newPartition,
                                     oldArticulationPoints,
                                     modifiedClusters = NULL)
{
  # Search for clusters that have been modified
  # If donorCluster and receiverCluster are provided by the user,
  # we assume that only these clusters have been modified
  if (is.null(modifiedClusters))
  {
    mask <- oldPartition != newPartition
    modifiedClusters <-
      unique(c(oldPartition[mask], newPartition[mask]))
  }

  if (length(modifiedClusters) == 0L)
    return(oldArticulationPoints)

  # Preserve the articulation points of unchanged clusters
  keptArticulationPts <-
    setdiff(oldArticulationPoints, which(oldPartition %in% modifiedClusters))

  # Search for articulation points in the modified clusters
  newArticulationPts <- articulation_pts_clusters(contiguity,
                                                  newPartition,
                                                  modifiedClusters)


  return(c(keptArticulationPts, newArticulationPts))
}
