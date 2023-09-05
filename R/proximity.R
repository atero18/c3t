neutralContiguityNN <- matrix(TRUE, nrow = 0L, ncol = 0L)

#' @title Find the nearest neighbor of some points
#' @description For a set of points (`subsetPoints`), find the nearest
#' neighbor in another set (`subsetNeighbors`). Distances are given in a
#' `p' x n'` matrix (`distances`), with `p' >= p` (p the number of points)
#' and `n' >= n` (n the number of neighbors). If points and neighbors come from
#' the same set, put inner to `TRUE` for assuring the nearest neighbor of a
#' point will not be itself.
#' @param distances distances between points and neighbors. A row is for a
#' point and a column a neighbor. If a vector of length `k > 0` is given,
#' it's consider like a matrix of distances with `p' = k` points and `n' = 1`
#' neighbor.
#' @param subsetPoints A vector of strictly positive integers giving the
#' position of the points the nearest neighbor must be found. Default is
#' all points (`p = p'`). If `distances` has row names can be those instead
#' (as characters). If no points are given, return an empty vector.
#' Values can be duplicated but it will increase complexity.
#' @param subsetNeighbors A vector of strictly positive integers giving the
#' position of the elements of the set that can be consider has a neighbor
#' for each point of `subsetPoints`. Default is all neighbors (`n = n'`).
#' If `distances` has column names can be those instead (as characters).
#' If no neighbors are given an error is generated. Values are expected
#' to be unique.
#' @param inner Flag indicating if points and neighbors come from the same
#' set and `distance` is the set matrix distance. Diagonal of `distance`
#' will then not be considered, assuring that the nearest neighbor of a point is
#' not itself. `FALSE` by default (different sets).
#' @param contiguity In the case there are connectivity constraint
#' between points. If an element of `subsetPoints` and one of
#' `subsetNeighbors` are not contiguous then the second element cannot be a
#' neighbor of the first one. `NULL` by default (no connectivity constraint).
#' If precised, must be a logical matrix with the same dimensions
#' as `distances`.
#' @returns a vector of length `p` with for each point the position of its
#' nearest neighbor (`NaN` if it doesn't exist). A point do not have a nearest
#' neighbor if each distance with its neighbors are missing. If `distances` has
#' column names position of neighbors will be replaced by their names. The
#' returned vector is named with the IDs of the points, or their name if
#' `distances` has row names.
#' @name nearest_neighbor
#' @export
#' @importFrom checkmate assertMatrix assertVector assertIntegerish
#' @importFrom checkmate assertFlag
nearest_neighbor <- function(distances, # nolint: cyclocomp_linter
                             inner = FALSE,
                             contiguity = NULL,
                             subsetPoints = seq_len(nrow(distances)),
                             subsetNeighbors = seq_len(ncol(distances)))
{

  # Ckeck arguments

  if (is.vector(distances))
    distances <- matrix(distances, nrow = 1L)

  if (length(subsetPoints) == 0L)
    return(integer(0L))

  assertMatrix(distances,
               mode = "numeric",
               any.missing = TRUE, all.missing = FALSE)

  if (length(subsetNeighbors) == 0L || ncol(distances) == 0L)
  {
    stop("no neighbors are given")
  }

  if (is.character(subsetPoints))
  {
    assertVector(subsetPoints, any.missing = FALSE, all.missing = FALSE)
    subsetPoints <- match(subsetPoints, rownames(distances))
    if (anyNA(subsetPoints))
    {
      stop("`subsetPoints` make references to `distances` row names but at least one of them doesn't exist.") # nolint: line_length_linter
    }
  }
  else
  {
    assertIntegerish(subsetPoints,
                     lower = 1L,
                     upper = nrow(distances),
                     any.missing = FALSE,
                     all.missing = FALSE)
  }

  if (is.character(subsetNeighbors))
  {
    assertVector(subsetNeighbors, any.missing = FALSE, all.missing = FALSE)
    subsetNeighbors <- match(subsetNeighbors, colnames(distances))
    if (anyNA(subsetNeighbors))
    {
      stop("`subsetNeighbors` make references to `distances` column names but at least one of them doesn't exist.") # nolint: line_length_linter
    }

    if (anyDuplicated(subsetNeighbors))
    {
      stop("Elements of `subsetNeighbors` are supposed to be unique")
    }
  }
  else
  {
    assertIntegerish(subsetNeighbors,
                     lower = 1L,
                     upper = ncol(distances),
                     unique = TRUE,
                     any.missing = FALSE,
                     all.missing = FALSE)
  }

  if (is_igraph(contiguity))
  {
      contiguity <- graph_to_contiguity_matrix(contiguity)
  }
  if (!is.null(contiguity))
    assertMatrix(contiguity, mode = "logical",
                 any.missing = FALSE,
                 all.missing = FALSE,
                 nrows = nrow(distances),
                 ncols = ncol(distances))
  else
    contiguity <- neutralContiguityNN

  assertFlag(inner)

  eval(body(.nearest_neighbor), envir = environment())

  if (!is.null(colnames(distances)))
    nns[] <- colnames(distances)[nns]

  if (!is.null(rownames(distances)))
    names(nns) <- rownames(distances)[as.integer(names(nns))]

  nns
}

#' @importFrom igraph is_igraph
.nearest_neighbor <- function(distances,
                              inner = FALSE,
                              contiguity = neutralContiguityNN,
                              subsetPoints = seq_len(nrow(distances)),
                              subsetNeighbors = seq_len(ncol(distances)))
{
  nbPoints <- length(subsetPoints)
  nbNeighbors <- length(subsetNeighbors)

  if (inner && nrow(distances) != ncol(distances))
  {
    stop("Looking for inner neighbors needs a square distance matrix")
  }

  if (nbNeighbors == 1L)
    nns <- rep(subsetNeighbors, nbPoints)
  else
    nns <- nearest_neighbor_matrix(distances, subsetPoints, subsetNeighbors,
                                   contiguity, inner)

  mask <- nns == 0L
  if (any(mask))
    nns[mask] <- NaN

  names(nns) <- subsetPoints

  nns
}

#' @importFrom checkmate checkNull checkString checkFunction
#' @importFrom checkmate assert assertInt
update_nearest_neighbor <- function(newDistances,
                                    cluster1,
                                    cluster2,
                                    oldNearestNeighbor,
                                    linkage = NULL,
                                    newContiguities = NULL)
{
  assertInt(cluster1, lower = 1L, upper = nrow(newDistances))
  assertInt(cluster2, lower = 1L, upper = nrow(newDistances))

  if (cluster1 == cluster2)
  {
    stop("`cluster1` and `cluster2` must be different")
  }

  assert(checkNull(linkage), checkString(linkage), checkFunction(linkage))

  if (!is.null(newContiguities) && length(newContiguities) == 0L)
    newContiguities <- NULL

  availableClusters <- seq_len(nrow(newDistances)) # nolint: object_usage_linter

  eval(body(.update_nearest_neighbor), envir = environment())
}

.update_nearest_neighbor <- function(newDistances,
                                     cluster1,
                                     cluster2,
                                     oldNearestNeighbor,
                                     linkage = NULL,
                                     newContiguities = NULL,
                                     availableClusters =
                                       which(!is.nan(oldNearestNeighbor)))
{
  # Note : nn = nearest neighbor

  oldNearestNeighbor[cluster2] <- NaN

  diag(newDistances) <- NaN

  availableClusters <- availableClusters[availableClusters != cluster2]
  toUpdate <- availableClusters


  newNearestNeighbor <- oldNearestNeighbor


  if (length(toUpdate) <= 1L)
  {
    newNearestNeighbor[cluster1] <- NaN
    return(newNearestNeighbor)
  }


  # Case when some clusters weren't contiguous with cluster1 but with cluster2.
  # They are now contiguous with cluster1+cluster2 and their nearest neighbor
  # will be in general their actual or cluster1 if their nn weren't cluster2
  # or anyone otherwise.

  # Special cases
  # -- SLINK (single linkage) adapted for contiguity
  if (is_single_linkage(linkage))
  {
    # Cluster different from cluster1 that had cluster1 or cluster2 as
    # nearest neighbor now have cluster1 (if the maximum size constraint is
    # still respected)
    oldNearestNeighbor[cluster1] <- -1L
    nnKeptc1c2 <-
      toUpdate[oldNearestNeighbor[toUpdate] %in% c(cluster1, cluster2)]
    nnKeptc1c2 <- nnKeptc1c2[!is.nan(newDistances[cluster1, nnKeptc1c2])]
    oldNearestNeighbor[cluster1] <- cluster2

    if (length(nnKeptc1c2) > 0L)
    {
      newNearestNeighbor[nnKeptc1c2] <- cluster1
      toUpdate <- setdiff(toUpdate, nnKeptc1c2)
    }

    rm(nnKeptc1c2)

    # Cluster different from cluster1 that didn't have cluster1 and cluster2
    # as nearest neighbor and don't have a new contiguity with cluster1+cluster2
    # keep their nearest neighbor
    if (length(toUpdate > 1L))
    {
      # if there is no new contiguity (i.e. newContiguities is NULL) then
      # noNewContiguities = toUpdate # nolint: commented_code_linter
      noNewContiguities <- setdiff(toUpdate, newContiguities)

      nnKeptDifferent <-
        noNewContiguities[!(
          oldNearestNeighbor[noNewContiguities] %in% c(cluster1, cluster2))] # nolint: indendation_linter

      if (length(nnKeptDifferent) > 0L)
        toUpdate <- setdiff(toUpdate, nnKeptDifferent)

      rm(noNewContiguities, nnKeptDifferent)
    }
  }

  # -- CLINK (complete linkage) adapted for contiguity
  else if (is_complete_linkage(linkage))
  {
    # Clusters that had a nearest neighbor different from cluster1 and cluster2
    # keep their neighbor. Can be applied on clusters with new contiguity
    # with cluster1+cluster2 because their distance with cluster1+cluster2
    # will be higher or equal than their distance with cluster2 which is not
    # their lowest.
    nnKeptDifferent <-
      toUpdate[!(oldNearestNeighbor[toUpdate] %in% c(cluster1, cluster2))]

    if (length(nnKeptDifferent > 0L))
      toUpdate <- setdiff(toUpdate, nnKeptDifferent)

    rm(nnKeptDifferent)
  }

  # Clusters where distances must be compared with all available clusters.
  # In general when the nearest neighbor was cluster1 or cluster2.
  # For single and complete linkage the only cluster neededing a
  # full comparaison is cluster1+cluster2
  if (is_single_linkage(linkage))
    needCompleteCompare <- cluster1

  else
  {
    needCompleteCompare <-
      toUpdate[oldNearestNeighbor[toUpdate] %in% c(cluster1, cluster2)]
  }

  newNearestNeighbor[needCompleteCompare] <-
    .nearest_neighbor(newDistances, inner = TRUE,
                      subsetPoints = needCompleteCompare,
                      subsetNeighbors = availableClusters)


  # If all nn are updated
  if (length(toUpdate) == length(needCompleteCompare))
    return(newNearestNeighbor)

  toUpdate <- setdiff(toUpdate, needCompleteCompare)

  # Clusters where the check for nearest needs only a comparison between
  # the previous neighbor and cluster1+cluster2. In general cluster that
  # had a nearest neighbor different from cluster1 and cluster2. Especially
  # needed for cluster which have a new contiguity with cluster1+cluster2
  needUniqueCompare <- toUpdate

  # Positions of the distances between needUniCompare clusters and
  # cluster1+cluster2 in `newDistances`
  indexs <- cbind(oldNearestNeighbor[needUniqueCompare],
                  needUniqueCompare)

  maskChange <-
    newDistances[cluster1, needUniqueCompare] < newDistances[indexs]

  # If cluster1+cluster2 is near than the previous nn
  if (any(maskChange))
    newNearestNeighbor[needUniqueCompare[maskChange]] <- cluster1

  newNearestNeighbor
}
