#' Number of clusters in a partition (per connected component)
#'
#' Give the number of clusters in a partition. This number can be the
#' total number or the number in each connected component.
#' @param partition a partition vector.
#' @param components a vector giving for each element its connected component.
#' If `NULL` (default) components will be ignored. Otherwise must be a vector
#' of the same length of `partition`.
#' @returns If `components` is `NULL`, the total number of clusters in
#' `partition` (positive integer). Otherwise a named vector of positive integers
#' of length the number of connected components and giving the number of
#' clusters per connected component.
#' @importFrom checkmate assertVector
#' @examples nbClusters(c(1L, 2L, "x", 1L)) # 3
#' @export
nbClusters <- function(partition, components = NULL)
{
  # Checking arguments
  assertVector(partition, min.len = 1L,
               any.missing = FALSE, all.missing = FALSE)

  assertVector(components, len = length(partition),
               any.missing = FALSE, all.missing = FALSE,
               null.ok = TRUE)

  .nbClusters(partition, components)
}

#' @name nbClusters
#' @keywords internal
#' @noRd
.nbClusters <- function(partition, components = NULL)
{
  if (is.null(components))
    length(unique(partition))
  else
  {
    tapply(partition, components, .nbClusters,
           components = NULL, simplify = TRUE)
  }

}

setGeneric(".nbClusters")

#' @name nbClusters
#' @keywords internal
#' @noRd
.nbRegions <- .nbClusters


clustersIDs <- function(partition)
{
  unique(partition)
}

nbSingletons <- function(partition)
{
  as.integer(sum(table(partition) == 1L))
}

setGeneric("nbSingletons")

#' Standardize a partition
#'
#' Changes the cluster identifiers in the partition vector to form a
#' consecutive sequence.
#'
#' @param partition A non-empty vector representing the partition.
#'
#' @returns A strictly positive integer vector representing the updated
#' partition with consecutive cluster identifiers.
#' @export
#' @details This form is also used by [stats::cutree()].
#' @details Two partitions are equivalent if and only if their
#' standardised forms are equal.
#' @example inst/doc/ex_standardize_partition.R
standardize_partition <- function(partition)
{
  clustersIDs <- clustersIDs(partition)
  nbClusters <- length(clustersIDs)

  if (nbClusters == 1L)
    return(rep(1L, length(partition)))

  match(partition, clustersIDs)
}

#' Are two partitions equivalent?
#'
#' Check if two partitions are equivalent, i.e. represent the same clustering.
#' @param partition1,partition2 two non-empty vectors representing.
#' partitions to be compared
#'
#' @returns `TRUE` if both partitions are equivalent, `FALSE` otherwise.
#' @importFrom checkmate assertVector
#' @export
#'
#' @example inst/doc/ex_equivalent_partitions.R
equivalent_partitions <- function(partition1, partition2)
{
  assertVector(partition1, min.len = 1L,
               any.missing = FALSE, all.missing = FALSE)

  assertVector(partition2, min.len = 1L,
               any.missing = FALSE, all.missing = FALSE)

  length(partition1) == length(partition2) &&
    .nbClusters(partition1) == .nbClusters(partition2) &&
    all(standardize_partition(partition1) == standardize_partition(partition2))
}

#' Merges Partition Connected Components
#'
#' This function merges the partitions of connected components back together.
#'
#' @param components A vector indicating the connected component index of
#' each element.
#' @param partitions A list of partition vectors, where each vector represents
#' a partition of a connected component.
#'
#' @returns A merged partition vector.
#'
#' @export
#' @importFrom checkmate assertIntegerish assertList
merge_cc_partitions <- function(components, partitions)
{
  # Ensure components is an integer-like vector with at least one element
  assertIntegerish(components, min.len = 1L)

  # Ensure partitions is a list with the same length as the number of unique
  # components, and each element is a vector
  assertList(partitions, len = length(unique(components)), types = "vector")

  # Check if each partition has the same number of elements as the
  # corresponding connected component
  if (!all(table(components) == lengths(partitions)))
  {
    stop("Each partition should have a number of elements corresponding to the connected components") # nolint: line_length_linter
  }

  # If there is only one partition, return it after changing the IDs
  if (length(partitions) == 1L)
    return(standardize_partition(partitions[[1L]]))

  # Change IDs within each partition vector
  partitions <- lapply(partitions, standardize_partition)
  # Calculate the maximum number of clusters in each partition
  nbClustersPartitions <- vapply(partitions, max, integer(1L))
  # Calculate cumulative counts of clusters in each partition
  cumNbClustersPartitions <- c(0L, cumsum(nbClustersPartitions))
  cumNbClustersPartitions <-
    cumNbClustersPartitions[-length(cumNbClustersPartitions)]
  # Shift partition IDs according to cumulative counts in order to have
  # unique IDs
  partitions <- .mapply("+", dots = list(partitions, cumNbClustersPartitions),
                        MoreArgs = NULL)

  # Reconstruct the merged partition by splitting partitions based
  # on the components
  partition <- unsplit(partitions, components)

  # Change the IDs in the merged partition
  standardize_partition(partition)
}


#' @describeIn check_solution Verifies if a partition check all connectivity
#' and size constraints
#' @param m,M minimum and maximum size constraints. (positive values)
#' @param sizes sizes of the elements (positive values)
#' @export
#' @importFrom checkmate assertNumber
is_feasible_solution <- function(partition, contiguity = NULL,
                                 sizes = rep(1L, length(partition)),
                                 m = 0.0, M = Inf)
{
  if (inherits(partition, "Partition"))
    return(all(partition$contraintes))

  assertPartition(partition)
  if (!is.null(contiguity))
  {
    isRegion <- is_regionalisation(partition, contiguity)
    if (!isRegion)
      return(FALSE)
  }

  assertNumber(m, lower = 0.0, finite = TRUE)
  assertNumber(M, lower = m, finite = FALSE)
  clustersSizes <- .clusters_sizes(partition, sizes)

  all(m <= clustersSizes & clustersSizes <= M)
}

setGeneric("is_feasible_solution")
