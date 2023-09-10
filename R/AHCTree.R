#' @include partition_class.R
#' @include pbCon.R


# `AHCTree` class definition ----------------------------------------------


#' An `AHCTree` class representing a hierarchical clustering tree.
#'
#' This class contains fields representing a hierarchical clustering tree and
#' its partitions. It provides various methods for manipulating and
#' analyzing hierarchical clustering results.
#' @field pb An object of class `pbCon` representing the clustering problem.
#' @field partitions A list of Partition objects representing
#' the hierarchical partitions.
#' @keywords internal
AHCTree <- setRefClass(
  "AHCTree",
  fields = list(pb = "pbCon",
                partitions = "list")
)


#' Create an `AHCTree` object.
#'
#' This function creates an `AHCTree` object with the given clustering problem
#' and partitions. The partitions can be provided as a list of Partition
#' objects or as a single Partition object with the number of desired partitions
#' (n). If n is specified and greater than 1, it will create
#' n copies of the given Partition object.
#'
#' @param pb An object of class `pbCon` representing the clustering problem.
#' @param partitions A list of Partition objects or a single Partition object.
#' @param n Number of partitions to create if a single Partition object
#' is provided. Default is NULL.
#' @returns An `AHCTree` object with the specified clustering problem
#' and partitions.
#' @keywords internal
#' @importFrom checkmate assertClass assertList assertInt
constructor_AHCTree <- function(pb, partitions = list(), n = NULL)
{
  assertClass(pb, "pbCon")
  if (is.list(partitions))
    assertList(partitions, types = "Partition")

  else
  {
    assertClass(partitions, "Partition")
    nbClustersPartition <- partitions$k()

    if (is.null(n))
      n <- nbClustersPartition
    else
      assertInt(n, lower = 1L, upper = nbClustersPartition)


    if (n > 1L)
      partitions <- replicate(n, partitions$copy(shallow = TRUE),
                              simplify = FALSE)

    names(partitions) <- nbClustersPartition - 0L:(n - 1L)
  }

  AHCTree$new(pb = pb, partitions = partitions)
}

#' Properties of a `AHCTree` object
#' @param x A `AHCTree` object
#' @name AHCTree_properties
#' @rdname AHCTree_properties
#' @keywords internal
NULL

#' @describeIn AHCTree_properties Give the number
#' of partitions stored in the tree. (positive integer)
#' @keywords internal
setMethod(
  "length",
  signature(x = "AHCTree"),
  function(x) length(x$partitions)
)

#' @describeIn AHCTree_properties Give the names
#' of the different partitions stored in the tree. If
#' the tree have not been manually modified then it corresponds
#' to the number of clusters in each partition.
#' @keywords internal
setMethod("names", signature(x = "AHCTree"), function(x) names(x$partitions))

#' Does partitions of a `AHCTree` verifies constraints?
#' @param x A `AHCTree`
#' @name AHCTree_check_constraints
#' @rdname AHCTree_check_constraints
#' @keywords internal
NULL


# Access to AHCTree data --------------------------------------------------


#' Access to the partitions of a `AHCTree`
#' @name AHCTree_access
#' @rdname AHCTree_access
#' @param x A `AHCTree` object
#' @param i As logical, strictly positive integers or character
#' vector.
#' @returns for "[", a list of `Partition` objects ; for "[[" a `Partition`
#' object if `i` is of length 1, a list of `Partition` objects otherwise.
#' @keywords internal
NULL

#' @rdname AHCTree_access
#' @keywords internal
setMethod(
  "[",
  signature(x = "AHCTree", i = "logical"),
  function(x, i)
  {
    x$partitions <- x$partitions[i]
    return(x)
  }
)

#' @rdname AHCTree_access
#' @keywords internal
setMethod(
  "[", # nolint: indentation_linter
  signature(x = "AHCTree", i = "vector"),
  function(x, i)
  {
    x$partitions <- x$partitions[as.character(i)]
    return(x)
  }
)

#' @rdname AHCTree_access
#' @keywords internal
setMethod(
  "[[",
  signature(x = "AHCTree", i = "numeric"),
  function(x, i)
  {
    x[[as.character(i)]]
  }
)

#' @rdname AHCTree_access
#' @keywords internal
setMethod(
  "[[",
  signature(x = "AHCTree", i = "character"),
  function(x, i)
  {
    if (length(i == 1L))
      x$partitions[[i]]
    else
      x$partitions[i]
  }
)

#' @rdname AHCTree_access
#' @keywords internal
setMethod(
  "[[",
  signature(x = "AHCTree", i = "missing"),
  function(x, i) x$partitions
)

loadNamespace("stats")
setGeneric("cutree", stats::cutree)

setMethod("cutree",
          signature(tree = "AHCTree", k = "numeric", h = "ANY"),
          function(tree, k, h) tree[[k]]$partition)

# Definition of equality with other types of objects ----------------------

setMethod(
  "all.equal",
  signature(target = "AHCTree", current = "AHCTree"),
  function(target, current)
  {

    targetPartitions <- lapply(target$partitions, function(p) p$partition)
    currentPartitions <- lapply(current$partitions, function(p) p$partition)

    nbClustersTarget <- vapply(targetPartitions, nbClusters, integer(1L))
    nbClustersCurrent <- vapply(currentPartitions, nbClusters, integer(1L))

    if (length(targetPartitions) != length(currentPartitions))
      return(FALSE)

    names(targetPartitions) <- nbClustersTarget
    names(currentPartitions) <- nbClustersCurrent

    all.equal(targetPartitions, currentPartitions)
  }
)


# Conversion --------------------------------------------------------------


#' Convert an `AHCTree` object to a `tibble`.
#'
#' This function converts an `AHCTree` object to a
#' [tibble][tibble::tibble-package] format,
#' suitable for further analysis and visualization.
#'
#' @param x An `AHCTree` object.
#' @returns A `tibble` representing the `AHCTree` object's partitions.
#' @noRd
#' @importFrom tibble as_tibble
setGeneric("as_tibble", tibble::as_tibble)

#' @importFrom tibble as_tibble
setMethod("as_tibble", signature(x = "AHCTree"),
          function(x) Partition_list_to_tibble(x$partitions))
