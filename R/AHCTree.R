#' @include partition_class.R
#' @include pbCon.R

#' An `AHCTree` class representing a hierarchical clustering tree.
#'
#' This class contains fields representing a hierarchical clustering tree and
#' its partitions. It provides various methods for manipulating and
#' analyzing hierarchical clustering results.
#' @field pb An object of class `pbCon` representing the clustering problem.
#' @field partitions A list of Partition objects representing
#' the hierarchical partitions.
#' @keywords internal
AHCTree <- setRefClass("AHCTree",
                       fields = list(pb = "pbCon",
                                     partitions = "list"))

# ... (Rest of the class methods with appropriate roxygen comments)

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
AHCTree_constructor <- function(pb, partitions = list(), n = NULL)
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
setMethod("length", signature(x = "AHCTree"),
          function(x) length(x$partitions))

#' @describeIn AHCTree_properties Give the names
#' of the different partitions stored in the tree. If
#' the tree have not been manually modified then it corresponds
#' to the number of clusters in each partition.
#' @keywords internal
setMethod("names", signature(x = "AHCTree"),
          function(x) names(x$partitions))

setMethod(".nbClusters", signature(partition = "AHCTree"),
          function(partition)
          {
            if (length(partition) == 0L)
              return(NULL)

            else if (!is.null(names(partition)))
              return(as.numeric(names(partition)))

            else
              return(vapply(partition$partitions, .nbClusters, integer(1L)))
          }
)

setMethod("nbSingletons", signature(partition = "AHCTree"),
          function(partition)
          {
            if (length(partition) == 0L)
              return(NULL)

            else
              return(vapply(partition$partition, nbSingletons, integer(1L)))
          }
)

setMethod("checkContiguityConst", signature(x = "AHCTree"),
          function(x)
          {
            if (length(x) == 0L)
              return(NULL)

            else
              return(vapply(x$partitions,
                            checkContiguityConst,
                            logical(1L)))
          })

setMethod("checkMinSizeConst", signature(x = "AHCTree"),
          function(x)
          {
            if (length(x) == 0L)
              return(NULL)

            else
              return(vapply(x$partitions,
                            checkMinSizeConst,
                            logical(1L)))
          })

setMethod("checkMaxSizeConst", signature(x = "AHCTree"),
          function(x)
          {
            if (length(x) == 0L)
              return(NULL)

            else
              return(vapply(x$partitions, checkMaxSizeConst, logical(1L)))
          })

#' Does partitions of a `AHCTree` verifies constraints?
#' @param x A `AHCTree`
#' @name AHCTree_check_constraints
#' @rdname AHCTree_check_constraints
#' @keywords internal
NULL

setMethod("scoreMinSizeConst", signature(x = "AHCTree"),
          function(x)
          {
            if (length(x) == 0L)
              return(NULL)
            else
              return(vapply(x$partitions, scoreMinSizeConst,
                            numeric(1L)))
          })

setMethod("scoreMaxSizeConst", signature(x = "AHCTree"),
          function(x)
          {
            if (length(x) == 0L)
              return(NULL)

            else
              return(vapply(x$partitions, scoreMinSizeConst,
                            numeric(1L)))
          })

setMethod("scoreSizeConsts", signature(x = "AHCTree"),
          function(x)
          {
            if (length(x) == 0L)
              return(NULL)

            else
              return(vapply(x$partitions, scoreSizeConsts, numeric(1L)))
          })

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
setMethod("[", signature(x = "AHCTree", i = "logical"),
          function(x, i)
          {
            x$partitions <- x$partitions[i]
            return(x)
          }
)

#' @rdname AHCTree_access
#' @keywords internal
setMethod("[", signature(x = "AHCTree", i = "vector"),
          function(x, i)
          {
            x$partitions <- x$partitions[as.character(i)]
            return(x)
          })

#' @rdname AHCTree_access
#' @keywords internal
setMethod("[[", signature(x = "AHCTree", i = "numeric"),
          function(x, i)
          {
            x[[as.character(i)]]
          }
)

#' @rdname AHCTree_access
#' @keywords internal
setMethod("[[", signature(x = "AHCTree", i = "character"),
          function(x, i)
          {
            if (length(i) == 1L)
              return(x$partitions[[i]])

            else
              return(x$partitions[i])
          }
)

#' @rdname AHCTree_access
#' @keywords internal
setMethod("[[", signature(x = "AHCTree", i = "missing"),
          function(x, i)
          {
            return(x$partitions)
          }
)

#' Clusters sizes of a AHC tree partitions.
#' Return in a list the sizes of the clusters of each partition
#' of an `AHCTree`.
#' @inheritParams clusters_sizes
#' @param partition An `AHCTree`.
#' @returns a list with for each partition of the `AHCTree`
#' a vector containing the size of each cluster.
#' @keywords internal
setMethod("clusters_sizes", signature(partition = "AHCTree",
                                      sizes = "ANY"),
          function(partition, sizes)
          {
            if (length(partition) == 0L)
              return(NULL)

            else
              return(lapply(partition$partitions,
                            clusters_sizes))
          })


#' @importFrom methods setGeneric
setGeneric("cutree", stats::cutree)

setMethod("cutree", signature(tree = "AHCTree", k = "numeric", h = "ANY"),
          function(tree, k, h) tree[[k]]$partition)

setMethod("all.equal", signature(target = "AHCTree", current = "AHCTree"),
          function(target, current)
          {
            nbClustersTarget <- nbClusters(target)
            nbClustersCurrent <- nbClusters(current)

            if (!setequal(nbClustersCurrent, nbClustersTarget))
               return(FALSE)

            compPartitions <-
              vapply(seq_along(nbClustersCurrent),
                     function(i)
                     {
                       all.equal(target$partitions[[i]],
                                current$partitions[[i]]) # nolint: indentation_linter
                     },
                     logical(1L))
            all(compPartitions)
          }
)

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
#' @importFrom methods setGeneric
setGeneric("as_tibble", as_tibble)

#' @importFrom tibble as_tibble
setMethod("as_tibble", signature(x = "AHCTree"),
          function(x) Partition_list_to_tibble(x$partitions))
