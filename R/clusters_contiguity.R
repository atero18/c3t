#' @include pbCon.R

#' @title List Contiguous Classes
#'
#' @description Function returning the set of contiguous classes in a problem.
#'
#' @param partition Partition vector. Must be numeric, containing integers,
#' and non-empty. (numeric vector)
#' @param x A contiguity matrix or a contiguity graph.
#' The number of elements in the matrix/graph should correspond to the size
#' of partition.
#' (matrix / contiguity graph)
#'
#' @returns A 2-column matrix containing a list of contiguous classes.
#' By definition, a class is contiguous with itself, so this data is not stored
#' for memory efficiency. Similarly, only pairs (i,j) where i and j
#' are contiguous and i < j are stored.
#'
#' @importFrom igraph is_igraph
#' @importFrom methods validObject
#' @export
clusters_contiguity_list <- function(partition, x)
{
  if (is_pbCon(x))
  {
    x <- getmatContPbCon(x)
  }
  else if (is.matrix(x))
  {
    if (inherits(x, "ContiguityMat"))
      validObject(x)

    else
      assertContiguityMatrix(x, isComplete = FALSE)
  }
  else if (is_igraph(x))
    x <- graph_to_contiguity_matrix(x)

  if (inherits(x, "ContiguityMat"))
  {
    validObject(x)
    x <- x[]
  }

  else if (inherits(x, "Matrix"))
    x <- as.matrix(x)

  else if (!is.matrix(x))
  {
    stop("`x` must be a matrix or graph contiguity")
  }

  assertPartition(partition, nrow(x))

  n <- length(partition)


  clustersIDs <- clustersIDs(partition)

  # Cse of complete contiguity
  diag(x) <- TRUE
  if (all(x))
  {
    contiguousClusters <- as.matrix(expand.grid(clustersIDs,
                                                clustersIDs))

    contiguousClusters <- contiguousClusters[contiguousClusters[, 1L] <
                                               contiguousClusters[, 2L], ]
    return(contiguousClusters)
  }

  clustersList <- partition_to_list(partition)

  k <- length(clustersList)

  # If there are une cluster per element
  if (k == n)
  {
    contiguousClusters <- which(x[], arr.ind = TRUE)
    masque <- contiguousClusters[, 1L] < contiguousClusters[, 2L]
    contiguousClusters <- contiguousClusters[masque, , drop = FALSE]
    if (nrow(contiguousClusters) > 0L)
    {
      contiguousClusters[, 1L] <- partition[contiguousClusters[, 1L]]
      contiguousClusters[, 2L] <- partition[contiguousClusters[, 2L]]
    }
    return(contiguousClusters)
  }

  if (0L %in% clustersIDs)
  {
    stop("0 cannot be a cluster ID")
  }

  x <- x * matrix(partition, ncol = n, nrow = n, byrow = TRUE)

  res <- c3t_lapply(clustersIDs, function(c)
  {
    contiguousClusters <- setdiff(x[clustersList[[as.character(c)]], ],
                                  c(0L, c))
    contiguousClusters <- contiguousClusters[contiguousClusters > c]
    mat <- matrix(NA, nrow = length(contiguousClusters), ncol = 2L)
    mat[, 1L] <- c
    mat[, 2L] <- contiguousClusters
    return(mat)
  })

  res <- do.call("rbind", res)
  colnames(res) <- c("cluster1", "cluster2")

  return(res)
}


#' @describeIn clusters_contiguity_list Returns the contiguity matrix
#' instead of a list.
#' @keywords internal
#' @importFrom Matrix sparseMatrix
clusters_contiguity_matrix <- function(partition, elemsContiguity,
                                       diag = FALSE)
{
  contiguityList <- clusters_contiguity_list(partition, elemsContiguity)

  contiguite_classes_mat <-
    sparseMatrix(i = c(contiguityList[, 1L], contiguityList[, 2L]),
                 j = c(contiguityList[, 2L], contiguityList[, 1L]),
                 x = TRUE)

  diag(contiguite_classes_mat) <- diag

  contiguite_classes_mat
}
