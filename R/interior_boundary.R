#' @include pbCon.R

# Determines the set of elements from a cluster that are contiguous
# with at least one element from another cluster.
# @param cluster_int The identifier of the cluster (strictly positive integer)
interior_boundary_cluster <- function(contiguity,
                                      cluster_int,
                                      partition,
                                      removeSymmetry = TRUE)
{
  COLNAMES <- c("x1", "x2", "cluster1", "cluster2")
  EMPTYCONTIGUITIES <- matrix(nrow = 0L, ncol = 4L,
                              dimnames = list(NULL, COLNAMES))

  # Identify which elements belong to our cluster
  masque <-  partition == cluster_int

  # If the cluster is empty / does not exist, there is no boundary
  if (!any(masque))
  {
    warning("There is no element which are in `cluster_int`")
    return(EMPTYCONTIGUITIES)
  }

  if (is_pbCon(contiguity))
    contiguity <- getmatContPbCon(contiguity)

  diag(contiguity) <- FALSE

  # Case where there's only one element in the cluster
  if (sum(masque) == 1L)
  {
    elementCluster <- which(masque)
    contiguites <- which(contiguity[elementCluster, ])

    if (removeSymmetry)
      contiguites <- contiguites[contiguites < elementCluster]

    contiguites <- cbind(elementCluster, contiguites)
  }

  else
  {

    if (removeSymmetry)
      contiguity[lower.tri(contiguity)] <- FALSE

    # Mask the elements:
    # - Those for which j > i (idElement) since i < j will be handled
    #   with another cluster
    # - Rows not belonging to the cluster
    # - Columns belonging to the cluster
    contiguity[!masque, ] <- FALSE
    contiguity[, masque] <- FALSE

    contiguites <- which(contiguity, arr.ind = TRUE)
  }

  if (nrow(contiguites) > 0L)
  {
    contiguites <- cbind(contiguites, cluster_int,
                         partition[contiguites[, 2L]])

    rownames(contiguites) <- NULL
    colnames(contiguites) <- COLNAMES
  }
  else
    contiguites <- EMPTYCONTIGUITIES



  contiguites
}


#' Interior boundary of a partition
#'
#' Determines the set of elements, distributed in clusters, that are
#' contiguous to an element from another cluster.
#'
#' @param contiguity A contiguity matrix or a contiguity graph
#' (package `igraph`). Prefer matrix is possible.
#' (contiguity matrix or `igraph` graph)
#' @param partition The partition. (vector of strictly positive integers)
#' @param removeSymmetry `TRUE` (default) if symmetry should be removed
#' from the result. (flag)
#' @returns a matrix with 4 columns (`x1`, `x2`, `cluster1`, `cluster2`) with
#' as many rows as possible connections between an element of a cluster and
#' another. If `removeSymmetry = TRUE`, only rows where `x1 < x2` are kept.
#' @details The interior boundary of a partition, under a certain
#' contiguity, is defined as the set of elements contiguous to at least
#' one other element from a different cluster than its own.
#' @importFrom igraph ends E is_igraph
#' @export
interior_boundary <- function(contiguity,
                              partition,
                              removeSymmetry = TRUE)
{
  if (is_pbCon(contiguity))
    contiguity <- getmatContPbCon(contiguity)

  else if (is_igraph(contiguity))
    contiguity <- graph_to_contiguity_matrix(contiguity)

  diag(contiguity) <- FALSE

  if (removeSymmetry)
    contiguity[lower.tri(contiguity)] <- FALSE

  # Retrieve contiguities as a list, in a 2-column matrix
  elementsBoundary_mat <- which(contiguity, arr.ind = TRUE)

  if (!is.matrix(elementsBoundary_mat))
    elementsBoundary_mat <- matrix(elementsBoundary_mat, ncol = 2L,
                                   byrow = TRUE)

  if (nrow(elementsBoundary_mat) > 0L)
  {
    # Add 2 columns (cluster1 & cluster2) indicating the cluster of each
    # element in the different connections
    elementsBoundary_mat <- cbind(elementsBoundary_mat,
                                  partition[elementsBoundary_mat[, 1L]],
                                  partition[elementsBoundary_mat[, 2L]])

    # Remove rows where the connection involves 2 elements from
    # the same cluster
    masque <- elementsBoundary_mat[, 3L] != elementsBoundary_mat[, 4L]
    elementsBoundary_mat <- elementsBoundary_mat[masque, , drop = FALSE]
  }


  rownames(elementsBoundary_mat) <- NULL
  colnames(elementsBoundary_mat) <- c("x1", "x2", "cluster1", "cluster2")

  return(elementsBoundary_mat)
}
