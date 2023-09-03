#' Validate Size Constraints
#' @param sizes A numeric vector representing the sizes
#' of elements.
#' @param composantesConnexes_vec An optional numeric vector specifying
#' the connected components.
#' @name valid_cont
#' @noRd
NULL


#' @describeIn valid_cont Verifies the minimum size constraint for
#' a given vector and its connected components.
#'
#' @param m_vec A numeric vector with min constraints values to be tested.
#' @keywords internal
#' @noRd
#' @importFrom checkmate testNumeric
contrainteMinValide <-
  function(m_vec, sizes,
           composantesConnexes_vec = rep(1L, length(sizes)))
{
  tailleComposantesConnexes <-
    tapply(sizes, composantesConnexes_vec, sum, simplify = TRUE)

  testNumeric(m_vec, upper = min(tailleComposantesConnexes),
              any.missing = FALSE, all.missing = FALSE)
}

#' @describeIn valid_cont Verifies the maximum size constraint
#' for a given vector.
#' @param M_vec A numeric vector with max constraints to be tested.
#' @keywords internal
#' @noRd
#' @importFrom checkmate testNumeric
contrainteMaxValide <- function(M_vec, sizes)
{
  testNumeric(M_vec, lower = max(sizes),
              any.missing = FALSE, all.missing = FALSE)
}

#' Simplify Minimum Size Constraint
#'
#' Simplifies the minimum size constraint of a given vector.
#'
#' @param m_vec A numeric vector with min constraints to be simplified.
#' @inheritParams valid_cont
#' @name simp_cont
#' @keywords internal
#' @noRd
simplificationContrainteMin <- function(m_vec, sizes)
{
  masqueContrainteInutile <- m_vec <= min(sizes)
  if (any(masqueContrainteInutile))
    m_vec[masqueContrainteInutile] <- 0.0
  m_vec
}

#' Simplify Maximum Size Constraint
#'
#' @describeIn simp_cont Simplifies the maximum size constraint
#' of a given vector.
#'
#' @param M_vec A numeric vector with max constraints values to be simplified.
#' @inheritParams valid_cont
#' @keywords internal
#' @noRd
simplificationContrainteMax <-
  function(M_vec, sizes,
           composantesConnexes_vec = rep(1L, length(sizes)))
{
  tailleComposantesConnexes <-
    tapply(sizes, composantesConnexes_vec, sum, simplify = TRUE)
  masqueContrainteInutile <- M_vec >= max(tailleComposantesConnexes)
  if (any(masqueContrainteInutile))
    M_vec[masqueContrainteInutile] <- Inf
  M_vec
}

#' Cluster Sizes
#'
#' @description Determines the sizes of clusters for a given
#' partition and elements sizes.
#'
#' @param partition A numeric vector representing the partition of
#' the elements.
#' @param sizes A numeric vector representing the sizes of
#' individual elements.
#' @name clusters_sizes
#' @export
#' @importFrom checkmate assertNumeric
clusters_sizes <-
  function(partition, sizes = rep(1.0, length(partition)))
{
  # Argument verification
  assertPartition(partition)
  assertNumeric(sizes, lower = 0.0,
                finite = TRUE,
                len = length(partition), any.missing = FALSE)

  tapply(sizes, partition, sum)
}

#' @importFrom methods setGeneric
setGeneric("clusters_sizes", clusters_sizes)

#' Clusters Too Small
#'
#' @describeIn clusters_sizes Determines the clusters that are too small
#' based on a given size constraint.
#'
#' @param m_num A numeric value representing the minimum
#' cluster size constraint.
#' @keywords internal
#' @noRd
clusters_trop_petits <- function(m_num, partition, sizes)
{
  taillesClusters <- clusters_sizes(partition, sizes)
  tooSmallClusters <- names(taillesClusters)[taillesClusters < m_num]
  tryCatch(expr = as.numeric(tooSmallClusters),
           warning = function(cond) tooSmallClusters)
}


#' @importFrom methods setGeneric
setGeneric("clusters_trop_petits", clusters_trop_petits)

#' Number of Clusters Too Small
#'
#' @describeIn clusters_sizes Determines the number of clusters that
#' are too small based on a given size constraint.
#' @keywords internal
#' @noRd
nb_clusters_trop_petits <- function(m_num, partition, sizes)
{
  as.integer(sum(clusters_sizes(partition, sizes) < m_num))
}

#' Clusters Too Large
#'
#' @describeIn clusters_sizes Determines the clusters that are too large
#' based on a given size constraint.
#'
#' @param M_num A numeric value representing the maximum cluster
#' size constraint.
#' @keywords internal
#' @noRd
clusters_trop_gros <- function(M_num, partition, sizes)
{
  taillesClusters <- clusters_sizes(partition, sizes)
  tooBigClusters <- names(taillesClusters)[taillesClusters > M_num]
  tryCatch(expr = as.integer(tooBigClusters),
           warning = function(cond) tooBigClusters)
}


#' @importFrom methods setGeneric
setGeneric("clusters_trop_gros", clusters_trop_gros)

#' Number of Clusters Too Large
#'
#' @describeIn clusters_sizes Determines the number of clusters that are
#' too large based on a given size constraint.
#' @keywords internal
#' @noRd
nb_clusters_trop_gros <- function(M_num, partition, sizes)
{
  as.integer(sum(clusters_sizes(partition, sizes) > M_num))
}
