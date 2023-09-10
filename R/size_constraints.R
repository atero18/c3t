# Checking size constraints arguments -------------------------------------

#' Validate Size Constraints
#' @param sizes A numeric vector representing the sizes
#' of elements.
#' @param connectedComponents An optional numeric vector specifying
#' the connected components.
#' @name valid_cont
NULL


#' @describeIn valid_cont Verifies the minimum size constraint for
#' a given vector and its connected components.
#'
#' @param m_vec A numeric vector with min constraints values to be tested.
#' @keywords internal
#' @importFrom checkmate assertFlag assertVector checkNumeric assertCount
checkMinSizeConstraints <- function(m_vec, sizes, connectedComponents = NULL,
                                    authorizeNegative = FALSE,
                                    len = NULL)
{

  assertCount(len, positive = TRUE, null.ok = TRUE)

  if (!is.null(connectedComponents))
  {
    assertVector(connectedComponents, min.len = 1L,
                 any.missing = FALSE, all.missing = FALSE,
                 null.ok = FALSE)

    assertSizes(sizes, len = length(connectedComponents))

    nbComponents <- length(unique(connectedComponents))

  }
  else
  {
    assertSizes(sizes)
    nbComponents <- 1L
  }


  if (nbComponents == 1L)
    upper <- sum(sizes)

  else
  {

    sizesConnectedComp <-
      connectedComponents_sizes(connectedComponents, sizes)

    upper <- min(sizesConnectedComp)
  }

  assertFlag(authorizeNegative)
  lower <- ifelse(authorizeNegative, -Inf, 0.0)

  checkNumeric(m_vec, min.len = 1L,
               lower = lower, upper = upper, finite = FALSE,
               null.ok = FALSE, len = len)
}

#' @importFrom checkmate makeAssertionFunction
assertMinSizeConstraints <- makeAssertionFunction(checkMinSizeConstraints)

#' @importFrom checkmate makeTestFunction
testMinSizeConstraints <- makeTestFunction(checkMinSizeConstraints)

#' @describeIn valid_cont Verifies the maximum size constraint
#' for a given vector.
#' @param M_vec A numeric vector with max constraints to be tested.
#' @keywords internal
#' @importFrom checkmate assertFlag assertVector checkNumeric assertCount
checkMaxSizeConstraints <- function(M_vec, sizes, len = NULL)
{
  assertCount(len, positive = TRUE, null.ok = TRUE)
  assertSizes(sizes)
  checkType <- checkDouble(M_vec,
                           min.len = 1L, len = len,
                           lower = 0.0, finite = FALSE,
                           any.missing = FALSE, all.missing = FALSE,
                           null.ok = FALSE)

  if (!isTRUE(checkType))
    return(checkType)



  tooLowValues <- vapply(M_vec, function(M) any(sizes > M), logical(1L))
  if (any(tooLowValues))
  {
    return(paste("following maximum constraints values are too low:",
                 unique(M_vec[tooLowValues])))
  }

  TRUE

}

#' @importFrom checkmate makeAssertionFunction
assertMaxSizeConstraints <- makeAssertionFunction(checkMaxSizeConstraints)

#' @importFrom checkmate makeTestFunction
testMaxSizeConstraints <- makeTestFunction(checkMaxSizeConstraints)



# Simplify size constraint values -----------------------------------------


#' Simplify Minimum Size Constraint
#'
#' Simplifies the minimum size constraint of a given vector.
#'
#' @param m_vec A numeric vector with min constraints to be simplified.
#' @inheritParams valid_cont
#' @name simp_cont
#' @keywords internal
simplify_minSizeConst <- function(m_vec, sizes = NULL)
{
  if (is.null(m_vec))
    return(0.0)

  if (anyNA(m_vec))
    m_vec[is.na(m_vec)] <- 0.0

  if (!is.null(sizes))
  {
    masqueContrainteInutile <- m_vec <= min(sizes)
    if (any(masqueContrainteInutile))
      m_vec[masqueContrainteInutile] <- 0.0
  }


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
simplify_maxSizeConst <- function(M_vec, sizes, connectedComponents = NULL)
{
  if (is.null(M_vec))
    return(Inf)

  if (is.null(connectedComponents))
    upper <- sum(sizes)

  else
  {
    sizeConnectedComp <- connectedComponents_sizes(connectedComponents, sizes)
    upper <- max(sizeConnectedComp)
  }

  maskUselessConstraint <- M_vec >= upper
  if (any(maskUselessConstraint))
    M_vec[maskUselessConstraint] <- Inf

  M_vec
}


# Clusters sizes ----------------------------------------------------------

.clusters_sizes <- function(partition, sizes)
{
  tapply(sizes, partition, sum, simplify = TRUE)
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
clusters_sizes <-
  function(partition, sizes = NULL)
{

  # Argument verification
  assertPartition(partition)
  n <- length(partition)
  if (is.null(sizes))
    sizes <- one_person_per_case(n, "integer")
  else
    assertSizes(sizes, len = n)

  .clusters_sizes(partition, sizes)
}

setGeneric("clusters_sizes", clusters_sizes)


#' @importFrom checkmate assertVector
connectedComponents_sizes <- function(components, sizes = NULL)
{
  clusters_sizes(components, sizes)
}


#' Clusters Too Small
#'
#' @describeIn clusters_sizes Determines the clusters that are too small
#' based on a given size constraint.
#'
#' @param m_num A numeric value representing the minimum
#' cluster size constraint.
#' @keywords internal
#' @noRd
too_small_clusters <- function(m_num, partition, sizes)
{
  taillesClusters <- clusters_sizes(partition, sizes)
  tooSmallClusters <- names(taillesClusters)[taillesClusters < m_num]
  tryCatch(expr = as.numeric(tooSmallClusters),
           warning = function(cond) tooSmallClusters)
}


setGeneric("too_small_clusters", too_small_clusters)

#' Number of Clusters Too Small
#'
#' @describeIn clusters_sizes Determines the number of clusters that
#' are too small based on a given size constraint.
#' @keywords internal
#' @noRd
nb_too_small_clusters <- function(m_num, partition, sizes)
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
too_big_clusters <- function(M_num, partition, sizes)
{
  taillesClusters <- clusters_sizes(partition, sizes)
  tooBigClusters <- names(taillesClusters)[taillesClusters > M_num]
  tryCatch(expr = as.integer(tooBigClusters),
           warning = function(cond) tooBigClusters)
}


setGeneric("too_big_clusters", too_big_clusters)

#' Number of Clusters Too Large
#'
#' @describeIn clusters_sizes Determines the number of clusters that are
#' too large based on a given size constraint.
#' @keywords internal
#' @noRd
nb_too_big_clusters <- function(M_num, partition, sizes)
{
  as.integer(sum(clusters_sizes(partition, sizes) > M_num))
}
