#' Give a score of constraint size respect
#'
#' Those functions give a score that indicates if a partition respects
#' or not some size constraint(s). A score equals to 0 means that the constraint
#' is completely respected. A strictly positive value means that some clusters
#' do not respect the constraint size.
#' @param pb A `pbCon` object.
#' @param clustersSizes the size of the different clusters. If equals to `NULL`
#' (default) sizes will be calculated.
#' @param partition a partition for the problem. Used only is `clustersSizes`
#' is `NULL`.
#' @param f a function that will take a vector of difference between cluster
#' sizes and constraint value. Must be vectorial, give 0 when values are
#' negative and positive otherwise.
#' @keywords internal
#' @name score_size_constraints
#' @rdname score_size_constraints
NULL

NEEDCLUSTERSIZESORPART <- "clusters sizes or partition must be given."

#' @rdname score_size_constraints
#' @keywords internal
score_constraints_min <- function(pb, clustersSizes = NULL,
                                  partition = NULL,
                                  f = identity)
{
  if (is.null(clustersSizes) && is.null(partition))
  {
    stop(NEEDCLUSTERSIZESORPART)
  }

  if (is.null(clustersSizes))
    clustersSizes <- clusters_sizes(partition, pb$sizes)

  .score_constraints_min(pb, clustersSizes, f)
}

#' @rdname score_size_constraints
#' @keywords internal
#' @importFrom checkmate assertNumber
.score_constraints_min <- function(pb, clustersSizes, f = identity)
{
  if (!pb$hasMinConstraint())
    return(0.)

  petitsClusters <- clustersSizes < pb$m

  if (any(petitsClusters))
  {
    score <- sum(f(pb$m - clustersSizes[petitsClusters]))
    assertNumber(score, lower = 0.0)
    return(score)
  }

  else
    return(0.)
}

#' @rdname score_size_constraints
#' @keywords internal
score_constraints_max <- function(pb, clustersSizes = NULL,
                                  partition = NULL,
                                  f = identity)
{
  if (is.null(clustersSizes) && is.null(partition))
  {
    stop(NEEDCLUSTERSIZESORPART)
  }

  if (is.null(clustersSizes))
    clustersSizes <- clusters_sizes(partition, pb$sizes)

  .score_constraints_max(pb, clustersSizes, f)
}

#' @rdname score_size_constraints
#' @keywords internal
#' @importFrom checkmate assertNumber
.score_constraints_max <- function(pb, clustersSizes, f = identity)
{
  if (!pb$hasMaxConstraint())
    return(0.0)

  grosClusters <- clustersSizes > pb$M

  if (any(grosClusters))
  {
    score <- sum(f(clustersSizes[grosClusters] - pb$M))
    assertNumber(score, lower = 0.0)
    return(score)
  }

  return(0.0)
}

#' @rdname score_size_constraints
#' @keywords internal
score_constraints_table <- function(scoreContraintesMin,
                                    scoreContraintesMax)
{
  total <- scoreContraintesMin + scoreContraintesMax
  scoreContraintes <- c(scoreContraintesMin,
                        scoreContraintesMax,
                        total)

  names(scoreContraintes) <- c("min", "max", "total")
  return(scoreContraintes)
}

#' @rdname score_size_constraints
#' @param details flag indicating if details must be given or only
#' the final score.
#' @keywords internal
#' @importFrom checkmate assertFlag
score_constraints <- function(pb, clustersSizes = NULL,
                              partition = NULL,
                              f = identity, details = FALSE)
{
  assertFlag(details)

  if (is.null(clustersSizes) && is.null(partition))
  {
    stop(NEEDCLUSTERSIZESORPART)
  }

  if (is.null(clustersSizes))
    clustersSizes <- clusters_sizes(partition, pb$sizes)

  scoreContraintesMin <- .score_constraints_min(pb, clustersSizes, f)
  scoreContraintesMax <- .score_constraints_max(pb, clustersSizes, f)


  if (details)
    return(score_constraints_table(scoreContraintesMin,
                                   scoreContraintesMax))

  else
    return(scoreContraintesMin + scoreContraintesMax)
}