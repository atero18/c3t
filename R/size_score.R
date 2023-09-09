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
#' @returns a positive double value indicating how much the partition
#' respect the size constraint(s). If it is equals to 0 then all constraint
#' are respected by each cluster. The higher the score is the less the
#' constraint(s) are respected. (positive double value)
#' @keywords internal
#' @name score_size_constraints
#' @rdname score_size_constraints
NULL

NEEDCLUSTERSIZESORPART <- "clusters sizes or partition must be given."

#' @rdname score_size_constraints
#' @keywords internal
#' @importFrom checkmate assertNumber assertFunction
score_constraints_min <- function(clustersSizes =
                                    clusters_sizes(partition, sizes),
                                  m = 0.0,
                                  sizes = NULL,
                                  partition = NULL,
                                  f = identity)
{
  assertNumber(m, lower = 0.0, finite = TRUE, null.ok = FALSE)
  assertFunction(f)

  if (!"clustersSizes" %in% names(as.list(match.call())))
  {
    if (is.null(sizes) || is.null(partition))
    {
      stop(NEEDCLUSTERSIZESORPART)
    }
    assertPartition(partition)
    n <- length(partition)
    assertSizes(sizes, len = n)
    clustersSizes <- clusters_sizes(partition, sizes)
  }


  .score_constraints_min(m, clustersSizes, f)
}

#' @rdname score_size_constraints
#' @keywords internal
#' @importFrom checkmate assertNumber
.score_constraints_min <- function(m, clustersSizes, f = identity)
{
  if (m == 0.0)
    return(0.0)

  petitsClusters <- clustersSizes < m

  if (any(petitsClusters))
  {
    score <- sum(f(m - clustersSizes[petitsClusters]))
    assertNumber(score, lower = 0.0)
    return(as.double(score))
  }

  else
    return(0.0)
}

#' @rdname score_size_constraints
#' @keywords internal
#' @importFrom checkmate assertNumber assertFunction
score_constraints_max <- function(clustersSizes =
                                    clusters_sizes(partition, sizes),
                                  M = Inf,
                                  partition = NULL,
                                  f = identity)
{
  assertNumber(M, lower = 0.0, finite = FALSE, null.ok = FALSE)
  assertFunction(f)
  if (!"clustersSizes" %in% names(as.list(match.call())))
  {
    if (is.null(sizes) || is.null(partition))
    {
      stop(NEEDCLUSTERSIZESORPART)
    }
    assertPartition(partition)
    n <- length(partition)
    assertSizes(sizes, len = n, M = M)
    clustersSizes <- clusters_sizes(partition, sizes)
  }
  .score_constraints_max(M, clustersSizes, f)
}

#' @rdname score_size_constraints
#' @keywords internal
#' @importFrom checkmate assertNumber
.score_constraints_max <- function(M, clustersSizes, f = identity)
{
  if (is.infinite(M))
    return(0.0)

  grosClusters <- clustersSizes > M

  if (any(grosClusters))
  {
    score <- sum(f(clustersSizes[grosClusters] - M))
    assertNumber(score, lower = 0.0)
    return(as.double(score))
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
score_constraints <- function(clustersSizes = clusters_sizes(partition, sizes),
                              m = 0.0, M = Inf,
                              partition = NULL,
                              f = identity, details = FALSE)
{
  assertFlag(details)

  if (!"clustersSizes" %in% names(as.list(match.call())))
  {
    if (is.null(sizes) || is.null(partition))
    {
      stop(NEEDCLUSTERSIZESORPART)
    }
    assertPartition(partition)
    n <- length(partition)
    assertSizes(sizes, len = n, M = M)
    clustersSizes <- clusters_sizes(partition, sizes)
  }

  scoreContraintesMin <- .score_constraints_min(m, clustersSizes, f)
  scoreContraintesMax <- .score_constraints_max(M, clustersSizes, f)


  if (details)
    return(score_constraints_table(scoreContraintesMin,
                                   scoreContraintesMax))

  else
    return(scoreContraintesMin + scoreContraintesMax)
}
