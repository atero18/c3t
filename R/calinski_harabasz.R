#' Calinski-Harabasz Index (CHI / ICH)
#'
#' Methods to determine the value of the Calinski-Harabasz Index of a
#' partition with elements variable(s) compatible with euclidean distance.
#' This index is defined by between variance divided by within variance,
#' each divided by their degree of liberty.
#' Depending of `valueOnly`, the value returned by `calinski_harabasz` if
#' only the index or a list with more data. It can be useful for instance to
#' update the index after some change instead of start from the beginning.
#' @param data A data frame with one row per element, containing context.
#' Each variable included must be compatible with the euclidean distance
#' (and so have to be quantitative) because CHI use the variance.
#' @param partition The partition of the problem: a vector with a length
#' equal to the number of rows of `data`.
#' @param valueOnly `TRUE` if the only value needed is the index. `FALSE`
#' if the user wants more information.
#'
#' @returns Depends on `valueOnly`.
#' If `valueOnly = TRUE`, returns only the Calinski-Harabasz index
#' (positive double value).
#' If `valueOnly = FALSE`, returns a list containing the following elements:
#' * "`k`": number of clusters in the partition (strictly positive integer)
#' * "`nbElementsClusters`": strictly positive integers vector of size `k`
#' giving the number of elements of each cluster.
#' * "`withinClusters`": A positive double vector of size `k` containing
#' the within sum of squares of each cluster.
#' * "`W`": The within sum of squares of the partition, i.e. the sum of the
#' empirical cluster variances weighted by the number of elements of each
#' clusters. (positive real value)
#' * "`betweenClusters`": A positive double vector of size `k` containing for
#' each cluster the square euclidean norm between their and the global centroid
#' (i.e. ||c_i - c||²)
#' * "`B`" : The between sum of squares of the partition, i.e. the sum of the
#' square euclidean norm given in "`betweenClusters`" weighted by the number
#' of elements of each cluster.
#' * "`CHI`" : The Calinski-Harabasz Index (positive real value)
#' @details Some special cases can be considered:
#' * `n = 1` : their is only one possible partition so it is optimal by
#' definition. Hence we take `CHI = Inf`.
#' * `k = 1`: the unique cluster of the partition is the entire set.
#' Here the between sum of squares (B) is considered equal to 0. In the index
#' formula `B` is divided by `k-1` which equals `0` here. Except when `n = 1`
#' this partition is of no interest so we assign `CHI = 0`.
#' * `k = n` : each cluster is made of one element. Consequently the empirical
#' biased variance of each cluster is equal to 0, so is `W`. In the index
#' formula `W` is divided by `n - k` which is equal to 0 here. In that case
#' compactness is optimal but separation is at its worth. The index is set to
#' `NaN`.
#'
#' @seealso `calinhara()` from package
#' [fpc](https://cran.r-project.org/web/packages/fpc/index.html)
#' which calculates the CHI and return only its value. This function can be 2
#' or 3 times faster than `calinski_harabasz` but take a lot more of memory
#' (empirically the ratio "memory `fpc`" / "memory `c3t`" seems to increase
#' linearly with `n`, being around 80 at `n = 3500`).
#' @references Calinski, T., and Harabasz, J. (1974) A Dendrite Method for
#' Cluster Analysis, Communications in Statistics, 3, 1-27.
#' @keywords internal
#' @importFrom stats var
#' @importFrom checkmate assertFlag checkDataFrame checkMatrix assert
calinski_harabasz <- function(data, partition, valueOnly = TRUE) # nolint: cyclocomp_linter
{
  assert(checkDataFrame(data,
                        min.rows = 1L, min.cols = 1L,
                        any.missing = FALSE, all.missing = FALSE),
         checkMatrix(data,
                     min.rows = 1L, min.cols = 1L,
                     any.missing = FALSE, all.missing = FALSE))

  # Number of elements
  n <- nrow(data)

  # Number of variables
  p <- ncol(data)


  assertPartition(partition, n)
  k <- nbClusters(partition)

  assertFlag(valueOnly)
  if (valueOnly)
  {
    specialCase <- special_cases_criterion("CHI", n = n, k = k)
    if (!is.na.not.nan(specialCase))
      return(specialCase)
  }


  # Number of elements of each cluster
  nbElementsClusters <- table(partition)
  clustersNames <- names(nbElementsClusters)
  nbElementsClusters <- as.vector(nbElementsClusters)
  names(nbElementsClusters) <- clustersNames


  # Within Variance

  # If `k == 1` the inner variance is equal to the variance of
  # the entire data
  if (k == 1L)
  {
    Tot <- sum_squares(data, aggregate = TRUE)
    W <- withinClusters <- Tot
    names(withinClusters) <- partition[1L]
  }
  # If there is one element per cluster (i.e. k = n)
  # then the inner variance is equal to 0
  else if (k == n)
  {
    withinClusters <- double(n)
    names(withinClusters) <- rownames(data)
    W <- 0.0
  }
  else
  {
    # Clusters with one element have a variance of empirical biased variance
    # of 0
    # withinClusters <- vapply(splittedData, sum_squares,
    #                          double(1L), aggregate = TRUE)
    withinClusters <- as.vector(by(data, partition, sum_squares,
                                   aggregate = TRUE))

    W <- sum(withinClusters)
  }


  # Between variance
  if (valueOnly)
  {
    # If there is only one cluster then the empirical biased between
    # variance is equal to 0
    if (k == 1L)
      B <- 0.0

    # We use the relation between the total and betwwen + within sums of
    # squares
    else
    {
      Tot <- sum_squares(data, aggregate = TRUE)

      # Tot = B + W # nolint: commented_code_linter
      B <- Tot - W
    }
  }
  else
  {
    # Centroid of the entire data
    centroid <- colMeans(data)

    # Centroid of each cluster
    #
    # If there is only one cluster then its centroid is equal to
    # to global centroid
    if (k == 1L)
    {
      centroids <- centroid
      B <- betweenClusters <- 0.0
    }

    else
    {
      # If k == n the centroid of each cluster is equal to
      # their unique element
      if (k == n)
        centroids <- data[partition, ]
      else
      {
        centroids <- by(data, partition, colMeans)
        centroids <- do.call("rbind", centroids)
      }

      # Sum of squares between cluster centroids and global centroid
      if (p > 1L)
        betweenClusters <- apply(centroids, 1L,
                                 function(c) sum((c - centroid)^2L))

      else
        betweenClusters <- sum((centroids - centroid)^2L)

      B <- sum(nbElementsClusters * betweenClusters)
    }
  }

  # Special cases
  if (n == 1L)
    CHI <- Inf
  else if (k == 1L)
    CHI <- 0.0
  else if (k == n)
    CHI <- NaN


  else
    CHI <- (B / (k - 1L)) / (W / (n - k))

  if (valueOnly)
    return(CHI)

  names(betweenClusters) <- clustersNames
  names(withinClusters) <- clustersNames


  list(nbElementsClusters = nbElementsClusters,
       withinClusters = withinClusters,
       W = W,
       centroid = centroid,
       centroids = centroids,
       betweenClusters = betweenClusters,
       B = B,
       k = k,
       CHI = CHI)

}

#' Update the Calinski-Harabasz index
#' @inheritParams update_criterion_params
#' @keywords internal
update_calinski_harabasz <- function(dataCriterion, donor, receiver,
                                     givenElements,
                                     partitionBefore,
                                     dataElements)
{

  if (length(givenElements) == 0L)
    return(dataCriterion)

  n <- length(partitionBefore)

  donor <- as.character(donor)
  receiver <- as.character(receiver)

  oldDonorCardinal <- dataCriterion$nbElementsClusters[donor]
  oldReceiverCardinal <- dataCriterion$nbElementsClusters[receiver]
  nbGivenElements <- length(givenElements)
  sumDataGiven <- colSums(dataElements[givenElements, ])

  donorStillExists <- oldDonorCardinal > nbGivenElements

  # Update partition
  if (donorStillExists)
    elementsNewDonor <- setdiff(which(partitionBefore == donor), givenElements)

  elementsNewReceiver <- c(which(partitionBefore == receiver), givenElements)


  # Update within variance
  dataCriterion$W <- dataCriterion$W - dataCriterion$withinClusters[donor]

  if (donorStillExists)
  {
    if (oldDonorCardinal - nbGivenElements == 1L)
      dataCriterion$withinClusters[donor] <- 0.0
    else
    {
      dataCriterion$withinClusters[donor] <-
        sum_squares(dataElements[elementsNewDonor, ], aggregate = TRUE)

    }

    dataCriterion$W <- dataCriterion$W + dataCriterion$withinClusters[donor]
  }
  else
    dataCriterion$withinClusters <-
    dataCriterion$withinClusters[names(dataCriterion$withinClusters) != donor]

  dataCriterion$W <- dataCriterion$W -
    dataCriterion$withinClusters[receiver]

  dataCriterion$withinClusters[receiver] <-
    sum_squares(dataElements[elementsNewReceiver, ], aggregate = TRUE)


  dataCriterion$W <- dataCriterion$W + dataCriterion$withinClusters[receiver]

  dataCriterion$W <- unname(dataCriterion$W)

  # Update cardinals
  if (donorStillExists)
  {
    dataCriterion$nbElementsClusters[donor] <-
      dataCriterion$nbElementsClusters[donor] - nbGivenElements
  }
  else
  {
    dataCriterion$nbElementsClusters <-
      dataCriterion$nbElementsClusters[
        names(dataCriterion$nbElementsClusters) != donor]
  }

  dataCriterion$nbElementsClusters[receiver] <-
    dataCriterion$nbElementsClusters[receiver] + nbGivenElements


  # update between variance
  dataCriterion$B <- dataCriterion$B -
    oldDonorCardinal * dataCriterion$betweenClusters[donor]

  if (donorStillExists)
  {
    dataCriterion$centroids[donor, ] <-
      oldDonorCardinal * dataCriterion$centroids[donor, ] -
      sumDataGiven

    dataCriterion$centroids[donor, ] <-
      dataCriterion$centroids[donor, ] /
      (oldDonorCardinal - nbGivenElements)

    dataCriterion$betweenClusters[donor] <-
      sum((dataCriterion$centroids[donor, ] - dataCriterion$centroid)^2L)
    dataCriterion$B <- dataCriterion$B +
      (oldDonorCardinal - nbGivenElements) *
      dataCriterion$betweenClusters[donor]

  }
  else
  {
    dataCriterion$centroids <-
      dataCriterion$centroids[rownames(dataCriterion$centroids) != donor, ]
    dataCriterion$betweenClusters <-
      dataCriterion$betweenClusters[names(dataCriterion$betweenClusters) !=
                                      donor]
  }

  dataCriterion$B <- dataCriterion$B -
    oldReceiverCardinal * dataCriterion$betweenClusters[receiver]

  dataCriterion$centroids[receiver, ] <-
    oldReceiverCardinal * dataCriterion$centroids[receiver, ] +
    sumDataGiven

  dataCriterion$centroids[receiver, ] <-
    dataCriterion$centroids[receiver, ] /
    (oldReceiverCardinal + nbGivenElements)

  dataCriterion$betweenClusters[receiver] <-
    sum((dataCriterion$centroids[receiver, ] - dataCriterion$centroid)^2L)

  dataCriterion$B <- dataCriterion$B +
    (oldReceiverCardinal + nbGivenElements) *
    dataCriterion$betweenClusters[receiver]

  dataCriterion$B <- unname(dataCriterion$B)

  if (!donorStillExists)
    dataCriterion$k <- dataCriterion$k - 1L

  if (dataCriterion$k %in% c(1L, n))
    dataCriterion$CHI <- Inf
  else
    dataCriterion$CHI <- (dataCriterion$B / (dataCriterion$k - 1L)) /
    (dataCriterion$W / (n - dataCriterion$k))

  dataCriterion
}

#' Variance and sum of square
#' Tools for calculating variance and sum of squares
#' with a matrix or a data.frame.
#' @param data matrix or data.frame containing numeric values.
#' @param aggregate flag indicating if returned value should
#' be the sum of the values obtained per variable. Ignored
#' if `data` has only one column.
#' @name variance_tools
NULL

#' @describeIn variance_tools Return the empirical variance for each
#' column of `data`.
#' @param unbiased flag indicating if returned empirical variance
#' should be unbiased (sum of squares divided by n - 1, default) or
#' biased (divided by n).
#' @keywords internal
#' @importFrom checkmate assertFlag
colVars <- function(data, unbiased = TRUE, aggregate = FALSE)
{
  assertFlag(unbiased)

  # Number of elements
  n <- nrow(data)

  # Number of variables
  p <- ncol(data)


  if (n == 0L || p == 0L)
  {
    stop("`data` does not contain any value")
  }


  if (p == 1L)
    aggregate <- FALSE
  else
    assertFlag(aggregate)

  if (n == 1L)
  {
    if (unbiased)
      vars <- rep(NaN, p)
    else
      vars <- double(p)

    unbiased <- TRUE
  }

  else if (is.data.frame(data))
    vars <- vapply(data, var, double(1L))

  else if (is.matrix(data))
    vars <- apply(data, 2L, var)

  else
  {
    stop("unrecognized type for `data`")
  }


  if (aggregate)
    vars <- sum(vars)

  if (!unbiased)
    vars <- vars * (n - 1L) / n

  vars
}

#' @describeIn variance_tools Return the sum of squares for each column
#' of `data`.
#' @keywords internal
sum_squares <- function(data, aggregate = TRUE)
{
  n <- nrow(data)

  if (n == 1L)
    return(0.0)

  # `var` R function uses the unbiased estimator so correction
  # needs to be made
  (n - 1L) * colVars(data, unbiased = TRUE, aggregate = aggregate)
}
