#' @include utils.R
#' @include criteria.R
#' @include size_constraints.R
#' @include partition_tools.R

#' Create a reference class for representing a Partition.
#' @keywords internal

Partition <-
  setRefClass("Partition",
              fields = list(pb = "pbCon", partition = "numeric",
                            contraintes = "logical",
                            method = "character",
                            scoreSizeConsts = "numeric",
                            methodDetails = "list",
                            criteria = "list"),
              methods = list(
                #' Initialize a new Partition object.
                initialize = function(pb = NULL, partition = NULL,
                                      contraintes = NULL,
                                      scoreSizeConsts = NULL,
                                      contiguity = FALSE,
                                      method = "unk",
                                      methodDetails = list(),
                                      criteria = list())
                {
                  if (is.null(pb))
                    return()

                  .self$pb <- pb
                  .self$partition <- partition

                  .self$contraintes <- contraintes
                  .self$method <- method
                  .self$methodDetails <- methodDetails
                  .self$criteria <- criteria

                }
              ))


#' Verify the integrity of the `Partition` object.
#' @param object A `Partition` object.
#' @returns A character vector with an error message,
#' if any integrity check fails.
#' @noRd
#' @importFrom checkmate testString
verif_partition <- function(object)
{

  if (!testPartition(object$partition))
    return("partition must be a partition vector")

  if (!testString(object$method))
    return("method must be a character string")

  if (!is.vector(object$contraintes))
    return("contraintes must be a logical vector")

  if (!all(c("connectivity", "min", "max") %in% names(object$contraintes)))
    return(paste0("contraintes must contain (at least) states for the checks ",
                  "of connectivity, min, and max"))

}

#' Create a `Partition` using the provided parameters.
#' @param pb A `pbCon` object.
#' @param partition A numeric vector representing the partition.
#' @param method A character string representing the method used
#' for partitioning.
#' @param methodDetails A list containing details about the method used
#' for partitioning.
#' @param contiguity A logical value indicating if the partition satisfies
#' the contiguity constraint.
#' @param minConstraint A logical value indicating if the partition satisfies
#' the min size constraint.
#' @param maxConstraint A logical value indicating if the partition satisfies
#' the max size constraint.
#' @param scoreMinSizeConst A numeric value representing the minimum
#' score for size constraints. If `NA` this value is calculated.
#' @param scoreMaxSizeConst A numeric value representing the maximum
#' score for size constraints. If `NA` this value is calculated.
#' @param criteria A list containing criteria for the partition.
#' @returns A new Partition object.
#' @keywords internal
#' @importFrom checkmate assertClass assertFlag assertNumber
partition <- function(pb, partition, method, methodDetails = list(),
                      contiguity = NA, minConstraint = NA_real_,
                      maxConstraint = NA_real_,
                      scoreMinSizeConst = NA_real_,
                      scoreMaxSizeConst = NA_real_, criteria = list())
{

  # Handling constraints (validation or not)
  assertClass(pb, "pbCon")
  assertFlag(contiguity, na.ok = TRUE)
  assertFlag(minConstraint, na.ok = TRUE)
  assertFlag(maxConstraint, na.ok = TRUE)
  assertNumber(scoreMinSizeConst, lower = 0.0, na.ok = TRUE)
  assertNumber(scoreMaxSizeConst, lower = 0.0, na.ok = TRUE)


  partition <- standardize_partition(partition)

  # Handling constraints (validation or not)
  if (is.na(contiguity))
    contiguity <- pb$isRegionalisation(partition)

  if (is.na(minConstraint))
    minConstraint <- pb$nbTooSmallClusters(partition) == 0L

  if (minConstraint)
    scoreMinSizeConst <- 0.0


  if (is.na(maxConstraint))
  {
    maxConstraint <- pb$nbTooBigClusters(partition) == 0L
  }

  if (maxConstraint)
    scoreMaxSizeConst <- 0.0

  contraintes <- c(contiguity, minConstraint, maxConstraint)
  names(contraintes) <- c("connectivity", "min", "max")

  # Handling score for size constraints
  if (!minConstraint && is.na(scoreMinSizeConst))
    scoreMinSizeConst <-
    score_constraints_min(pb, partition = partition)

  if (!maxConstraint && is.na(scoreMaxSizeConst))
    scoreMaxSizeConst <-
    score_constraints_max(pb, partition = partition)


  scoreSizeConsts <-
    score_constraints_table(scoreMinSizeConst,
                            scoreMaxSizeConst)

  Partition$new(pb = pb, partition = partition, contraintes = contraintes,
                method = method, methodDetails = methodDetails,
                scoreSizeConsts = scoreSizeConsts,
                criteria = criteria)
}

#' @describeIn verif_cont_Partition Verify if the contiguity
#' constraint is satisfied.
#' @returns A logical value indicating if the constraint is satisfied.
#' @noRd
#' @name verif_cont_Partition
checkContiguityConst <- function(x)
{
  x$contraintes["connectivity"]
}

#' @importFrom methods setGeneric
setGeneric("checkContiguityConst")

#' @describeIn verif_cont_Partition Verify if the "min" constraint is
#' satisfied.
#' @noRd
checkMinSizeConst <- function(x)
{
  x$contraintes["min"]
}

#' @importFrom methods setGeneric
setGeneric("checkMinSizeConst")

#' @describeIn verif_cont_Partition Verify if the "max" constraint is
#' satisfied.
#' @noRd
checkMaxSizeConst <- function(x)
{
  x$contraintes["max"]
}

#' @importFrom methods setGeneric
setGeneric("checkMaxSizeConst")

#' Verify if constraints are satisfied for a `Partition` object
#'
#' @description Verify if all constraints are satisfied
#'
#' @param x A Partition object.
#' @name verif_cont_Partition
#' @noRd
verifContraintes <- function(x)
{
  all(x$contraintes)
}

setMethod(
  ".nbClusters",
  signature(partition = "Partition"),
  function(partition) .nbClusters(partition$partition)
)

Partition$methods(
  k = function() length(unique(partition))
)

#' Get the number of singleton clusters in the `Partition` object.
#' @param x A Partition object.
#' @returns The number of singleton clusters.
#' @noRd
setMethod(
  "nbSingletons",
  signature(partition = "Partition"),
  function(partition) nbSingletons(partition$partition)
)

#' @describeIn score_contraintes_Partition Get the min score.
#' @noRd
scoreMinSizeConst <- function(x)
{
  x$scoreSizeConsts["min"]
}

#' @importFrom methods setGeneric
setGeneric("scoreMinSizeConst")

#' @describeIn score_contraintes_Partition Get the max score.
#' @noRd
scoreMaxSizeConst <- function(x)
{
  x$scoreSizeConsts["max"]
}

#' @importFrom methods setGeneric
setGeneric("scoreMaxSizeConst")

#' Get the scores for size constraints in the `Partition` object.
#' @description Get the total score.
#' @param x A `Partition` object.
#' @name score_contraintes_Partition
#' @noRd
#' @keywords internal
scoreSizeConsts <- function(x)
{
  x$scoreSizeConsts["total"]
}

#' @importFrom methods setGeneric
setGeneric("scoreSizeConsts")

#' Get the sizes of clusters of the `Partition` object.
#' @inheritParams clusters_sizes
#' @param partition A Partition object.
#' @returns A numeric vector representing the sizes of clusters.
#' @keywords internal
setMethod(
  "clusters_sizes",
  signature(partition = "Partition", sizes = "ANY"),
  function(partition, sizes)
  {
    clusters_sizes(partition$partition,
                   sizes = partition$pb$sizes)
  }
)


setMethod(
  "all.equal",
  signature(target = "Partition", current = "Partition"),
  function(target, current) all.equal(target$partition, current$partition)
)

#' Number of elements in a partition
#'
#' Give the number of elements in a partition
#' @param x A `Partition` object
#' @returns the number of elements in the partition. (positive integer)
#' @name partition_length
#' @rdname partition_length
#' @keywords internal
NULL

#' @rdname partition_length
#' @keywords internal
setMethod("length", signature(x = "Partition"), function(x) length(x$partition))

#' Convert a list of partitions to a [tibble][tibble::tibble-package]
#' @param partitions A list of `Partition` objects.
#' @returns A [tibble][tibble::tibble-package] containing information
#' about each partition.
#' @keywords internal
#' @importFrom tibble tibble
Partition_list_to_tibble <- function(partitions)
{
  tibble(partition  = c3t_lapply(partitions, "[[", "partition"),
         nbClusters = c3t_sapply(partitions, .nbClusters, simplify = TRUE),
         checkContiguityConst = c3t_sapply(partitions,
                                           checkContiguityConst),
         checkMinSizeConst  = c3t_sapply(partitions,
                                         checkMinSizeConst),
         checkMaxSizeConst  = c3t_sapply(partitions,
                                         checkMaxSizeConst),
         scoreMinSizeConst = c3t_sapply(partitions,
                                        scoreMinSizeConst),
         scoreMaxSizeConst = c3t_sapply(partitions,
                                        scoreMaxSizeConst),
         scoreSizeConsts    = c3t_sapply(partitions,
                                         scoreSizeConsts)
  )
}
