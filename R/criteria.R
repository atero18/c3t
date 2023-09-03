#' @include arguments.R
#' @include tools.R
#' @include calinski_harabasz.R
#' @include dunn_index.R

#' @importFrom tibble tibble
CRITERIA <- tibble(criterion = c("CHI", "Dunn"), # nolint: object_name_linter
                   package = "c3t",
                   names = list(c("chi", "ich", "ch",
                                  "calinski\u002Dharabasz", "calinski"),
                                c("dunn", "dunn_index", "dunn_score")),
                   fun = c(calinski_harabasz, dunn_index),
                   updateFun = c(update_calinski_harabasz, update_dunn_index),
                   criterionListName = c("CHI", "D"),
                   lower = 0.0,
                   upper = Inf,
                   needsElemsDist = c(FALSE, TRUE),
                   needsd = c(FALSE, FALSE),
                   needsLinkage = c(FALSE, FALSE),
                   needsData = c(TRUE, FALSE),
                   needsQuantitativeData = c(TRUE, FALSE),
                   oneElemValue = Inf,
                   oneClusterValue = 0.0,
                   trivialPartitionValue = NaN,
                   optimality = "max",
                   handleAbstractSymMat = FALSE)


#' @importFrom tibble add_row
#' @importFrom checkmate assertString assertCharacter assertFunction
#' @importFrom checkmate assertNumber assertFlag
#' @keywords internal
add_criterion <- function(criterion,
                          package, names,
                          fun, updateFun,
                          criterionListName,
                          lower,
                          upper,
                          needsElemsDist,
                          needsd,
                          needsLinkage,
                          needsData,
                          needsQuantitativeData,
                          oneElemValue = NA_real_,
                          oneClusterValue = NA_real_,
                          trivialPartitionValue = NA_real_,
                          optimality = "max",
                          handleAbstractSymMat = FALSE)
{
  assertString(criterion, min.chars = 1L)
  criterion <- tolower(criterion)
  substr(criterion, 1L, 1L) <- toupper(substr(criterion, 1L, 1L))

  if (criterion %in% CRITERIA$criterion)
  {
    stop("This criterion already exists")
  }

  assertString(package, min.chars = 1L)
  assertCharacter(names, min.chars = 1L, min.len = 1L)
  names <- tolower(unique(names))

  assertFunction(fun)
  assertFunction(updateFun, null.ok = TRUE)

  if (!is.null(updateFun))
    assertString(criterionListName, min.chars = 1L)
  else
    criterionListName <- NA_character_

  assertNumber(lower)
  assertNumber(upper, lower = lower)
  assertFlag(needsElemsDist)
  assertFlag(handleAbstractSymMat)

  assertFlag(needsd)
  assertFlag(needsLinkage)

  if (needsd)
    needsData <- TRUE
  else
    assertFlag(needsData)

  if (needsData)
    assertFlag(needsQuantitativeData, na.ok = TRUE)
  else
    needsQuantitativeData <- FALSE

  assertNumber(oneElemValue,
               lower = lower, upper = upper,
               na.ok = TRUE)
  assertNumber(oneClusterValue,
               lower = lower, upper = upper,
               na.ok = TRUE)
  assertNumber(trivialPartitionValue,
               lower = lower, upper = upper,
               na.ok = TRUE)

  assertStringChoice(optimality, c("min", "max"))

  assign("CRITERIA",
         add_row(CRITERIA,
                 criterion = criterion,
                 package = package,
                 names = list(names),
                 fun = list(fun),
                 updateFun = list(updateFun),
                 criterionListName = criterionListName,
                 lower = lower,
                 upper = upper,
                 needsElemsDist = needsElemsDist,
                 needsd = needsd,
                 needsLinkage = needsLinkage,
                 needsData = needsData,
                 needsQuantitativeData = needsQuantitativeData,
                 oneElemValue = oneElemValue,
                 oneClusterValue = oneClusterValue,
                 trivialPartitionValue = trivialPartitionValue,
                 optimality = optimality),
         envir = asNamespace("c3t"))
}


if (requireNamespace("cluster", quietly = TRUE))
{
  #' silhouette_score
  #'
  #' Returns the Silhouette score of a partition using the inter-element
  #' distance matrix.
  #' The score used is the mean of all Silhouette coefficients
  #' (faster than aggregating by cluster).
  #'
  #' @param partition The partition of interest.
  #' @param distances The distance matrix (distance matrix)
  #' between elements.
  #' @param elements An optional subset of elements, for example, to speed up
  #' the calculation. Defaults to the complete set of elements.
  #'
  #' @returns The associated Silhouette score (floating-point value in [-1,1]).
  #' @noRd
  silhouette_score <- function(partition,
                               distances,
                               elements = seq_along(partition), ...)
  {
    n <- length(partition)
    nbClusters <- .nbClusters(partition)

    if (nbClusters %in% c(1L, n))
      return(1.0)

    if (length(elements) < n)
    {
      partition <- partition[elements]
      distances <- distances[elements, elements, drop = FALSE]
    }

    if (inherits(distances, "pbCon"))
      distances <- distances[]

    mean(cluster::silhouette(partition, dmatrix = distances)[, 3L],
         na.rm = TRUE)
  }

  add_criterion(criterion = "Silhouette",
                package = "cluster",
                names = c("silhouette", "s",
                          "silhouette_score", "silhouette_index"),
                fun = function(partition, distances)
                  silhouette_score(partition, distances),
                updateFun = NULL,
                criterionListName = NA_character_,
                lower = -1.0,
                upper = 1.0,
                optimality = "max",
                needsElemsDist = TRUE,
                needsd = FALSE,
                needsLinkage = FALSE,
                needsData = FALSE,
                needsQuantitativeData = FALSE,
                oneElemValue = 1.0,
                oneClusterValue = -1.0,
                trivialPartitionValue = 1.0,
                handleAbstractSymMat = FALSE)

}

#' @importFrom checkmate assertCharacter
corresponding_criterion <- function(criteria)
{
  assertCharacter(criteria, min.len = 1L)

  corresponding <- ifelse(criteria %in% CRITERIA$criterion,
                          criteria, NA_character_)


  if (anyNA(corresponding))
  {
    corresponding[is.na(corresponding)] <-
      CRITERIA$criterion[
        match_str(criteria[is.na(corresponding)], CRITERIA$names)]
  }

  corresponding
}

get_data_criterion <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  CRITERIA[CRITERIA$criterion == criterion, ]
}

get_function_criterion <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  CRITERIA$fun[CRITERIA$criterion == criterion][[1L]]
}

needs_elemDistances <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)
  CRITERIA$needsElemsDist[CRITERIA$criterion == criterion][[1L]]
}

creation_hdl_AbstractSymMat <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  !CRITERIA$needsElemsDist[CRITERIA$criterion == criterion] ||
    CRITERIA$handleAbstractSymMat[CRITERIA$criterion == criterion]
}

get_lower_criterion <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  CRITERIA$lower[CRITERIA$criterion == criterion]

}

get_upper_criterion <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  CRITERIA$upper[CRITERIA$criterion == criterion]

}

get_optimality_criterion <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  CRITERIA$optimality[CRITERIA$criterion == criterion]
}

get_update_function_criterion <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  CRITERIA$updateFun[CRITERIA$criterion == criterion][[1L]]
}

updatable_criterion <- function(criterion)
{
  !is.null(get_update_function_criterion(criterion))
}

get_list_name_criterion <- function(criterion)
{
  criterion <- corresponding_criterion(criterion)

  CRITERIA$criterionListName[CRITERIA$criterion == criterion]
}

#' Available clustering criteria
#'
#' Give what criteria are available actually.
#' @details Silhouette score will be available if `cluster` package
#' is installed.
#' @family available parameters
#' @export
available_criteria <- function()
{
  CRITERIA$criterion
}

#' @importFrom checkmate assertCount assertInt
special_cases_criterion <- function(criterion,
                                    partition,
                                    n = length(partition),
                                    k = .nbClusters(partition))
{
  criterion <- corresponding_criterion(criterion)
  assertCount(n, positive = TRUE)
  assertInt(k, lower = 1L, upper = n)

  dataCriterion <- CRITERIA[CRITERIA$criterion == criterion, ]


  if (n == 1L)
    dataCriterion$oneElemValue
  else if (k == 1L)
    dataCriterion$oneClusterValue
  else if (k == n)
    dataCriterion$trivialPartitionValue
  else
    NA_real_
}


#' @importFrom tibble tibble
#' @importFrom tidyr expand_grid
#' @importFrom dplyr distinct
criteria_linkages_grid <- function(criteria, linkages)
{
  criteria <- corresponding_criterion(criteria)
  criteria <- unique(criteria)

  if (!any(CRITERIA$needsLinkage[CRITERIA$criterion %in% criteria]))
  {
    grid <- tibble(criterion = criteria,
                   linkage = replicate(length(criteria), NULL,
                                       simplify = FALSE))
  }
  else
  {
    linkages <- corresponding_linkage(linkages)
    linkages <- unique(linkages)

    grid <- expand_grid(criterion = criteria, linkage = linkages)
    if (!all(CRITERIA$needsLinkage[CRITERIA$criterion %in% criteria]))
    {
      noLinkageCriteria <-
        setdiff(criteria, CRITERIA$criterion[!CRITERIA$needsLinkage])

      grid[grid$criterion %in% noLinkageCriteria, "linkage"] <- NULL
    }
  }

  grid <- distinct(grid)
  grid$name <- grid$criterion
  noLinkageLines <- vapply(grid$linkage, is.null, logical(1L))

  if (!all(noLinkageLines))
  {
    grid[!noLinkageLines, "name"] <- paste(grid$criterion[!noLinkageLines],
                                           grid$linkage[!noLinkageLines],
                                           sep = "_")
  }

  grid
}

#' Evaluating quality of partitions
#'
#' This function allow to evaluate some internal clustering
#' criteria on one or several partitions. Connectivity can be taken into
#' account for reduce the impact of contiguity / connectivity constraint.
#'
#' @param partitions a partition or a list of partitions of the same set
#' the criterion will be calculated on. (vector or list of vectors)
#' @param criterion one of the available criterion. Use [available_criteria()]
#' to see what criterion can be applied. (string)
#' @inheritParams arguments_problem
#' @inheritParams arguments_distMat
#' @param d Distance function between elements. Some criteria require this value
#' .If present, `data` must also be
#' specified. Some classical distances are available, it is recommended to use
#' them rather than a personal function for optimization reasons :
#' * "`euclidean`": euclidean distance.
#' * "`manhattan`" : manhattan distance.
#' * "`minkowski`" : minkowski distance. In that case a value for p >= 1
#' must be specified.
#' (function or string)
#' @param data  a data.frame where each row represents data related to an
#' element. This can be omitted if `d` is omitted, but might be necessary for
#' some criteria (e.g. Calinski-Harabasz). The present variables can
#' be quantitative or qualitative. If qualitative variables are present,
#' some distances and criteria may not be used. Possibility of standardising
#' variables and transforming qualitative variables into binary
#' variables (one-hot encoding)
#' using `standardQuant` and `binarQual`. (data.frame)
#' @param linkage a distance linkage. Can be a string
#' (see [available_linkages()]) or a user function. Used for some of the
#' criteria (e.g. Dunn).
#' @param connected a flag equals to `TRUE` if the criterion should be
#' calculated in its connected form, i.e. the mean of its value on each
#' connected component.
#' @param ... Arguments specific for the criterion.
#'
#' @importFrom checkmate assertFlag assert
#' @importFrom igraph is_igraph is_connected
#' @seealso [available_criteria()]
#' @export
clustering_criterion <- function(partitions, # nolint: cyclocomp_linter
                                 criterion,
                                 distances = NULL,
                                 d = NULL,
                                 data = NULL,
                                 standardQuant = FALSE,
                                 binarQual = FALSE,
                                 linkage = NULL,
                                 contiguity = NULL,
                                 connected = FALSE,
                                 ...)
{
  assertFlag(connected)

  if (!is.null(data))
  {
    assertFlag(standardQuant)
    assertFlag(binarQual)

    if (standardQuant || binarQual)
      data <- normalize_df(data, standardQuant, binarQual)
  }

  if (is.list(partitions))
  {
    lapply(partitions, assertPartition)
    nbPartitions <- length(partitions)
    if (nbPartitions == 1L)
    {
      partition <- partitions[[1L]]
      n <- length(partition)
    }

    else
    {
      lengthsPartitions <- lengths(partitions)
      n <- lengthsPartitions[1L]
      if (any(lengthsPartitions != n))
      {
        stop("All partitions must have the same number of elements")
      }
    }
  }

  else
  {
    nbPartitions <- 1L
    partition <- partitions
    assertPartition(partition)
    n <- length(partition)
  }

  assertCompatibleCriterion(criterion,
                            distances,
                            d, data,
                            linkage)


  if (needs_elemDistances(criterion) && anyNA(distances))
    diag(distances) <- 0.0

  if (needs_elemDistances(criterion) && anyNA(distances))
  {

    if (inherits(distances, "DistMat"))
      distances$calcul_missing_distances()
    else
    {
      indexesNA <- which_na(distances,
                            removeSymmetry = TRUE,
                            removeDiagonal = TRUE)

      naDistances <- calculate_distances(indexesNA, d, data)
      distances[indexesNA] <- naDistances
      distances[indexesNA[, 2L:1L, drop = FALSE]] <- naDistances
    }
  }

  if (!creation_hdl_AbstractSymMat(criterion) &&
      inherits(distances, "AbstractSymMat"))
  {
    if (inherits(distances, "DistMat"))
      distances <- distances$distances

    if (inherits(distances, "SymMMat"))
      distances <- distances$values
    else
      distances <- distances[]
  }

  if (connected)
  {
    assert(checkContiguityMatrix(contiguity, isComplete = TRUE, nrows = n),
           checkContiguityGraph(contiguity, n))


    if (!is_igraph(contiguity))
      contiguity <- contiguity_matrix_to_graph(contiguity)

    if (!is_connected(contiguity))
    {
      components <- components(contiguity)$membership
      nbConnectedComponents <- max(components)
      indexesCC <- split(seq_len(n), components)

      if (nbPartitions == 1L)
        partitionsCC <- split(partition, components)
      else
      {
        partitionsCC <- lapply(partitions, split, f = components)
        changeOrderList <- function(i) lapply(partitionsCC, "[[", i)
        partitionsCC <- lapply(seq_len(nbConnectedComponents), changeOrderList)
      }



      if (!is.null(distances))
      {
        distancesCC <- lapply(indexesCC,
                              function(indexes) distances[indexes,
                                                          indexes,
                                                          drop = FALSE])
      }
      else
        distancesCC <- replicate(nbConnectedComponents, NULL, simplify = FALSE)

      if (!is.null(data))
        dataCC <- split(data, components)

      else
        dataCC <- replicate(nbConnectedComponents, NULL, simplify = FALSE)

      criterionCC <-
        function(i) clustering_criterion(partitionsCC[[i]],
                                         criterion,
                                         distancesCC[[i]],
                                         d,
                                         dataCC[[i]],
                                         FALSE,
                                         FALSE,
                                         linkage,
                                         contiguity = NULL,
                                         connected = FALSE)

      res <- vapply(seq_len(nbConnectedComponents),
                    criterionCC,
                    double(nbPartitions))

      if (nbPartitions == 1L)
        return(mean(res))
      else
        return(colMeans(res))

    }
  }

  if (nbPartitions > 1L)
  {
    return(vapply(partitions, clustering_criterion,
                  double(1L),
                  criterion = criterion,
                  distances = distances,
                  d = d,
                  data = data,
                  standardQuant = FALSE,
                  binarQual = FALSE,
                  linkage = linkage,
                  contiguity = NULL,
                  connected = FALSE,
                  ...))
  }

  kwargs <- list(...)

  specialCase <- special_cases_criterion(criterion, partition)
  if (!is.na(specialCase))
    return(specialCase)

  partition <- standardize_partition(partition)

  funCrit <- get_function_criterion(criterion)


  if ("valueOnly" %in% names(kwargs))
    valueOnly <- !isFALSE(kwargs[["valueOnly"]]) # nolint: object_usage_linter
  else
    valueOnly <- TRUE # nolint: object_usage_linter

  eval(body(funCrit), envir = environment())
}

#' Update criterion parameters
#' @name update_criterion_params
#' @param partitionBefore the partition before the exchange of elements.
#' @param dataCriterion a list containing all data about
#' the actual criterion value.
#' @param donor the id of the donor cluster.
#' @param receveir the id of the receiver cluster.
#' @param givenElement element of `donor` given to
#' `receveir`. (strictly positive integer)
#' @param givenElements set of the elements of `donor` given to
#' `receveir`. (vector of length at least 1 made of strictly positive integers)
#' @param dataElements data about the elements. Corresponds to `data` in
#' other `c3t` functions.
#' @keywords internal
NULL

#' @inheritParams update_criterion_params
#' @inheritParams arguments_problem
#' @keywords internal
update_criterion <- function(partitionBefore,
                             criterion,
                             dataCriterion,
                             donor, receiver,
                             givenElements,
                             distances = NULL,
                             d = NULL,
                             dataElements = NULL,
                             linkage = NULL)

{

  # Can be use in special cases ont onlu one element can be transferred
  givenElement <- givenElements # nolint: object_use_linter

  funUpdate <- get_update_function_criterion(criterion)


  eval(body(funUpdate), envir = environment())
}


#' Sort partitions by a clustering criterion
#'
#' Returns the order of the best partitions according to a certain criterion.
#'
#' @param criterionValues Vector of numerical data indicating
#' the `criterion` values for different partitions. (numeric vector)
#' @param criterion Name of the clustering criterion. (character string)
#' @param byName `TRUE` (default) if what should be returned is the names of
#' partitions rather than their identifier. `FALSE` otherwise. The partition
#' names can be specified via the `names` variable. If it's `NULL`, then the
#' names of the `criterionValues` vector will be used.
#' If none are available, `byName` is set to `FALSE`. (boolean)
#' @param names Names of the partitions. Used when `byName` is `TRUE`,
#' ignored otherwise. If used, it must be the same size as
#' `criterionValues`. (vector or `NULL`)
#' @param keepNA flag indicating if NA values should be kept in the result.
#' If `TRUE`, they are kept at the end of the vector. If `FALSE`, they are
#' removed.
#' @returns A vector of the same size as `criterionValues`, indicating
#' the order (from best to worst) of partitions relative to the
#' requested criterion. The indicated data are the identifiers of partitions
#' if `byName = FALSE`, their names otherwise. (vector)
#' @export
#' @importFrom checkmate assertFlag assertDouble
optimal_partitions <- function(criterionValues, criterion, byName = TRUE,
                               names = NULL, keepNA = TRUE)
{
  if (length(criterionValues) == 0L)
    return(integer(0L))

  # Checking aruments
  assertCriterion(criterion)
  assertFlag(byName)

  lowerValueCrit <- get_lower_criterion(criterion)
  upperValueCrit <- get_upper_criterion(criterion)


  assertDouble(criterionValues,
               lower = lowerValueCrit,
               upper = upperValueCrit,
               any.missing = TRUE,
               all.missing = FALSE,
               null.ok = FALSE)

  if (anyNA(criterionValues))
    assertFlag(keepNA)

  if (byName && is.null(names) && !is.null(names(criterionValues)))
    names <- names(criterionValues)

  byName <- byName && !is.null(names)

  # Ordonnancement depending of the criterion
  optimalityCrit <- get_optimality_criterion(criterion)

  decreasing <- switch(optimalityCrit,
                       "max" = TRUE,
                       "min" = FALSE,
                       stop("unrecognized optimality type for the criterion"))

  na.last <- ifelse(isTRUE(keepNA), TRUE, NA)
  order <- order(criterionValues, na.last = na.last, decreasing = decreasing)

  if (byName)
    order <- names[order]

  order
}
