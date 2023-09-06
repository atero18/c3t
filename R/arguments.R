# Class unions
#' @importFrom methods setClassUnion
setClassUnion("numericOrMissing", members = c("numeric", "missing"))
setClassUnion("functionOrNULL", members = c("function", "NULL"))
setClassUnion("dfOrNULL", members = c("data.frame", "NULL"))
setClassUnion("numericOrLogicalOrNULL",
              members = c("numeric", "logical", "NULL"))
setClassUnion("numericOrLogical", members = c("numeric", "logical"))
setClassUnion("matriceOrNULL", members = c("matrix", "NULL"))
setClassUnion("numericOrNULL", members = c("numeric", "NULL"))
setClassUnion("logicalOrNULL", members = c("logical", "NULL"))

#' @name check_general
#' @param .var.name Name of the checked object to print in assertions.
#' Defaults to the heuristic implemented in [checkmate::vname()].
#' @returns Please see [checkmate::checkMatrix()] help page (for instance)
#' for further details about the difference between check_, test_ and assert_
#' functions results.
#' @details Extension of the `checkmate` package.
NULL

# Find what element of `x` are true missing
#' values, i.e. `NA` but not `NaN`.
#' @keywords internal
is.na.not.nan <- function(x) is.na(x) & !is.nan(x)

#' @keywords internal
#' @importFrom checkmate assertFlag
checkInf <- function(x, positiveOnly = FALSE)
{
  assertFlag(positiveOnly)

  if (length(x) != 1L)
    return("argument must have length one.")

  if (!isTRUE(is.infinite(x)))
    return("argument is not infinite.")

  if (positiveOnly && x < 0.0)
    return("argument is negative.")

  TRUE
}

#' @keywords internal
checPosInf <- function(x) checkInf(x, positiveOnly = TRUE)

#' @keywords internal
#' @importFrom checkmate assert checkCount
assertCountOrInfinite <- function(x, positiveOnly = TRUE, ...)
{
  assert(checkInf(x, positiveOnly), checkCount(x, ...))
}

#' @importFrom checkmate checkDouble
checkProportions <- function(x, # nolint: cyclocomp_linter
                             zeros.ok = TRUE,
                             ones.ok = TRUE,
                             onlyZeros.ok = TRUE,
                             onlyOnes.ok = TRUE,
                             onlyExtremes.ok = TRUE,
                             any.missing = FALSE, all.missing = FALSE,
                             len = NULL, min.len = NULL, max.len = NULL,
                             unique = FALSE, sorted = FALSE,
                             null.ok = FALSE)
{

  if (null.ok && is.null(x))
    return(TRUE)

  assertFlag(zeros.ok)
  assertFlag(ones.ok)

  if (zeros.ok || ones.ok)
  {
    assertFlag(onlyExtremes.ok)
    if (onlyExtremes.ok)
    {
      assertFlag(onlyZeros.ok)
      assertFlag(onlyOnes.ok)
    }
    else
    {
      onlyZeros.ok <- FALSE
      onlyOnes.ok <- FALSE
    }
  }

  checkProp <- checkDouble(x, lower = 0.0, upper = 1.0,
                           any.missing = any.missing,
                           all.missing = all.missing,
                           len = len, min.len = min.len, max.len = max.len,
                           unique = unique, sorted = sorted,
                           null.ok = null.ok)

  if (!isTRUE(checkProp))
    return(checkProp)


  if ((all.missing && all(is.na(x))) || length(x) == 0L)
    return(TRUE)

  if (!zeros.ok || !onlyZeros.ok)
  {
    maskZeros <- x == 0.0
    if (!zeros.ok && any(maskZeros, na.rm = TRUE))
      return("argument contains at least one zero")
    if (all(maskZeros, na.rm = TRUE))
      return("argument is made only of zeros")

  }

  if (!ones.ok || !onlyOnes.ok)
  {
    maskOnes <- x == 1.0
    if (all(maskOnes, na.rm = TRUE))
      return("argument is made only of ones")
    else if (!ones.ok && any(maskOnes, na.rm = TRUE))
      return("argument contains at least one one")
  }

  if (!onlyExtremes.ok && all(maskZeros | maskOnes, na.rm = TRUE))
    return("argument is made only of zeros and ones")

  TRUE
}

#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertProportions <- makeAssertionFunction(checkProportions)

#' @keywords internal
#' @importFrom checkmate makeTestFunction
testProportions <- makeTestFunction(checkProportions)

#' Check if an argument is a subset of a character vector
#' @details Main differences with [checkmate::checkSubset()] is that
#' it specific to characters and allow before comparison to decide
#' if case and accents must be considered.
#' @param x The argument to be tested.
#' @param choices The non-empty character vectors `x` is supposed to
#' come.
#' @param len integer indicating the size `x` is supposed to have.
#' Default to `NULL` (do verification made).
#' @param keepCase `TRUE` if case must be maintained. If `FALSE` (default)
#' each letter in `x` strings are lowered, as well as those in `choices`.
#' @param keepAccents `TRUE` if accents must be maintained. If `FALSE`
#' (default) accents are removed in `x` and `choices` strings.
#' @inheritParams check_general
#' @inherit check_general return
#' @inherit check_general details
#' @name checkCharacterSubset
#' @rdname checkCharacterSubset
#' @keywords internal
NULL

#' @rdname checkCharacterSubset
#' @keywords internal
#' @importFrom checkmate assertCharacter checkCharacter checkSubset
checkCharacterSubset <- function(x, choices,
                                 keepCase = FALSE, keepAccents = FALSE,
                                 len = NULL, empty.ok = FALSE)
{
  if (!is.null(len))
    assertCount(len)

  empty.ok <- empty.ok && (is.null(len) || len == 0L)
  if (length(x) == 0L)
  {
    return(ifelse(empty.ok, TRUE, "argument must contains at least one value"))
  }

  assertCharacter(choices, min.len = 1L,
                  any.missing = FALSE, all.missing = FALSE)

  check <- checkCharacter(x,
                          len = len,
                          any.missing = FALSE, all.missing = FALSE)

  if (!isTRUE(check))
    return(check)

  if (!keepCase || !keepAccents)
  {
    x <- character_transformation(x, keepCase, keepAccents)
    choices <- character_transformation(choices, keepCase, keepAccents)
  }
  checkSubset(x, choices)
}

#' @rdname checkCharacterSubset
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertCharacterSubset <- makeAssertionFunction(checkCharacterSubset)

#' @rdname checkCharacterSubset
#' @keywords internal
#' @importFrom checkmate makeTestFunction
testCharacterSubset <- makeTestFunction(checkCharacterSubset)

#' @describeIn checkCharacterSubset Shortcut for when `len = 1L`,
#' i.e. `x` must be one element of `choices`.
#' @keywords internal
checkStringChoice <- function(x, choices,
                              keepCase = FALSE,
                              keepAccents = FALSE)
{
  checkCharacterSubset(x, choices,
                       keepCase = FALSE, keepAccents = FALSE,
                       len = 1L)
}


#' @rdname checkCharacterSubset
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertStringChoice <- makeAssertionFunction(checkStringChoice)

#' @rdname checkCharacterSubset
#' @keywords internal
#' @importFrom checkmate makeTestFunction
testStringChoice <- makeTestFunction(checkStringChoice)

#' Check if an argument is a contiguity matrix
#'
#' Check if an argument is a matrix with contiguity matrix properties.
#' @inheritParams checkDistanceMatrix
#' @inheritParams checkContiguityGraph
#' @param nrows,ncols exact number of rows and columns. If `NULL` check is
#' ignored. If `isComplete` is `TRUE` then `nrows` must be equals to `ncols`
#' @inheritParams check_general
#' @inherit check_general return
#' @inherit check_general details
#' @name checkContiguityMatrix
#' @rdname checkContiguityMatrix
#' @family arguments checkers
#' @export
#' @importFrom checkmate assertFlag assertCount checkMatrix
checkContiguityMatrix <- function(x, # nolint: cyclocomp_linter
                                  isComplete = TRUE,
                                  any.missing = FALSE,
                                  all.missing = FALSE,
                                  nrows = NULL,
                                  ncols = nrows)
{
  # Checking arguments
  assertFlag(isComplete)
  assertCount(nrows, null.ok = TRUE)
  assertCount(ncols, null.ok = TRUE)

  if (isComplete && (!is.null(nrows) && !is.null(ncols) && nrows != ncols))
  {
    stop("Complete matrix is expected but user ask for different number of rows and columns") # nolint: line_length_linter
  }


  if (inherits(x, "Matrix"))
    x <- as.matrix(x)

  check <- checkMatrix(x, mode = "logical",
                       any.missing = any.missing,
                       all.missing = all.missing)

  if (!isTRUE(check))
    return(check)

  if (isComplete && !isSymmetric(x))
    return("The argument must be symmetric when `isComplete = TRUE`")


  TRUE

}

#' @rdname checkContiguityMatrix
#' @export
#' @importFrom checkmate makeAssertionFunction
assertContiguityMatrix <- makeAssertionFunction(checkContiguityMatrix)

#' @rdname checkContiguityMatrix
#' @export
#' @importFrom checkmate makeTestFunction
testContiguityMatrix <- makeTestFunction(checkContiguityMatrix)

#' Check if an argument is a contiguity graph
#'
#' Check if an argument is a `ìgraph` object from `igraph` package
#' with contiguity graph properties, i.e. being simple and undirected.

#' @param g The argument to try.
#' @param n The order (number of vertices) `g` is supposed to have. By
#' default this number is arbitrary.
#' @param add Collection to store assertion messages.
#' See [checkmate::AssertCollection].
#' @inheritParams check_general
#' @inherit check_general return
#' @inherit check_general details
#' @name checkContiguityGraph
#' @rdname checkContiguityGraph
#' @family arguments checkers
#' @export
#' @importFrom igraph is_igraph is.simple is.directed gorder gsize
checkContiguityGraph <- function(g, n = gorder(g))
{
  if (!is_igraph(g))
    return("`g` is not an `igraph` object from `igraph` package.")

  else if (gorder(g) != n)
    return("`g` does not have the expected order.")

  else if (gorder(g) == 0L)
    return(TRUE)

  else if (gsize(g) == 0L)
    return(TRUE)

  else if (!is.simple(g))
    return("`g` must be a simple graph")

  else if (is.directed(g))
    return("`g` must be undirected")

  TRUE
}


#' @rdname checkContiguityGraph
#' @keywords internal
#' @export
#' @importFrom checkmate makeAssertionFunction
assertContiguityGraph <- makeAssertionFunction(checkContiguityGraph)

#' @rdname checkContiguityGraph
#' @keywords internal
#' @export
#' @importFrom checkmate makeTestFunction
testContiguityGraph <- makeTestFunction(checkContiguityGraph)

#' Check if an argument is a distance matrix
#'
#' @param x Matrix to be tested
#' @param isComplete A flag set to `TRUE` if `x` is supposed to be
#' a distance matrix of a complete set, i.e. must be symmetric with 0
#' at the diagonal (can be NA if `any.missing = TRUE`).
#' @param any.missing Are missing values allowed? Default is `FALSE`.
#' @param all.missing Are matrices with only missing values allowed?
#' Default is `FALSE`.
#' @name checkDistanceMatrix
#' @rdname checkDistanceMatrix
#' @inheritParams check_general
#' @inheritParams checkContiguityGraph
#' @inherit check_general return
#' @inherit check_general details
#' @family arguments checkers
#' @export
#' @importFrom checkmate checkMatrix
checkDistanceMatrix <- function(x, # nolint: cyclocomp_linter
                                isComplete = TRUE,
                                any.missing = FALSE,
                                all.missing = FALSE)
{

  if (inherits(x, "DistMat"))
  {
    if (!any.missing && anyNA(x))
      return("Argument contains missing values.")

    return(TRUE)
  }

  else if (inherits(x, "AbstractSymMat"))
  {
    if (!is.numeric(x$values))
      return("argument contains not numeric values")

    if (any(x$values < 0.0, na.rm = TRUE))
      return("All values are not positive.")

    if (!isTRUE(all(diag(x) == 0.0, na.rm = any.missing)))
      return("diagonal contains non 0 values")

    if (!any.missing && anyNA(x))
      return("Argument contains missing values.")

    return(TRUE)
  }

  check <- checkMatrix(x, mode = "numeric",
                       any.missing = any.missing,
                       all.missing = all.missing)

  if (!isTRUE(check))
    return(check)

  if  (any(x < 0.0, na.rm = TRUE))
    return("All values are not positive.")

  if (isTRUE(isComplete))
  {
    if (!isSymmetric(x))
      return("Matrix is not symetric")

    if (!isTRUE(all(diag(x) == 0.0, na.rm = any.missing)))
      return("diagonal contains non 0 values")
  }

  TRUE
}

#' @rdname checkDistanceMatrix
#' @keywords internal
#' @export
#' @importFrom checkmate makeAssertionFunction
assertDistanceMatrix <- makeAssertionFunction(checkDistanceMatrix)

#' @rdname checkDistanceMatrix
#' @keywords internal
#' @export
#' @importFrom checkmate makeTestFunction
testDistanceMatrix <- makeTestFunction(checkDistanceMatrix)


#' Check if an argument is a partition
#'
#' A partition is a Partition object or a non-empty integer-ish vector
#' with no missing values.
#' @param x Argument to be tested
#' @param n the expected length of the partition.
#' @name checkPartition
#' @rdname checkPartition
#' @inheritParams check_general
#' @inheritParams checkContiguityGraph
#' @inherit check_general return
#' @inherit check_general details
#' @family arguments checkers
#' @export
#' @importFrom checkmate checkIntegerish
checkPartition <- function(x, n = length(x))
{

  assertCount(n, positive = TRUE)

  if (length(x) != n)
    return("Size of `x` if different from the one asked.")

  if (inherits(x, "Partition"))
    return(TRUE)

  else if (!is.vector(x))
    return("argument must be a Partition object or a vector")


  checkIntegerish(x, lower = 1L,
                  any.missing = FALSE,
                  all.missing = FALSE,
                  min.len = 1L)
}

#' @rdname checkPartition
#' @keywords internal
#' @export
#' @importFrom checkmate makeAssertionFunction
assertPartition <- makeAssertionFunction(checkPartition)

#' @rdname checkPartition
#' @keywords internal
#' @export
#' @importFrom checkmate makeTestFunction
testPartition <- makeTestFunction(checkPartition)


#' Check if an argument is a regionalisation
#'
#' A regionalisation if a partition with all clusters being connected
#' (there only have one connected component each), regarding to a contiguity
#' relation.
#' @param regionalisation The argument to check.
#' @param contiguity A contiguity matrix or a contiguity graph (`igraph`).
#' @name checkRegionalisation
#' @rdname checkRegionalisation
#' @inheritParams check_general
#' @inheritParams checkContiguityGraph
#' @inherit check_general return
#' @inherit check_general details
#' @family arguments checkers
#' @export
#' @importFrom igraph gorder is_igraph
checkRegionalisation <- function(regionalisation, contiguity)
{

  if (is.matrix(contiguity))
    n <- nrow(contiguity)

  else if (is_igraph(contiguity))
    n <- gorder(contiguity)

  check <- checkPartition(regionalisation, n)

  if (!isTRUE(check))
    return(check)

  if (!isTRUE(is_regionalisation(regionalisation, contiguity)))
    return("at least one cluster is not connected")

  TRUE
}

#' @rdname checkRegionalisation
#' @keywords internal
#' @export
#' @importFrom checkmate makeAssertionFunction
assertRegionalisation <- makeAssertionFunction(checkRegionalisation)

#' @rdname checkRegionalisation
#' @keywords internal
#' @export
#' @importFrom checkmate makeTestFunction
testRegionalisation <- makeTestFunction(checkRegionalisation)

#' Checking is an argument is a criterion

#' Check if an argument is a correct clustering criterion.
#' @inheritParams check_general
#' @inheritParams checkContiguityGraph
#' @inherit check_general return
#' @inherit check_general details
#' @name checkCriteria
#' @rdname checkCriteria
NULL

#' @describeIn checkCriteria Check if criteria are actually available
#' @keywords internal
#' @importFrom checkmate assertFlag
checkCriteria <- function(criteria, unique = FALSE)
{
  assertFlag(unique)

  if (unique && length(criteria) != 1L)
    return("Argument must contains only one value.")

  else if (length(criteria) == 0L)
    return(TRUE)

  checkCharacterSubset(criteria,
                       c(CRITERIA$criterion, do.call("c", CRITERIA$names)))
}

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertCriteria <- makeAssertionFunction(checkCriteria)

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeTestFunction
testCriteria <- makeTestFunction(checkCriteria)

#' @describeIn checkCriteria Check if ONE criterion is actually available
#' @keywords internal
checkCriterion <- function(criterion) checkCriteria(criterion, unique = TRUE)

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertCriterion <- makeAssertionFunction(checkCriterion)

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeTestFunction
testCriterion <- makeTestFunction(checkCriterion)

#' @describeIn checkCriteria Check if criteria are compatible with given data.
#' @keywords internal
checkCompatibleCriteria <- function(criteria, # nolint: cyclocomp_linter
                                    distances = NULL,
                                    d = NULL,
                                    data = NULL,
                                    linkage = NULL,
                                    unique = FALSE)
{
  check <- checkCriteria(criteria, unique)

  if (!isTRUE(check))
    return(check)

  if (inherits(distances, "pbCon"))
  {
    data <- distances$data
    d <- distances$d
  }

  else if (!is.null(distances))
  {
    checkDistances <-
      checkDistanceMatrix(distances,
                          any.missing = TRUE,
                          all.missing = TRUE,
                          isComplete = TRUE)

    if (!isTRUE(checkDistances))
      return(checkDistances)
  }

  if (!is.null(d) && !is.function(d))
    return("When given `d` must be a function.")

  if (!is.null(data) && !is.data.frame(data))
    return("When given `data` must be a data frame.")

  if (!is.null(linkage))
  {
    checkLinkage <- checkLinkages(linkage)
    if (!isTRUE(checkLinkage))
      return(checkLinkage)

    linkage <- corresponding_linkage(linkage)
    allNA <- all(vapply(linkage, is.na, logical(1L)))
    if (allNA)
      linkage <- NULL
  }

  criteria <- unique(corresponding_criterion(criteria))

  dataCriteria <- CRITERIA[CRITERIA$criterion %in% criteria, ]

  if (any(dataCriteria$needsElemsDist))
  {
    if (is.null(distances) && (is.null(d) || is.null(data)))
    {
      return(paste0("Elements distances are needed for ",
                    toString(criteria[dataCriteria$needsd])))
    }

    if (!is.null(distances) &&
        anyNA(distances) &&
        (is.null(d) || is.null(data)))
    {
      return(paste0("Some elements distances are missing, can't be ",
                    "calculated and are needed for ",
                    toString(criteria[dataCriteria$needsd])))
    }
  }
  if (is.null(d) && any(dataCriteria$needsd))
  {
    return(paste0("elements distance function is needed for ",
                  toString(criteria[dataCriteria$needsd])))
  }


  if (is.null(data) && any(dataCriteria$needsData))
  {
    return(paste0("elements data is needed for ",
                  toString(criteria[dataCriteria$needsData])))
  }

  if (any(dataCriteria$needsQuantitativeData) && !all_numeric_variables(data))
  {
    return(paste0("elements data must be quantitative for ",
                  toString(criteria[dataCriteria$needsQuantitativeData])))

  }

  if (is.null(linkage) && any(dataCriteria$needsLinkage))
  {
    return(paste0("linkage is needed for ",
                  toString(criteria[dataCriteria$needsLinkage])))
  }

  TRUE
}

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertCompatibleCriteria <- makeAssertionFunction(checkCompatibleCriteria)

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeTestFunction
testCompatibleCriteria <- makeTestFunction(checkCompatibleCriteria)

#' @rdname checkCriteria
#' @keywords internal
checkCompatibleCriterion <-
  function(criterion, distances = NULL,
           d = NULL, data = NULL, linkage = NULL)
  {
    checkCompatibleCriteria(criterion, distances, d, data,
                            linkage, unique = TRUE)
  }

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertCompatibleCriterion <- makeAssertionFunction(checkCompatibleCriterion)

#' @rdname checkCriteria
#' @keywords internal
#' @importFrom checkmate makeTestFunction
testCompatibleCriterion <- makeTestFunction(checkCompatibleCriterion)


#' @importFrom checkmate assertFlag checkCharacter checkList
checkLinkages <- function(linkages, unique = FALSE) # nolint: cyclocomp_linter
{
  if (is.null(linkages))
    return(TRUE)

  assertFlag(unique)

  if (length(linkages) == 0L)
    return("No content given.")

  if (unique && length(linkages) > 1L)
  {
    return("Suppose to have only one linkage")

  }

  if (is.vector(linkages))
  {
    check <- checkCharacter(linkages,
                            min.len = 1L,
                            any.missing = TRUE,
                            all.missing = TRUE)

    if (!isTRUE(check))
      return(check)

  }

  else
  {
    check <- checkList(linkages,
                       types = c(character(1L), "function"),
                       min.len = 1L,
                       any.missing = TRUE,
                       all.missing = TRUE)

    if (!isTRUE(check))
      return(check)
  }

  corresponding <- corresponding_linkage(linkages, simplify = FALSE)

  nanMask <- vapply(corresponding, is.nan, logical(1L))

  if (any(nanMask))
    return("At least one linkage has not been recognized")

  TRUE
}

#' @importFrom checkmate makeAssertionFunction
assertLinkages <- makeAssertionFunction(checkLinkages)
#' @importFrom checkmate makeTestFunction
testLinkages <- makeTestFunction(checkLinkages)

checkLinkage <-
  function(linkage) checkLinkages(linkage, unique = TRUE)

#' @importFrom checkmate makeAssertionFunction
assertLinkage <- makeAssertionFunction(checkLinkage)
#' @importFrom checkmate makeTestFunction
testLinkage <- makeTestFunction(checkLinkage)
