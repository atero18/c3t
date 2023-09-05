#' @include abstractSymMat.R


#' Create a reference class for `DistMat`, representing the distance matrix
#'
#' @field distances An `AbstractSymMat` object representing the distance matrix.
#' @field d A function to calculate distance values.
#' @field data A dataframe containing the data points.
#' @keywords internal
distMat <-
  setRefClass("DistMat",
              fields = list(distances = "AbstractSymMat",
                            d = "functionOrNULL",
                            data = "dfOrNULL"))

distMat$methods(
  which_na = function(subsetRow = seq_len(nrow(.self)),
                      subsetCol = seq_len(nrow(.self)))
  {
    if (is.list(subsetRow))
    {
      if (length(subsetRow) == 1L)
      {
        subsetCol <- subsetRow[[1L]][[2L]]
        subsetRow <- subsetRow[[1L]][[1L]]
      }
      else if (length(subsetRow) == 0L)
        return(matrix(nrow = 0L, ncol = 2L))
    }

    if (is.list(subsetRow))
    {
      subsetCol <- lapply(subsetRow, "[[", 2L)

      subsetRow <- lapply(subsetRow, "[[", 1L)

      na <- lapply(seq_along(subsetCol),
                   function(k) distances$which_na(subsetRow[[k]],
                                                  subsetCol[[k]],
                                                  nanIsNA = FALSE,
                                                  simplify = FALSE))

      na <- do.call("rbind", na)
      if (nrow(na) == 0L)
        return(na)

      maskSuperior <- na[, 1L] > na[, 2L]
      if (any(maskSuperior))
        na[maskSuperior, ] <- na[maskSuperior, 2L:1L, drop = FALSE]

      return(unique(na, MARGIN = 1L))
    }

    distances$which_na(subsetRow, subsetCol, nanIsNA = FALSE, simplify = TRUE)
  }
)

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("names", signature = "DistMat", function(x) names(x$distances))
# Few shortcuts to access to the property of the distances values

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("length", signature = "DistMat", function(x) length(x$distances))

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("dim", signature = "DistMat", function(x) dim(x$distances))
DISTANCSEPROPERTIES <- c("names", "length", "dim", "nrow", "ncol")

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("nrow", signature = "DistMat", function(x) nrow(x$distances))

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("ncol", signature = "DistMat", function(x) ncol(x$distances))



#' Check if the function `d` is present to calculate distance values
#' @param x A `DistMat` object
#' @returns TRUE if the function d is present, FALSE otherwise
#' @noRd
presenceFonctionDistance <- function(x) !is.null(x$d)


distMat$methods(
  calcul_needed_distances = function(subsetRow, subsetCol = NULL)
  {
    na <- which_na(subsetRow, subsetCol)

    if (nrow(na) == 0L)
      return(invisible(TRUE))

    .self[na] <- calculate_distances(na, d, data)

    invisible(TRUE)
  },
  calcul_missing_distances = function()
  {
    na <- which_na()

    if (nrow(na) == 0L)
      return(invisible(TRUE))

    .self[na] <- calculate_distances(na, d, data)

    invisible(TRUE)
  }
)

#' Access to `DistMat` data`
#'
#' Rules are the following :
#' * when rows and columns are given a submatrix is returned with the
#' asked rows and columns
#' * when only rows (resp. columns are given) a submatrix is return with
#' the asked rows (resp. columns) and the entire columns (resp. rows)
#' * When rows and columns are missing the entire matrix is sent back.
#' Depending on `drop` (a flag) is there is only one row and / or
#' one column given the result will be simplified in a vector form or not.
#' If some values are missing they will be calculated if possible with the `d`
#' distance function with the `data` data frame.
#' @name distMat-access
#' @rdname distMat-access
#' @param i A vector composed of the indices of the wished rows.
#' @param j A vector composed of the indices of the wished columns
#' @param x An `AbstractSymMat` or inheriting object
#' @returns Depending on `drop` and the dimensions of `i` and `j`, a matrix or
#' a vector.
NULL

#' @rdname distMat-access
#' @keywords internal
setMethod(
  "[",
  signature(x = "DistMat", i = "matrix", j = "missing", drop = "ANY"),
  function(x, i, j, drop)
  {
    if (missing(drop))
      drop <- TRUE
    else
      drop <- isTRUE(drop)

    res <- x$distances[i, drop = drop]

    posNA <- is.na(res)
    if (any(posNA))
    {
      naIndexes <- i[posNA, , drop = FALSE]
      naDistances <- calculate_distances(naIndexes,
                                         x$d, x$data)

      x[naIndexes] <- naDistances
      res[posNA] <- naDistances
    }

    res

  }
)

distMat_vectorial_access <- function(x, i, j, drop)
{
  drop <- ifelse(missing(drop), TRUE, isTRUE(drop))

  if (missing(i))
    i <- seq_len(nrow(x))

  if (missing(j))
    j <- seq_len(nrow(x))

  res <- x$distances[i, j, drop = drop]

  if (presenceFonctionDistance(x))
  {
    posNA <- is.na(res)
    if (any(posNA))
    {

      if (length(res) == 1L)
        res <- calculate_distances(matrix(c(i, j), ncol = 2L),
                                   x$d, x$data)

      else
      {
        if (length(i) == 1L || length(j) == 1L)
        {
          posNA <- is.na(res)
          indexs <- matrix(NA, nrow = sum(posNA), ncol = 2L)
          if (length(i) == 1L)
          {
            indexs[, 1L] <- i
            indexs[, 2L] <- j[posNA]
          }
          else
          {
            indexs[, 1L] <- j
            indexs[, 2L] <- i[posNA]
          }
        }
        else
        {
          posNA  <- which(is.na(res), arr.ind = TRUE)
          indexs <- matrix(NA, nrow = nrow(posNA), ncol = 2L)
          indexs[, 1L] <- i[posNA[, 1L]]
          indexs[, 2L] <- j[posNA[, 2L]]
        }

        naDistances <- calculate_distances(indexs, x$d, x$data)
        x[indexs] <- naDistances
        res[posNA] <- naDistances
      }
    }
  }

  return(res)
}

#' @rdname distMat-access
#' @keywords internal
setMethod(
  "[",
  signature(x = "DistMat", i = "numericOrMissing",
            j = "numericOrMissing", drop = "ANY"),
  distMat_vectorial_access
)

#' @rdname distMat-access
#' @keywords internal
setMethod(
  "[",
  signature(x = "DistMat", i = "list", j = "missing", drop = "ANY"),
  function(x, i, j, drop)
  {
    drop <- ifelse(missing(drop), TRUE, isTRUE(drop))

    k <- length(i)
    distances <- x$distances$values
    presenceNA <- anyNA(distances)

    matrices <- lapply(i, function(l) distances[l$i, l$j, drop = FALSE])

    if (presenceNA)
    {
      elementsNA <-
        lapply(matrices, function(sM) which(is.na(sM), arr.ind = TRUE))

      sousMatricesAvecNA <- which(lengths(elementsNA) > 0L)
      for (k in sousMatricesAvecNA)
      {
        indexs <- elementsNA[[k]]
        indexs[, 1L] <- (i[[k]]$i)[indexs[, 1L]]
        indexs[, 2L] <- (i[[k]]$j)[indexs[, 2L]]
        naDistances <- calculate_distances(indexs, x$d, x$data)
        x[indexs] <- naDistances
        matrices[[k]][elementsNA[[k]]] <- naDistances
      }
    }

    return(matrices)
  }
)

#' `DistMat` values replacement
#' Call the `AbstractSymMat` replacement method. Check that values
#' are all positive and if `NA` replace them by `NA_real_`
#' @keywords internal
#' @name distMat_replace_values
NULL

#' @describeIn distMat_replace_values case when row(s) and column(s)
#' are given
#' @keywords internal
#' @importFrom checkmate assertDouble
setReplaceMethod(
  "[",
  signature(x = "DistMat", i = "numeric", j = "numeric", value = "numeric"),
  function(x, i, j, value) {
    assertDouble(value, lower = 0.0,
                 any.missing = TRUE,
                 all.missing = TRUE)

    if (anyNA(value) &&
        sum(is.na(value) & !is.nan(value)) > 0L)
    {
      value[is.na(value) & !is.nan(value)] <- NA_real_
    }

    x$distances[i, j] <- value


    return(x)
    }
)

#' @describeIn distMat_replace_values case when only row(s) are
#' given
#' @keywords internal
#' @importFrom checkmate assertDouble
setReplaceMethod(
  "[",
  signature(x = "DistMat", i = "numeric", j = "missing", value = "numeric"),
  function(x, i, j, value)
  {
    assertDouble(value, lower = 0.0,
                 any.missing = TRUE,
                 all.missing = TRUE)

    if (anyNA(value) &&
        sum(is.na(value) & !is.nan(value)) > 0L)
    {
      value[is.na(value) & !is.nan(value)] <- NA_real_
    }

    x$distances[i, ] <- value
    return(x)
  }
)

#' @describeIn distMat_replace_values case when only columns are given
#' @keywords internal
#' @importFrom checkmate assertDouble
setReplaceMethod(
  "[",
  signature(x = "DistMat", i = "missing", j = "numeric", value = "numeric"),
  function(x, i, j, value)
  {
    assertDouble(value, lower = 0.0,
                 any.missing = TRUE,
                 all.missing = TRUE)

    if (anyNA(value) &&
        sum(is.na(value) & !is.nan(value)) > 0L)
    {
      value[is.na(value) & !is.nan(value)] <- NA_real_
    }

    x$distances[, j] <- value
    return(x)
  }
)

#' @describeIn distMat_replace_values case when indexes are given
#' by pairs in a matrix
#' @keywords internal
#' @importFrom checkmate assertDouble
setReplaceMethod(
  "[",
  signature(x = "DistMat", i = "matrix", j = "missing", value = "numeric"),
  function(x, i, j, value)
  {
    assertDouble(value, lower = 0.0,
                 any.missing = TRUE,
                 all.missing = TRUE)

    if (anyNA(value) &&
        sum(is.na(value) & !is.nan(value)) > 0L)
    {
      value[is.na(value) & !is.nan(value)] <- NA_real_
    }

    x$distances[i] <- value
    return(x)
  }
)

#' @describeIn distMat_replace_values when the entire distance matrix is
#' replaced by another matrix
#' @keywords internal
#' @importFrom checkmate assertDouble
setReplaceMethod(
  "[",
  signature(x = "DistMat", i = "missing", j = "missing", value = "matrix"),
  function(x, i, j, value)
  {
    assertDouble(value, lower = 0.0,
                 any.missing = TRUE,
                 all.missing = TRUE)

    if (anyNA(value) &&
        sum(is.na(value) & !is.nan(value)) > 0L)
    {
      value[is.na(value) & !is.nan(value)] <- NA_real_
    }

    x$distances[] <- value
    return(x)
  }
)


#' @importFrom methods setGeneric
setGeneric("diag<-")

#' @rdname abstractSymMat_replace_diag
#' @keywords internal
setReplaceMethod(
  "diag",
  signature(x = "DistMat", value = "ANY"),
  function(x, value)
  {
    diag(x$distances) <- value
    return(x)
  }
)

#' @importFrom methods setGeneric
setGeneric("diag")

#' Access to diagonal of a `DistMat` object. The result is constantly
#' a double vector of length the number of elements.
#' @keywords internal
#' @name distMat_access_diag
setMethod("diag", signature(x = "DistMat"), function(x) double(nrow(x)))

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("isSymmetric", signature = "DistMat", function(object, ...) TRUE)

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("is.matrix", signature = "DistMat", function(x) TRUE)

#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("is.numeric", signature = "DistMat", function(x) TRUE)


#' Is there missing values in the distance matrix ?
#' @param x a `DistMat` object
#' @name distMat_NA
#' @rdname distMat_NA
#' @keywords internal
NULL

#' @describeIn distMat_NA Check if the `DistMat` object contains any NA values
#' @keywords internal
setMethod("anyNA", signature = "DistMat", function(x) anyNA(x$distances))

#' @describeIn distMat_NA Check what values are NA
#' @keywords internal
setMethod("is.na", signature = "DistMat", function(x) is.na(x$distances))



#' Custom all.equal method for comparing `DistMat` objects with matrices
#' @param target,current two variables, one of them at least being a `DistMat`.
#' @returns TRUE if the target and current are equal, FALSE otherwise
#' @name distMat_alleq
#' @keywords internal
#' @importFrom methods validObject setMethod
all.equal.DistMat.Matrix <- function(target, current)
{
  validObject(target)
  all(target[] == current)
}

#' @name distMat_alleq
#' @keywords internal
setMethod(
  "all.equal",
  signature(target = "DistMat", current = "matrix"),
  all.equal.DistMat.Matrix
)

#' @describeIn distMat_alleq matrix + `DistMat`
setMethod(
  "all.equal",
  signature(target = "matrix", current = "DistMat"),
  function(target, current, ...) all.equal.DistMat.Matrix(current, target)
)

#' @describeIn distMat_alleq `DistMat` + `DistMat`
#' @keywords internal
setMethod(
  "all.equal",
  signature(target = "DistMat", current = "DistMat"),
  function(target, current, ...) all.equal(target$distances, current$distances)
)

#' Convert a matrix to a `DistMat` object
#' @param from A matrix
#' @returns A `DistMat` object with distances stored as a symmetric matrix
#' @noRd
#' @importFrom methods as setAs
setAs(
  "matrix", "DistMat",
  function(from)
  {
    distances <- as(from, "SymMMat")
    distMat$new(distances = distances, d = NULL, data = NULL)
  }
)

#' Convert a `dist` object to a `DistMat` object
#' @param from A `dist` object
#' @returns A `DistMat` object with distances stored as a vector
#' @noRd
#' @importFrom stats dist
#' @importFrom methods setAs
setAs(
  "dist", "DistMat",
  function(from)
  {
    distances <- as(from, "SymVMat")
    distMat$new(distances = distances, d = NULL, data = NULL)
  }
)


#' Convert a `DistMat` object to a `dist` object
#' @param from A `DistMat` object
#' @returns A `dist` object with distances extracted from the `DistMat` object
#' @importFrom methods as setAs
#' @importFrom stats dist
#' @name distMat_to_dist
#' @noRd
setAs("DistMat", "dist", function(from) as(from$distances, "dist"))

#' @describeIn distMat_to_dist Convert a `DistMat` object to a dist object
#' with options for diag and upper.
#' @param m A `DistMat` object
#' @param diag, upper, ... Additional parameters
#' @noRd
#' @importFrom stats as.dist
#' @importFrom methods as
setMethod("as.dist",
          signature(m = "DistMat", diag = "ANY", upper = "ANY"),
          function(m, diag, upper) as(m, "dist"))

#' Constructor function for creating a `DistMat` object
#'
#' @param distances The distance matrix of the problem. This can be omitted
#' if a distance function `d` and data context `data` are provided. If only
#' `distances` is provided, all distances must be present. (distance matrix)
#' @param d Distance function between elements. This can be omitted if
#' `distances` is already indicated. If present, `data` must also be
#' specified. Some classical distances are available, it is recommended to use
#' them rather than a personal function for optimisation reasons :
#' * "`euclidean`": Euclidean distance.
#' * "`manhattan`" : Manhattan distance.
#' * "`minkowski`" : Minkowski's distance. In that case a value for p >= 1
#' must be specified.
#'
#' (function or string)
#' @param data A data.frame where each row represents data related to an
#' element. This can be omitted if `d` is omitted. Present variables can
#' be quantitative or qualitative. If qualitative variables are present,
#' some distances may not be used. Possibility of standardising variables and
#' transforming qualitative variables into binary variables (one-hot encoding)
#' using `standardQuant` and `binarQual`. (data.frame)
#' @param standardQuant `TRUE` if the variables in `data` should be
#' standardised (i.e., centered and scaled), `FALSE` (default) otherwise.
#' Standardisation is applied after the possible binarization of qualitative
#' variables (see `binarQual`). (flag)
#' @param binarQual `TRUE` if qualitative variables should be binarized (one-hot
#' encoding), for example, to make the data set compatible with common distances
#' or to standardize these variables. `FALSE` (default) otherwise. (flag)
#' @param storageMode Indicates how distance data should be stored. Do not
#' modify for now. (string)
#' @param names A character vector with names for the elements in
#' the distance matrix.
#' @param p An integer value for the Minkowski distance parameter.
#' @returns A `DistMat` object representing the distance matrix.
#' @keywords internal
#' @name arguments_distMat
#' @importFrom methods as
#' @importFrom checkmate testString
constructor_DistMat <- function(distances = NULL, d = NULL, data = NULL, # nolint: cyclocomp_linter
                                standardQuant = FALSE, binarQual = FALSE,
                                storageMode = "vector", names = NULL,
                                p = 2L)
{

  if (is.null(distances) && is.null(d))
  {
    stop("At least the distance matrix or the distance function must be provided") # nolint: line_length_linter
  }

  if (!is.null(d) && !is.data.frame(data))
  {
    stop("Data must be provided and must be a dataframe when 'd' is a distance function") # nolint: line_length_linter
  }

  if (!is.null(distances) && !is.null(data) &&
      nrow(distances) != nrow(data))
  {
    stop("The number of data points must be the same for `distances` and `data`") # nolint: line_length_linter
  }


  if (!is.null(distances))
  {
    dim <- nrow(distances)
    if (inherits(distances, "AbstractSymMat"))
      distances <- distances

    else if (storageMode == "vector")
      distances <- as(distances, "SymVMat")

    else if (storageMode == "matrix")
      distances <- as(distances, "SymMMat")

    else
    {
      stop("Unrecognized storage type")
    }

    if (is.null(names) && !is.null(names(distances)))
      names <- names(distances)

  }
  else
  {
    dim <- nrow(data)
    if (is.null(names))
      names <- seq_len(dim)

    if (storageMode == "vector")
    {
      values <- rep(NA_real_, as.integer((dim - 1L) * dim / 2L))
      distances <- symVMat$new(values = values, dim = dim,
                               names = names, defaultDiag = 0.0)
    }

    else if (storageMode == "matrix")
    {
      values <- matrix(NA_real_, nrow = dim, ncol = dim)
      diag(values) <- 0.0
      distances <- as(values, "SymMMat")
    }
  }

  if (!is.null(d))
  {
    INCORRECTDSTR <- gettext("The provided distance function `d` is incorrect")
    if (testString(d))
    {
      d <- tolower(d)
      if (d == "euclidean")
        d <- euclidean_distance
      else if (d == "manhattan")
        d <- manhattan_distance
      else if (d == "minkowski")
        d <- function_distance_minkowski(p)
      else
      {
        stop(INCORRECTDSTR)
      }

    }
    else if (!is.function(d))
    {
      stop(INCORRECTDSTR)
    }

    if (nrow(data) == 0L || ncol(data) == 0L)
    {
      stop("`data` must have rows and columns")
    }


    if (standardQuant || binarQual)
      data <- normalize_df(data, standardQuant, binarQual)
  }

  distances$defaultDiag <- 0.0
  return(distMat$new(distances = distances, d = d, data = data))
}


#' Extract a sub-matrix in `DistMat` type from a `DistMat` object
#' based on given indices
#'
#' @param x A `DistMat` object
#' @param indices A numeric vector of row/column indices for the sub-matrix.
#' @returns A new `DistMat` object representing the sub-matrix.
#' @noRd
#' @importFrom checkmate assertIntegerish
sousMatrice_distMat <- function(x, indices)
{
  assertIntegerish(indices, lower = 1L, upper = nrow(x))

  distances <- sousMatCarree(x$distances, indices)

  d <- NULL
  data <- NULL

  if (!is.null(x$d))
  {
    d <- x$d
    data <- x$data[indices, ]
  }

  distMat$new(distances = distances, d = d, data = data)

}
