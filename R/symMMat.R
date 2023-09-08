#' @include abstractSymMat.R

# Création de la classe

#' @keywords internal
symMMat <- setRefClass("SymMMat",
                       fields = list(defaultDiag = "numericOrLogicalOrNULL",
                                     values = "matrix"),
                       contains = "AbstractSymMat")


symMMat$methods(
  which_na = function(subsetRow = seq_len(nrow(.self)),
                      subsetCol = seq_len(nrow(.self)),
                      nanIsNA = TRUE,
                      simplify = TRUE)
  {
    na <- is.na(values[subsetRow, subsetCol, drop = FALSE])

    if (!nanIsNA)
      na <- na & !is.nan(values[subsetRow, subsetCol, drop = FALSE])

    if (!any(na))
      return(matrix(nrow = 0L, ncol = 2L))

    na <- which(na, arr.ind = TRUE)
    colnames(na) <- NULL

    na[, 1L] <- subsetRow[na[, 1L]]
    na[, 2L] <- subsetCol[na[, 2L]]

    if (simplify)
    {
      maskSuperior <- na[, 1L] > na[, 2L]
      if (any(maskSuperior))
      {
        na[maskSuperior, 2L:1L] <- na[maskSuperior, 1L:2L, drop = FALSE]
        if (sum(maskSuperior) < nrow(na))
          na <- unique(na, MARGIN = 1L)
      }
    }

    na
  }
)

#' @rdname abstractSymMat_access
#' @keywords internal
setMethod(
  "[",
  signature(x = "SymMMat", i = "matrix", j = "missing", drop = "ANY"),
  function(x, i, j, drop)
  {
    if (nrow(i) == 0L)
      return(NULL)

    else if (ncol(i) != 2L)
    {
      stop("`i` must have exactly two columns")
    }

    return(x$values[i])
  }
)

#' @rdname abstractSymMat_access
#' @keywords internal
setMethod(
  "[",
  signature(x = "SymMMat", i = "numeric", j = "numeric", drop = "ANY"),
  function(x, i, j, drop)
  {
    drop <- ifelse(missing(drop), TRUE, isTRUE(drop))

    x$values[i, j, drop = drop]
  }
)


# Assignation des données

assignation_symMat <- function(x, indexs, values)
{
  if (!is.null(x$defaultDiag))
  {
    masqueEgalite <- any(indexs[, 1L] == indexs[, 2L])
    if (any(masqueEgalite) && any(values[masqueEgalite] != x$defaultDiag))
    {
      stop("Default value is already defined for the diagonal")
    }
  }


  x$values[indexs] <- values
  x$values[indexs[, 2L:1L, drop = FALSE]] <- values

  return(x)
}

# nocov start
#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod(
  "[",
  signature(x = "SymMMat", i = "numeric", j = "numeric", value = "numeric"),
  function(x, i, j, value) {
    if (any(i > x$dim) || any(j > x$dim))
    {
      stop("Specified indexes are incorrect")
    }


    x$values[i, j] <- x$values[j, i] <- value

    return(x)
  }
)

#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod(
  "[",
  signature(x = "SymMMat", i = "numeric", j = "missing", value = "numeric"),
  function(x, i, j, value)
  {
    x$values[i, ] <- value
    x$values[, i] <- value
    return(x)
  }
)

#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod(
  "[",
  signature(x = "SymMMat", i = "missing", j = "numeric", value = "numeric"),
  function(x, i, j, value) x[j, ] <- value
)

#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod(
  "[",
  signature(x = "SymMMat", i = "matrix", j = "missing", value = "numeric"),
  function(x, i, j, value) assignation_symMat(x, i, value)
)

#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod(
  "[",
  signature(x = "SymMMat", i = "missing", j = "missing", value = "matrix"),
  function(x, i, j, value)
  {
    if (!isSymmetric(value))
    {
      stop("Matrix must be symmetric")
    }

    if (nrow(value) != nrow(x))
    {
      stop("Matrices must have the same size")
    }

    if (!is.null(x$defaultDiag) &&
        any(diag(value) != x$defaultDiag))
    {
      stop("Defaut value for the diagonal isn't respected")
    }

    x$values <- value

    return(x)
  }
)

# nocov end

# Conversions
setAs(
  "matrix", "SymMMat",
  function(from)
  {
    names <- names(from)
    if (is.null(names))
      names <- seq_len(nrow(from))
    res <- symMMat$new(values = from, dim = nrow(from),
                       names = names, defaultDiag = NULL)
    return(res)
  }
)

setAs(
  "dist", "SymMMat",
  function(from)
  {
    res <- as(as.matrix(dist), "SymMMat")
    res$defaultDiag <- 0.0
    res
  }
)

setMethod(
  "sousMatCarree",
  signature(x = "SymMMat", indexs = "numeric"),
  function(x, indexs)
  {
    data <- x[indexs, indexs]
    res <- as(data, "SymMMat")
    res$defaultDiag <- x$defaultDiag
    res$names <- x$names[indexs]
    res
  }
)
