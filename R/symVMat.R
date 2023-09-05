#' @include abstractSymMat.R

#' @keywords internal
symVMat <- setRefClass(
  "SymVMat",
  fields = list(defaultDiag = "numericOrLogicalOrNULL", values = "ANY"),
  contains = "AbstractSymMat"
)

#' @importFrom checkmate assert checkNumber checkFlag
symVMat$methods(
  setDefautDiag = function(defautD)
  {
    if (length(defautD) > 1L)
      defautD <- defautD[1L]

    assert(checkNumber(defautD, null.ok = TRUE), checkFlag(defautD))
    defaultDiag <<- defautD # nolint: undesirable_operator_linter
  }
)
position_element_symVMat <- function(x, i, j = NULL)
{
  if (is.null(j))
    indexs <- i
  else
    indexs <- matrix(c(i, j), ncol = 2L)


  if (!is.null(x$defaultDiag) && any(indexs[, 1L] == indexs[, 2L]))
  {
    stop("Element is not contained in values because is part of the diagonal")
  }

  masqueEchange <- indexs[, 1L] > indexs[, 2L]
  if (any(masqueEchange))
    indexs[masqueEchange, c(1L, 2L)] <- indexs[masqueEchange, c(2L, 1L)]

  ##indexs = cbind(pmin.int(indexs[, 1L], indexs[, 2L]),
  ##               pmax.int(indexs[, 1L], indexs[, 2L]))

  pos <- integer(nrow(indexs))

  # Moving with the row
  masqueLigne <- indexs[, 1L] >= 2L
  if (any(masqueLigne))
    pos[masqueLigne] <-
    (indexs[masqueLigne, 1L] - 1L) *
    (x$dim - (indexs[masqueLigne, 1L] - 2L) / 2L - !is.null(x$defaultDiag))

  # Moving with the column
  pos <- pos + (indexs[, 2L] - indexs[, 1L]) + is.null(x$defaultDiag)

  return(pos)
}

#' @importFrom utils tail
#' @importFrom checkmate assertIntegerish
indexs_elements_symVMat <- function(x, pos)
{
  assertIntegerish(pos, lower = 1L, upper = length(x@values))

  indexs <- matrix(nrow = length(pos), ncol = 2L)
  posInitLigne <- 1L +
    c(0L, cumsum(x$dim - as.integer(!is.null(x$defaultDiag)):1L))
  indexs[, 1L] <-  vapply(pos,
                          function(p) tail(which(posInitLigne <= p), n = 1L),
                          integer(1L))

  indexs[, 2L] <-  pos - posInitLigne[indexs[, 1L]] +
    1L + !is.null(x$defaultDiag)

  return(indexs)

}

#' @rdname abstractSymMat_access
#' @keywords internal
#' @importFrom checkmate assertIntegerish
setMethod(
  "[",
  signature(x = "SymVMat", i = "matrix", j = "missing", drop = "ANY"),
  function(x, i, j, drop)
  {
    if (ncol(i) != 2L)
    {
      stop("`i` must exactly have 2 columns")
    }

    if (nrow(i) == 0L)
      return(NULL)


    assertIntegerish(i[, 1L],
                     lower = 1L, upper = nrow(x),
                     any.missing = FALSE, all.missing = FALSE)

    assertIntegerish(i[, 2L],
                     lower = 1L, upper = ncol(x),
                     any.missing = FALSE, all.missing = FALSE)

    grille <- i
    donnees <- rep(NA, nrow(grille))


    if (!is.null(x$defaultDiag))
    {
      masqueEgalite <- grille[, 1L] == grille[, 2L]
      masque <- !masqueEgalite
      if (any(masqueEgalite))
      {
        donnees[masqueEgalite] <- x$defaultDiag
        grille <- grille[masque, ]
        if (sum(masque) == 1L)
          grille <- matrix(grille, ncol = 2L)
      }
    }

    else
      masque <- rep(TRUE, nrow(grille))

    if (any(masque))
    {
      positions <- position_element_symVMat(x, grille)
      donnees[masque] <- x$values[positions]
    }

    return(donnees)

  }
)


#' @rdname abstractSymMat_access
#' @keywords internal
#' @importFrom checkmate testScalar
setMethod(
  "[",
  signature(x = "SymVMat", i = "numeric", j = "numeric", drop = "ANY"),
  function(x, i, j, drop)
  {
    if (missing(drop))
      drop <- TRUE

    else
      drop <- isTRUE(drop)

    grille <- vecs_pos_vers_indexs_pos(i, j)

    if (!drop || !testScalar(i) && !testScalar(j))
      return(data_indexs_pos_vers_matrice(x[grille, ],
                                          names(x)[i],
                                          names(x)[j]))

    else
      return(x[grille, ])
  }
)

# Data affectation
#' @importFrom checkmate assertNumber
assignationSimple_symVMat <- function(x, i, j, value)
{
  if (!is.null(x$defaultDiag) && i == j)
  {
    if (value != x$defaultDiag)
    {
      stop("A default value is defined for the diagonal")
    }
    else
      return(x)
  }


  assertNumber(value)

  pos <- position_element_symVMat(x, i, j)

  x$values[pos] <- value


  return(x)
}

#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod(
  "[",
  signature(x = "SymVMat", i = "numeric", j = "numeric", value = "numeric"),
  function(x, i, j, value) {
    if (any(i > x$dim) || any(j > x$dim))
    {
      stop("Specified indices are incorrect")
    }

    grille <- expand.grid(i, j)
    apply(grille, 1L,
          function(couple) assignationSimple_symVMat(x, couple[1L],
                                                     couple[2L], value))

    return(x)
  }
)

#' @rdname abstractSymMat_replace
#' @keywords internal
#' @importFrom checkmate assertIntegerish
setReplaceMethod(
  "[",
  signature(x = "SymVMat", i = "matrix", j = "missing", value = "numeric"),
  function(x, i, j, value)
  {

    assertIntegerish(i[, 1L],
                     lower = 1L, upper = nrow(x),
                     any.missing = FALSE, all.missing = FALSE)

    assertIntegerish(i[, 2L],
                     lower = 1L, upper = ncol(x),
                     any.missing = FALSE, all.missing = FALSE)

    if (length(value) != nrow(i))
    {
      stop("The size of `value` doesn't correspond to the number of indexes")
    }


    lapply(seq_len(nrow(i)),
           function(k) assignationSimple_symVMat(x, i[k, 1L],
                                                 i[k, 2L], value[k]))

    return(x)
  }
)


#' @rdname abstractSymMat_replace_diag
#' @keywords internal
setReplaceMethod(
  "diag",
  signature(x = "SymVMat", value = "ANY"),
  function(x, value)
  {
    x$setDefautDiag(value)
    return(x)
  }
)

#' @importFrom methods setAs
setAs("dist", "SymVMat", function(from)
{
  dim <- as.integer((1L + sqrt(1L + 8L * length(from))) / 2L)
  names <- names(from)
  if (is.null(names))
    names <- seq_len(dim)

  symVMat$new(values = from, dim = dim, names = names, defaultDiag = 0.0)
})

#' @importFrom methods as setAs
#' @importFrom stats as.dist
setAs("matrix", "SymVMat", function(from)
{
  res <- as(as.dist(from), "SymVMat")
  res$defaultDiag <- from[1L, 1L]
  res
})

#' @importFrom methods as
setMethod(
  "sousMatCarree",
  signature(x = "SymVMat", indexs = "numeric"),
  function(x, indexs)
  {
    data <- x[indexs, indexs]
    res <- as(data, "SymVMat")
    res@names <- x$names[indexs]
  }
)

#' @importFrom methods setAs
#' @importFrom stats dist as.dist
setAs(
  "SymVMat", "dist", # nolint: indentation_linter
  function(from)
  {
    if (inherits(from$values, "dist"))
      return(from$values)
    else
      return(as.dist(from[]))
  }
)
