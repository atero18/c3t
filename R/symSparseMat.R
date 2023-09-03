#' @include abstractSymMat.R

# Creation of `SymSparseMat` class
#' @keywords internal
symSparseMat <-
  setRefClass("SymSparseMat",
              fields = list(values = "ANY",
                            defaut = "numericOrLogical",
                            defaultDiag = "numericOrLogicalOrNULL"),
              contains = "AbstractSymMat")

# Function to retrieve a specific element from the distance matrix.
# Checks if the response is trivial (i = j), and if the element is present
# (for the pair (i,j) or (j,i) by symmetry).


#' @importFrom methods setGeneric
#' @importFrom checkmate assertInt
setGeneric("get_elem_symSparseMat", function(x, i, j)
{
  assertInt(i, lower = 1L, upper = nrow(x))
  assertInt(j, lower = 1L, upper = ncol(x))

  .get_elem_symSparseMat(x, i, j)

})

#' @importFrom methods setGeneric
setGeneric(".get_elem_symSparseMat", function(x, i, j)
{

  maskReduction <- x$values[, 1L] %in% c(i, j) & x$values[, 2L] %in% c(i, j)
  reducedList <- matrix(x$values[maskReduction, ], ncol = 3L)


  grille <- as.matrix(expand.grid(i, j), ncol = 2L)
  donnees <- apply(grille, 1L, function(couple) {

    if (couple[1L] == couple[2L] && !is.null(x$defaultDiag))
      return(x$defaultDiag)

    if (couple[1L] > couple[2L])
    {
      t <- couple[1L]
      couple[1L] <- min(t, couple[2L])
      couple[2L] <- max(t, couple[1L])
    }


    if (length(reducedList) > 0L)
    {
      posValeur <- which(reducedList[, 1L] == couple[1L] &
                         reducedList[, 2L] == couple[2L])

      if (length(posValeur) == 0L)
        return(x$defaut)

      else
        return(reducedList[posValeur[1L], 3L])
    }


    else
      return(x$defaut)
  })

  if (length(i) == length(j) && length(i) == 1L)
    return(donnees)

  else
  {
    donnees <- matrix(donnees, byrow = FALSE,
                      nrow = length(i),
                      ncol = length(j))

    rownames(donnees) <- i
    colnames(donnees) <- j
    return(donnees)
  }

})

#' @rdname abstractSymMat_access
#' @keywords internal
#' @importFrom checkmate assertVector assertIntegerish
setMethod("[", signature(x = "SymSparseMat", i = "numeric",
                         j = "numeric", drop = "ANY"),
          function(x, i, j, drop) {
            assertVector(i)
            assertIntegerish(i,
                             lower = 1L, upper = nrow(x),
                             any.missing = FALSE,
                             all.missing = FALSE)

            assertVector(j)
            assertIntegerish(j,
                             lower = 1L, upper = ncol(x),
                             any.missing = FALSE,
                             all.missing = FALSE)

            .get_elem_symSparseMat(x, i, j)
          })


# Data affectation

#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod("[", signature(x = "SymSparseMat", i = "numeric",
                                j = "numeric", value = "numeric"),
                 function(x, i, j, value) {
                   if (i > x$dim || j > x$dim)
                   {
                     stop("Given indexes are incorrect")
                   }

                   if (i == j)
                   {
                     if (!is.null(x$defaultDiag) && value != x$defaultDiag)
                     {
                       stop("For the diagonal a default value has been given")
                     }

                     else
                       return(x)
                   }

                   if (i > j)
                   {
                     t <- i
                     i <- min(t, j)
                     j <- max(t, j)
                   }

                   if (nrow(x$values) > 0L)
                   {
                     mask <- x$values[, 1L] == i & x$values[, 2L] == j

                     if (any(mask))
                       x$values[mask, 3L] <- value

                     return(x)
                   }

                   x$values <- rbind(x$values, c(i, j, value))

                   return(x)

                 })
