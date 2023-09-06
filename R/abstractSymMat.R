#' @include arguments.R
#'
#'
# -- `dist` object from `stats` package
#' @seealso [stats::dist()]
#' @noRd
#' @importFrom stats dist
setOldClass("dist", where = environment())

#' @title Abstract Class: `AbstractSymMat`
#' @description A class handling the storage of a symmetric square matrix.
#'
#' @field dim The dimension of the square matrix (numeric).
#' @field names A vector storing the names for rows and columns (vector).
#' @field defaultDiag The default values for diagonal elements
#' (numericOrLogicalOrNULL).
#' @keywords internal
#' @name abstractSymMat
setRefClass(
  "AbstractSymMat",
  fields = list(dim = "integer", names = "vector",
                defaultDiag = "numericOrLogicalOrNULL"),
  methods = list(
    anyNA = function()
    {
      if (exists("anyNA_bool", envir = environment(anyNA)) &&
          !anyNA_bool)
      {
        return(FALSE)
      } else
      {
        anyNA_bool <- base::anyNA(values) # nolint
        if (!anyNA_bool) {
          assign("anyNA_bool", FALSE, envir = environment(anyNA))
        }
        return(anyNA_bool)
      }
    },
    which_na = function(subsetRow = seq_len(nrow(.self)),
                        subsetCol = seq_len(nrow(.self)),
                        nanIsNA = TRUE,
                        simplify = TRUE) NULL
  )
)

#' Access properties of a `AbstractSymMat` object or subclasses
#' @param x An `AbstractSymMat` object.
#' @name abstractSymMat_properties
#' @rdname abstractSymMat_properties
#' @keywords internal
NULL

# @describeIn abstractSymMat_properties names of the elements
#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("names", signature(x = "AbstractSymMat"), function(x) x$names)

# @describeIn abstractSymMat_properties Number of elements
#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("length", signature(x = "AbstractSymMat"), function(x) x$dim^2L)

# @describeIn abstractSymMat_properties Dimensions of the object (2 positive
#' integers vector)
#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("dim", signature(x = "AbstractSymMat"), function(x) c(x$dim, x$dim))

#' @importFrom methods setGeneric
setGeneric("nrow")

# @describeIn abstractSymMat_properties Number of rows
#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("nrow", signature(x = "AbstractSymMat"), function(x) x$dim)

#' @importFrom methods setGeneric
setGeneric("ncol")

# @describeIn abstractSymMat_properties Number of columns
#' @rdname abstractSymMat_properties
#' @keywords internal
setMethod("ncol", signature(x = "AbstractSymMat"), function(x) x$dim)

#' Replace the names of elements in an `AbstractSymMat` or subclasses
#' @param x An `AbstractSymMat` (or subclasses) object.
#' @param value A vector of names to be assigned to rows and columns.
#' @returns The updated object.
#' @name abstractSymMat_replace_name
#' @rdname abstractSymMat_replace_name
#' @keywords internal
#' @importFrom methods setReplaceMethod
NULL

#' @rdname abstractSymMat_replace_name
#' @keywords internal
setReplaceMethod(
  "names",
  signature(x = "AbstractSymMat", value = "vector"),
  function(x, value)
  {

    if (length(value) != x$dim)
    {
      stop("The name vector must be of size `dim`")
    }


    if (anyDuplicated(value) == 0L)
      return("The given names must be unique")

    x$names <- value
    return(x)
  }
)


#' Inherents properties of an `AbstractSymMat
#' ` object
#' @name abstractSymMat
#' @rdname abstractSymMat
#' @param x,object a `AbrastSymMat` object
#' @keywords internal
NULL

#' @rdname abstractSymMat_properties
#' @returns For `isSymmetric` : always TRUE as any `AbstractSymMat`
#' is symmetric.
#' @keywords internal
setMethod(
  "isSymmetric",
  signature(object = "AbstractSymMat"),
  function(object, ...) TRUE
)

#' @rdname abstractSymMat_properties
#' @returns For `is.matrix` : always `TRUE` as any `AbstractSymMat` is a matrix.
#' @keywords internal
setMethod("is.matrix", signature(x = "AbstractSymMat"), function(x) TRUE)

#' Is there any NA is our symmetric matrix?
#' @param x An `AbstractSymMat` object.
#' @returns Logical value indicating if any NA values are present in the
#' `AbstractSymMat` object.
#' @keywords internal
setMethod("anyNA", signature(x = "AbstractSymMat"), function(x) x$anyNA())

#' @importFrom methods setGeneric
setGeneric("diag")

setMethod(
  "diag", # nolint: indentation_linter
  signature = "AbstractSymMat",
  function(x) x[cbind(seq_len(x$dim), seq_len(x$dim))]
)

#' Where are located NA is our symmetric matrix?
#' @param x An `AbstractSymMat` object.
#' @returns A logical matrix of the same dimensions as the `AbstractSymMat`
#' object, indicating if each element is NA.
#' @keywords internal
setMethod("is.na", signature(x = "AbstractSymMat"), function(x) is.na(x$values))

# Data access abstract rules
#' @importFrom methods selectMethod
accesMatriceR <-
  selectMethod(
    "[",
    signature(x = "matrix", i = "numeric", j = "numeric", drop = "logical")
  )
#' Access to `AbstractSymMat` data`
#'
#' Rules are the following:
#' * when rows and columns are given a submatrix is returned with the
#' asked rows and columns
#' * when only rows (resp. columns are given) a submatrix is return
#' with the asked rows (resp. columns) and the entire columns (resp. rows)
#' * When rows and columns are missing the entire matrix is sent back.
#' Depending on `drop` (a flag) is there is only one row and / or one
#' column given the result will be simplified in a vector form or not.
#' @name abstractSymMat_access
#' @rdname abstractSymMat_access
#' @param i A vector composed of the indices of the wished rows.
#' @param j A vector composed of the indices of the wished columns
#' @param x An `AbstractSymMat` or inheriting object
#' @returns Depending on `drop` and the dimensions of `i` and `j`, a matrix or
#' a vector.
#' @keywords internal
NULL

#' @rdname abstractSymMat_access
#' @keywords internal
setMethod(
  "[",
  signature(x = "AbstractSymMat", i = "numeric", j = "missing", drop = "ANY"),
  function(x, i, j, drop) x[i, seq_len(x$dim)]
)

#' @rdname abstractSymMat_access
#' @keywords internal
setMethod(
  "[",
  signature(x = "AbstractSymMat", i = "missing", j = "numeric", drop = "ANY"),
  function(x, i, j, drop) (x[seq_len(x$dim), j])
)

#' @rdname abstractSymMat_access
#' @keywords internal
setMethod(
  "[",
  signature(x = "AbstractSymMat", i = "missing", j = "missing", drop = "ANY"),
  function(x, i, j, drop) (x[seq_len(x$dim), seq_len(x$dim)])
)

# Equality properties
#' @importFrom methods validObject
setMethod(
  "all.equal",
  signature(target = "AbstractSymMat", current = "matrix"),
  function(target, current, ...)
  {
    validObject(target)
    all(target[] == current)
  }
)

#' @importFrom methods validObject
setMethod(
  "all.equal",
  signature(target = "matrix", current = "AbstractSymMat"),
  function(target, current, ...)
  {
    validObject(current)
    all.equal(current, target)
  }
)

setMethod(
  "all.equal",
  signature(target = "AbstractSymMat", current = "AbstractSymMat"),
  function(target, current, ...) all(current[] == target[])
)

#' Replace value of a `AbstractSymMat` object or subclasses
#'
#' Because an `AbstractSymMat` is symmetrical replacing values at some position
#' (i,j) will automatically set the same value at position (j,i).
#' @param x a `AbstractSymMat` object or subclasses
#' @param i,j strictly positive integers vectors giving
#' the position of the elements (rows for `i` and columns for `j`)
#' that will be replaced. If `i` is a matrix then it is supposed that
#' one row corresponds to one element (first column being the row and second
#' column being the column). 0, 1 or both of them can be missing.
#' @param values A vector of values to replace the elements of `x`. If
#' `i` and `j` are missing can be a matrix.
#' @returns the updated object `x`, even if modifications are made in place.
#' @name abstractSymMat_replace
#' @rdname abstractSymMat_replace
#' @keywords internal
NULL

#' Replace diagonal values of a `AbstractSymMat` object or subclasses
#' @param x a `AbstractSymMat` object or subclasses
#' @param value a 1 or `nrow(x)` length vector
#' @name abstractSymMat_replace_diag
#' @rdname abstractSymMat_replace_diag
#' @keywords internal
NULL

#' @param i A numeric vector representing row indices.
#' @param j A numeric vector representing column indices.
#' @returns A matrix containing all combinations of elements from
#' `i` and `j`.
#' @noRd
vecs_pos_vers_indexs_pos <- function(i, j)
{
  cbind(rep(i, times = length(j)),
        rep(j, each = length(i)))
}

#' @param data A vector of data to be converted to a matrix.
#' @param rowNames A vector of row names (optional).
#' @param colNames A vector of column names (optional).
#' @param drop Logical value indicating whether to drop dimensions of size 1.
#' @returns A matrix with the given data and optional row/column names.
#' @noRd
#' @importFrom checkmate testScalar
data_indexs_pos_vers_matrice <-
  function(data, rowNames = NULL,
           colNames = NULL, drop = FALSE)
  {
    if (drop && (testScalar(rowNames) || testScalar(colNames)))
      return(data)


    res <- matrix(data, byrow = FALSE,
                  nrow = length(rowNames),
                  ncol = length(colNames))
    rownames(res) <- rowNames
    colnames(res) <- colNames

    res
  }


#' @title Extract Submatrix from `AbstractSymMat`
#' @param x An `AbstractSymMat` object.
#' @param indexs A numeric vector representing row and column indices.
#' @returns A submatrix of `AbstractSymMat` based on the given row and
#' column indices.
#' @noRd
sousMatCarree <- function(x, indexs) x[indexs, indexs]

#' @importFrom methods setGeneric
setGeneric("sousMatCarree")

#' @importFrom stats as.dist
setMethod(
  "as.dist", # nolint: indentation_linter
  signature(m = "AbstractSymMat", diag = "ANY", upper = "ANY"),
  function(m, diag, upper) as(m, "dist"))

#' @importFrom stats as.dist
#' @importFrom methods setAs
setAs("AbstractSymMat", "dist", function(from) as.dist(from[]))
