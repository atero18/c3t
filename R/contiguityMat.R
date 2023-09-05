#' @include symSparseMat.R

# Create the class

#' `ContiguityMat`
#' A class representing a contiguity matrix.
#' @slot dim The size of the set.
#' @slot values A matrix of contiguity values
#' @slot names Names for the rows.
#' @slot defaut Logical value indicating if the default value
#'     should be considered for non-contiguous elements.
#' @slot defaultDiag Logical value indicating if
#'     the diagonal default value should be considered.
#' @keywords internal
contiguityMat <- setRefClass("ContiguityMat", contains = "SymSparseMat")

# Data retrieval

#' @rdname abstractSymMat_access
#' @keywords internal
#' @importFrom methods callNextMethod
setMethod(
  "[",
  signature(x = "ContiguityMat", i = "numeric", j = "numeric", drop = "ANY"),
  function(x, i, j, drop)
  {
    res <- callNextMethod()
    return(res == 1.0 | res)
  }
)

#' @rdname abstractSymMat_access
#' @keywords internal
#' @importFrom methods callNextMethod
setMethod(
  "[",
  signature(x = "ContiguityMat", i = "missing", j = "numeric", drop = "ANY"),
  function(x, i, j, drop)
  {
    res <- callNextMethod()
    return(res == 1.0 | res)
  }
)

#' @rdname abstractSymMat_access
#' @keywords internal
#' @importFrom methods callNextMethod
setMethod(
  "[",
  signature(x = "ContiguityMat", i = "numeric", j = "missing", drop = "ANY"),
  function(x, i, j, drop)
  {
    res <- callNextMethod()
    return(res == 1.0 | res)
  }
)

#' @rdname abstractSymMat_access
#' @keywords internal
#' @importFrom methods callNextMethod
setMethod(
  "[",
  signature(x = "ContiguityMat", i = "missing", j = "missing", drop = "ANY"),
  function(x, i, j, drop)
  {
    res <- callNextMethod()
    return(res == 1.0 | res)
  }
)

# Overloading igraph functions

#' @importFrom methods setGeneric
#' @importFrom igraph are_adjacent
setGeneric("are_adjacent", igraph::are_adjacent)

#' @importFrom igraph are_adjacent
setMethod(
  "are_adjacent",
  signature(graph = "ContiguityMat", v1 = "numeric", v2 = "numeric"),
  function(graph, v1, v2)
  {
    assertCount(v1, positive = TRUE)
    assertCount(v2, positive = TRUE)

    ifelse(v1 == v2, TRUE, graph[v1, v2])
  }
)

#' ContiguityMat
#'
#' Constructor for the ContiguityMat class.
#' @param M A contiguity matrix or a matrix of contiguity lists
#'          (2-column matrix consisting of elements being indexes).
#' @param dim The size of the set. Ignored if M is a contiguity matrix
#'            (strictly positive integer).
#'
#' @returns A contiguity matrix of type ContiguityMat
#' @examples
#' # Create a contiguity matrix
#' M <- matrix(c(F, T, T, F), ncol = 2)
#' dim <- 2
#' cont_mat <- ContiguityMat(M, dim)
#' @noRd
#' @importFrom checkmate testIntegerish assertCount
contiguityMat <- function(M, dim = NULL)
{
  # Checking arguments
  if (testContiguityMatrix(M))
  {
    couples <- expand.grid(seq_len(nrow(M)), seq_len(nrow(M)))
    couples <- couples[couples[, 1L] < couples[, 2L], ]
    res <- apply(couples, 1L, function(x) M[x[1L], x[2L]])
    values <- cbind(couples, res)
    values <- as.matrix(values[values[, 3L], ])
    return(new("ContiguityMat", dim = nrow(M), values = values,
               names = seq_len(nrow(M)), defaut = FALSE, defaultDiag = FALSE))
  }
  else if (is.matrix(M) && ncol(M) > 2L &&
           testIntegerish(M[, 1L], lower = 1L,
                          any.missing = FALSE, all.missing = FALSE) &&
           testIntegerish(M[, 2L], lower = 1L,
                          any.missing = FALSE, all.missing = FALSE))
  {
    assertCount(dim, positive = TRUE)

    if (any(M[, 1L:2L] > dim))
    {
      stop("At least one datum indicates a rank greater than the number of individuals indicated via `dim`") # nolint: line_length_linter
    }

    values <- cbind(M[, 1L:2L], TRUE)

    return(new("ContiguityMat", dim = dim, values = values,
               names = seq_len(nrow(M)), defaut = FALSE, defaultDiag = FALSE))
  }
  else
  {
    stop("`M` is not in the correct format")
  }



}

# Conversion from a contiguity matrix to a contiguity graph (as)

#' @importFrom methods setOldClass
setOldClass("igraph")

#' @importFrom methods setAs
setAs("ContiguityMat", "igraph",
      function(from) contiguity_matrix_to_graph(from[]))

# Conversion from a contiguity graph to a contiguity matrix
setAs("igraph", "ContiguityMat",
      function(from) graph_to_contiguity_matrix(from))
