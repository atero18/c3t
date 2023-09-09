#' @include pbCon.R

#' Converts contiguity between `igraph` and matrix
#' @param x object to be converted.
#' @name contiguity_conversion
#' @rdname contiguity_conversion
#' @keywords internal
#' @importFrom igraph graph_from_adjacency_matrix
contiguity_matrix_to_graph <- function(x)
{
  # Argument verification
  assertContiguityMatrix(x, isComplete = TRUE)

  graph <- graph_from_adjacency_matrix(x,
                                       mode = "undirected",
                                       diag = FALSE)

  return(graph)
}


#' @rdname contiguity_conversion
#' @keywords internal
graph_to_contiguity_matrix <- function(x)
{
  # If graph is a `pbCon`, use its contiguity graph
  if (is_pbCon(x))
  {
    x$contiguityMatrix <- as.matrix(x$contiguity[]) == 1.0
    diag(x$contiguityMatrix) <- TRUE
  }
  else
  {
    assertContiguityGraph(x)
    matrix <- as.matrix(x[]) == 1.0
    diag(matrix) <- TRUE
    return(matrix)
  }
}

#' Object representing complete contiguity
#'
#' There is complete contiguity when all elements in a set are directly
#' connected to the others.
#' @param x argument giving information about the number of elements.
#' Can be an strictly positive integer, a `pbCon` object or a matrix.
#' @name complete_contiguity
#' @rdname complete_contiguity
#' @keywords internal
#' @importFrom checkmate testCount
complete_contiguity_matrix <- function(x)
{
  if (testCount(x, positive = TRUE))
    n <- x

  else if (is_pbCon(x))
    n <- nrow(x)

  else if (is.matrix(x))
    n <- nrow(x)

  else
  {
    stop("`x` is not of the correct type")
  }

  mat <- matrix(TRUE, nrow = n, ncol = n)
  diag(mat) <- FALSE
  mat
}

#' @rdname complete_contiguity
#' @keywords internal
#' @seealso [igraph::make_full_graph()]
#' @importFrom igraph make_full_graph
#' @importFrom checkmate testCount
complete_contiguity_graph <- function(x)
{
  if (testCount(x, positive = TRUE))
    n <- x

  else if (is.matrix(x))
    n <- nrow(x)

  else
  {
    stop("`x` is not of the correct type")
  }

  make_full_graph(n, directed = FALSE, loops = FALSE)
}

#' @importFrom igraph is_connected
setGeneric("is_connected", igraph::is_connected)

#' Check Connectivity
#'
#' Verifies via the contiguity matrix whether a set is connected.
#' @param graph Contiguity matrix of the elements.
#' (contiguity matrix)
#' @returns TRUE if the set is connected, FALSE otherwise.
#' @seealso [igraph::is_connected()]
#' @keywords internal
#' @example inst/examples/is_connected.R
#' @export
#' @importFrom igraph is_connected
setMethod(
  "is_connected",
  signature(graph = "matrix", mode = "ANY"),
  function(graph, mode)
  {
    if (nrow(graph) <= 1L)
      return(TRUE)

    graph <- contiguity_matrix_to_graph(graph)

    is_connected(graph, mode = "weak")
  }
)

#' Is a partition checking constraints?
#' @param partition Partition vector of the set (numeric vector
#' composed of non-null integers).
#' @param contiguity contiguity matrix or graph. If `NULL` then the contiguity
#' constraint is not considered.
#' @returns A flag equals to `TRUE` if `partition` check the requested
#' constraints, `FALSE` otherwise.
#' @name check_solution
NULL

#' @describeIn check_solution  Verifies if all classes in a set are connected.
#' @importFrom igraph induced_subgraph is_connected
#' @export
is_regionalisation <- function(partition, contiguity)
{
  # Argument verification
  assertPartition(partition)
  if (is.matrix(contiguity) || inherits(contiguity, "Matrix"))
    contiguity <- contiguity_matrix_to_graph(contiguity)

  clustersIDs <- unique(partition)

  res <- vapply(clustersIDs, function(i)
  {
    clustersPoints <- which(partition == i)

    # A cluster of size 1 is, by definition, connected
    if (length(clustersPoints) == 1L)
      return(TRUE)

    g <- induced_subgraph(contiguity, clustersPoints)
    return(is_connected(g, mode = "weak"))
  }, logical(1L))

  return(all(res))
}


#' Connected Components
#'
#' Determines the set of connected components in a set via
#' its contiguity matrix.
#' @param contiguity The contiguity matrix of the `n` elements.
#' @returns A integer vector of length `n` with for each element the id of
#' its connected component.
#' @keywords internals
#' @importFrom igraph components
connected_components <- function(contiguity)
{
  graph <- contiguity_matrix_to_graph(contiguity)

  components(graph, mode = "weak")$membership

}

#' @describeIn connected_components Returns the number of
#' connected components in the set.
#' @keywords internal
#' @importFrom igraph count_components
nb_connected_components <- function(contiguity)
{
  graph <- contiguity_matrix_to_graph(contiguity)

  count_components(graph, mode = "weak")
}
