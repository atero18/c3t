#' @include medoids.R

min_distance <- function(distances, nearestNeighbor = NULL)
{
  if (!is.null(nearestNeighbor))
  {
    usefulClusters <- which(!is.nan(nearestNeighbor))
    indexes <- cbind(usefulClusters, nearestNeighbor[usefulClusters])
    values <- distances[indexes]
    return(unname(indexes[which.min(values), , drop = TRUE]))
  }

  diag(distances) <- NaN
  min <- min(distances, na.rm = TRUE)
  which(distances == min, arr.ind = TRUE)[1L, ]
}

# Distances between elements ----------------------------------------------


#' @importFrom rlang arg_match
#' @importFrom checkmate assertNumber
dist_avec_standard <- function(data,
                               distance = c("euclidean", "maximum",
                                            "manhattan", "canberra",
                                            "binary", " minkowski"), p = 2L,
                               standardQuant = FALSE, binarQual = FALSE)
{
  arg_match(distance)
  if (distance == "minkowski")
    assertNumber(p, lower = 1.0)

  if (standardQuant || binarQual)
    data <- normalize_df(data, standardQuant, binarQual)


  return(dist(data, method = distance, p = p))
}


#' Calculate the distance between two elements using the Minkowski distance.
#'
#' @param x First element (numeric vector of non-zero length).
#' @param y Second element. Must be of the same length as x
#' (numeric vector of non-zero length).
#' @param p Real number greater than or equal to one.
#'
#' @returns A positive real number representing the Minkowski
#' distance between x and y.
#' @noRd
distance_minkowski <- function(x, y, p = 2L)
{
  if (length(x) != length(y))
  {
    stop("Both vectors must have the same size")
  }


  if (p == 2.0)
    return(sqrt(sum(abs((x - y))^p)))

  else if (is.infinite(p))
    return(max(abs(x - y)))


  return(sum(abs((x - y))^p)^(1.0 / p))
}

#' @importFrom checkmate assertNumber
function_distance_minkowski <- function(p = 2L)
{
  assertNumber(p, lower = 1.0, finite = FALSE)
  distanceM <- function(x, y) distance_minkowski(x, y, p = p)
  class(distanceM) <- c(class(distanceM), "distCompatible")

  if (is.infinite(p))
    methode <- "maximum"

  else
  {
    methode <- switch(p,
                      "1" = "manhattan",
                      "2" = "euclidean",
                      "Inf" = "maximum",
                      "minkowski")
  }


  assign("dataDist", list(method = methode, p = p),
         envir = environment(distanceM))
  distanceM
}


#' @describeIn distance_minkowski Euclidean distance (p = 2)
#' @noRd
euclidean_distance <-
  distance_euclidienne <-
  function_distance_minkowski(p = 2L)

#' @describeIn distance_minkowski Manhattan distance (p = 1)
#' @noRd
manhattan_distance <-
  distance_manhattan <-
  function_distance_minkowski(p = 1L)


#' @describeIn distance_minkowski Supremum distance (p = + infinite)
#' @noRd
sup_distance <-
  distance_sup <-
  function_distance_minkowski(p = Inf)


#' Calculate the diameter of a set (or subset) of elements.
#'
#' Those functions calculates the diameter (i.e. the maximum distance between
#' two elements of the set) of a set (or subset) of elements
#' based on a distance matrix.
#' @param distances A distance matrix for the set of elements.
#' @name diameter_set
#' @returns The diameter of the set (or subset) as a positive real number.
NULL

#' @describeIn diameter_set get the diameter of an entire set.
#' @export
diameter_set <- function(distances)
{
  diameter_cluster(distances, seq_len(nrow(distances)))
}

#' @describeIn diameter_set calculates the diameter of a subset
#' defined by `indices`.
#'
#' @param indices Indices or names of the elements that
#' are members of the cluster.
#' @export
diameter_cluster <- function(distances, indices)
{
  max(distances[indices, indices])
}

#' @importFrom stats dist
#' @importFrom rlang inject
#' @importFrom checkmate assertNumeric
calculate_distances <- function(indexs, d, data)
{
  # If `d` is compatible with [stats::dist] it is preferable to use it
  # instead.
  if (inherits(d, "distCompatible"))
  {
    elementsConcernes <- unique(c(indexs[, 1L], indexs[, 2L]))
    parametresDist <- get("dataDist", envir = environment(d))
    resDist <- inject(dist(x = data[elementsConcernes, ], !!!parametresDist))

    posElementsDist <- attr(resDist, "Label")
    if (is.null(posElementsDist))
      posElementsDist <- seq_len(nrow(data))

    else
      posElementsDist <- as.integer(posElementsDist)

    indexsBis <- indexs
    indexsBis[, 1L] <- match(indexs[, 1L], posElementsDist)
    indexsBis[, 2L] <- match(indexs[, 2L], posElementsDist)
    resDist <- as.matrix(resDist)[indexsBis]

    return(resDist)
  }
  else
  {
    res <- c3t_apply(indexs, 1L, function(couple) d(data[couple[1L], ],
                                                    data[couple[2L], ]))

    assertNumeric(res, lower = 0.0,
                  any.missing = FALSE, all.missing = FALSE,
                  null.ok = FALSE)

    return(res)
  }

}
