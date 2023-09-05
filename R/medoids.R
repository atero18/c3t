#' @title Medoid Calculation
#'
#' @description Determines the medoid of one or several subsets of
#' elements using a distance matrix.
#' @param distances The distance matrix used to find the
#' medoid(s) (distance matrix).
#' @returns For `medoid`: depending on `byName`, the name or the
#' identifier of the medoid.
#' @name medoid
NULL


#' @describeIn medoid Determine the medoid of a set.
#' @param byName `TRUE` if the returned result should be the row name
#' corresponding to the medoid, `FALSE` if it should be the row number.
#' Ignored if `distances` is not named. (flag)
#' @export
medoid <- function(distances, byName = FALSE)
{
  if (length(distances) == 0L)
  {
    stop("the distance matrix is empty")
  }

  else if (nrow(distances) <= 2L)
    posMedoid <- 1L

  else
    posMedoid <- which.min(rowSums(distances))

  ifelse(byName && !is.null(rownames(distances)),
         rownames(distances)[posMedoid],
         posMedoid)
}



#' @describeIn medoid Calculation of the medoid for each cluster
#' in a partition.
#' @param partition The partition for which the medoids of each cluster
#' should be determined (vector).
#' @param clusters The set of clusters for which the medoid should be
#' calculated (all by default).
#' @returns For `medoids_partition`: a named numeric vector indicating,
#' for each cluster or a subset of them, the identifier of its medoid.
#' @export
medoids_partition <- function(distances, partition, clusters = NULL)
{

  listeClusters <- partition_to_list(partition)

  if (!is.null(clusters))
    listeClusters <- listeClusters[as.character(clusters)]

  nbElements <- lengths(listeClusters)

  medoids <- rep(NA_integer_, length(listeClusters))
  names(medoids) <- names(listeClusters)

  masqueInf2 <- nbElements <= 2L

  if (any(masqueInf2))
    medoids[masqueInf2] <- sapply(listeClusters[masqueInf2], "[[", 1L) # nolint

  if (anyNA(medoids))
  {
    masqueNA <- is.na(medoids)
    medoids[masqueNA] <-
      c3t_sapply(listeClusters[masqueNA],
                 function(indices)
                 {
                   indices[which.min(rowSums(distances[indices, indices]))]
                 }
      )
  }

  medoids
}
