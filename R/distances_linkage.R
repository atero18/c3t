#' @include distances_elements.R
#' @include pbCon.R
#'
#' Distance inter-clusters de saut
#'
#' Mesure de ressemblance entre deux classes en appliquant une fonction
#  f : R^(n1*n2) -> R (avec n1 et n2 le nombre d'éléments de chaque agrégat)
#  sur l'ensemble des distances entre les différents éléments de la première et
#  les différents éléments de la seconde classe.
#' @param elemsDistances matrice de distance entre les éléments
#' des deux classes.
#' @param comparateur une fonction (R_+)^(n1*n2) -> R_+
#' (avec n1 et n2 le nombre d'éléments de chaque agrégat).
#' défaut : max (fonction)
#'
#' @returns one positive real.
#' @keywords internal
#' @importFrom checkmate assertFunction assertNumber
distance_saut <- function(elemsDistances, comparateur = max)
{
  # Checking arguments
  assertDistanceMatrix(elemsDistances, isComplete = FALSE)
  assertFunction(comparateur)

  res <- .distance_saut(elemsDistances, comparateur)

  assertNumber(res, lower = 0.0, na.ok = FALSE, null.ok = FALSE)

  res
}

.distance_saut <- function(elemsDistances, comparateur = max)
{
  comparateur(elemsDistances)
}

#' Distances between clusters (linkage distances)
#'
#' Some functions for calculating distances between two sets
#' using the distances between elements.
#' @param distances the distance matrix between elements of the first and the
#' second cluster. (matrix)
#' @name linkage_distances
#' @returns The distance between the two clusters. (one positive value)
NULL

#' @describeIn linkage_distances do the single linkage, i.e. take the minimum
#' of the distances between the elements of the first and the elements of the
#' second cluster.
single_linkage <- function(distances) distance_saut(distances, min)
distance_saut_min <- single_linkage

#' @describeIn linkage_distances do the complete linkage, i.e. take the maximum
#' of the distances between the elements of the first and the elements of the
#' second cluster.
complete_linkage <- function(distances) distance_saut(distances, max)
distance_saut_max <- diametre <- complete_linkage


.distance_saut_min <- function(elemsDistances)
{
  .distance_saut(elemsDistances, min)
}

.distance_saut_max <- function(elemsDistances)
{
  .distance_saut(elemsDistances, max)
}

#' @describeIn linkage_distances do the centroid linkage distance. Distance
#' used to obtain a distance on classes from elements composed of
#' quantitative variables. The centroid of each cluster is determined, and the
#' inter-element distance is applied between these two centroids.
#' @param dataCluster1 Data related to elements of the first cluster.
#' Must contain at least one individual and at least one variable. All variables
#' must be quantitative. (data frame)
#' @param dataCluster2 Data related to elements of the second cluster.
#' The restrictions are the same as for `dataCluster1`. Must contain the same
#' variables as the first data frame. (data frame)
#' @param d The distance function between elements. Apart from the first two
#' arguments, all others (if any) must be optional. Must return a positive real
#' number. Default is the euclidean distance. (function)
#' @export
#' @importFrom checkmate assertDataFrame testNumber
centroid_linkage <- function(dataCluster1,
                             dataCluster2,
                             d = euclidean_distance)
{
  # Checking arguments
  assertDataFrame(dataCluster1, types = "numeric",
                  min.rows = 1L, min.cols = 1L)
  assertDataFrame(dataCluster2, types = "numeric",
                  min.rows = 1L, ncols = ncol(dataCluster1))

  if (any(colnames(dataCluster1) != colnames(dataCluster1)))
  {
    stop("both data.frames must have the same variables")
  }

  res <- d(colMeans(dataCluster1),
           colMeans(dataCluster2))

  if (!testNumber(res, lower = 0.0))
  {
    stop("the distance function provided a result that is not a positive number") # nolint: line_length_linter
  }


  return(res)
}

distance_Hausdorff <- function(elemsDistances)
{
  max(max(apply(elemsDistances, 1L, min)),
      max(apply(elemsDistances, 2L, min)))
}

# nocov start
#' What linkages are available?
#'
#' Give the list of implemented linkages.
#' @param description flag indicating if description of implemented
#' linkages should be printed.
#' @returns A character vector containing names of implemented linkages.
#' If `description = TRUE` return will me bade invisibly.
#' @family available parameters
#' @export
#' @keywords internal
#' @importFrom checkmate assertFlag
available_linkages <- function(description = FALSE)
{
  assertFlag(description)

  if (description)
  {
    for (i in seq_len(nrow(LINKAGES)))
    {
      cat(paste0(LINKAGES$linkage[[i]], ": ", LINKAGES$description[[i]]))
      if (i != nrow(LINKAGES))
        cat("\n")
    }

    return(invisible(LINKAGES$linkage))
  }


  LINKAGES$linkage
}
# nocov end

#' @importFrom checkmate assertFlag
corresponding_linkage <- function(linkages, simplify = TRUE)
{

  if (is.null(linkages))
    return(NA_character_)

  assertFlag(simplify)

  if (all(vapply(linkages, is.na, logical(1L))))
  {
    if (simplify)
      return(NA_character_)

    else
      return(linkages)
  }

  maskNA <- vapply(linkages, is.na, logical(1L)) |
    vapply(linkages, is.null, logical(1L))
  if (any(maskNA))
  {
    correspondingNA <- NA_real_
    linkages <- linkages[!maskNA]
  }
  correspondingNA <- NULL

  maskString <- vapply(linkages, is.character, logical(1L))
  maskFunction <- !maskString

  if (any(maskString))
  {
    linkagesStrings <- linkages[maskString]

    if (is.list(linkagesStrings))
      linkagesStrings <- unlist(linkagesStrings)

    correspondingStrings <- ifelse(linkagesStrings %in% LINKAGES$linkage,
                                   linkagesStrings, NaN)


    if (any(is.nan((correspondingStrings))))
    {
      maskNaN <- is.nan(correspondingStrings)
      correspondingStrings[maskNaN] <-
        LINKAGES$linkage[match_str(linkagesStrings[maskNaN], LINKAGES$names)]
    }

    if (simplify)
        correspondingStrings <- unique(correspondingStrings)
  }
  else
    correspondingStrings <- NULL

  if (any(maskFunction))
    correspondingFunctions <- linkages[maskFunction]
  else
    correspondingFunctions <- NULL

  c(correspondingNA, correspondingStrings, correspondingFunctions)
}

test_type_linkage <- function(type, linkage)
{
  !is.null(linkage) &&
    (linkage == type || corresponding_linkage(linkage) == type)
}


function_test_type_linkage <-
  function(type) function(linkage) test_type_linkage(type, linkage)

#' @describeIn linkage_distances Check if the given linkage distance is
#' the single linkage
#' @keywords internal
#' @importFrom purrr partial
is_single_linkage <- partial(test_type_linkage, type = "Single")


#' @describeIn linkage_distances Check if the given linkage distance is
#' the complete linkage
#' @keywords internal
#' @importFrom purrr partial
is_complete_linkage <- partial(test_type_linkage, type = "Complete")

#' @describeIn linkage_distances Check if the given linkage distance is
#' the single or the complete linkage
#' @keywords internal
is_single_or_complete_linkage <-
  function(linkage) is_single_linkage(linkage) || is_complete_linkage(linkage)

#' @describeIn linkage_distances Check if the given linkage distance is
#' the average linkage
#' @keywords internal
#' @importFrom purrr partial
is_average_linkage <- partial(test_type_linkage, type = "Average")

#' @describeIn linkage_distances Check if the given linkage distance is
#' the medoid linkage
#' @keywords internal
#' @importFrom purrr partial
is_medoid_linkage <- partial(test_type_linkage, type = "Medoid")

#' @describeIn linkage_distances Check if the given linkage distance is
#' the centroid linkage
#' @keywords internal
#' @importFrom purrr partial
is_centroid_linkage <- partial(test_type_linkage, type = "Centroid")

#' @describeIn linkage_distances Check if the given linkage distance is
#' the Hausdorff linkage
#' @keywords internal
#' @importFrom purrr partial
is_hausdorff_linkage <- partial(test_type_linkage, type = "Hausdorff")

get_linkage_function <- function(linkage)
{
  if (is.function(linkage))
    return(linkage)

  linkage <- corresponding_linkage(linkage)

  LINKAGES$fun[LINKAGES$linkage == linkage][[1L]]

}

#' @importFrom stats hclust
get_hclust_equivalent <- function(linkages)
{
  assertLinkages(linkages)
  linkages <- corresponding_linkage(linkages)

  LINKAGES$hclust[LINKAGES$hclust == linkages]
}

#' @importFrom checkmate testString
unitary_linkage <- function(linkage)
{
  if (!testString(linkage))
    return(FALSE)

  linkage <- corresponding_linkage(linkage)

  LINKAGES$unitary[LINKAGES$linkage == linkage]
}

#' @importFrom checkmate testString
get_Cppmode_linkage <- function(linkage)
{
  if (!testString(linkage))
    return(FALSE)

  if (linkage %in% LINKAGES$cppMode)
    return(linkage)

  linkage <- corresponding_linkage(linkage)
  cppMode <- LINKAGES$cppMode[LINKAGES$linkage == linkage]

  ifelse(testString(cppMode), cppMode, FALSE)
}


.calcul_distances_inter <- function(pb, partition, distancesACalculer, # nolint: cyclocomp_linter
                                    linkage, Cpp = TRUE, modeEval = 2L, ...)
{

  if (nrow(distancesACalculer) == 0L)
    return(double(0L))

  linkage <- corresponding_linkage(linkage)


  # Particular case when we have the trivial partition et the used distance
  # give for un deux clusters made of une element the distance between those
  # 2 elements.
  n <- length(partition)
  nbClusters <- .nbClusters(partition)
  if (n == nbClusters && unitary_linkage(linkage))
  {
    ordre <- order(partition)
    distancesACalculer[, 1L] <- ordre[distancesACalculer[, 1L]]
    distancesACalculer[, 2L] <- ordre[distancesACalculer[, 2L]]

    return(pb[distancesACalculer[, 1L:2L, drop = FALSE]])
  }


  if (Cpp)
  {
    compCpp <- get_Cppmode_linkage(linkage)
    Cpp <- !isFALSE(compCpp) && !anyNA(pb)
    ##Cpp <- !isFALSE(compCpp) && (!anyNA(pb) || inherits(pb, "DistMat"))
  }

  if (Cpp)
  {
    if (inherits(pb, "DistMat") && anyNA(pb))
    {

      argumentsSupp <- list(...)
      if ("listeClusters" %in% names(argumentsSupp))
        listeClusters <- argumentsSupp[["listeClusters"]]
      else
        listeClusters <- split(seq_len(n), partition)

      col1 <- as.character(distancesACalculer[, 1L])
      col2 <- as.character(distancesACalculer[, 2L])

      indexesMatrices <- lapply(seq_len(nrow(distancesACalculer)),
                                function(k) list(i = listeClusters[[col1[k]]],
                                                 j = listeClusters[[col2[k]]]))

      rm(col1, col2)
      pb$calcul_needed_distances(indexesMatrices)
    }

    return(calculDistancesInter_Cpp(pb, partition,
                                    distancesACalculer,
                                    compCpp))

  }


  else if (is_centroid_linkage(linkage))
  {
    argumentsSupp <- list(...)

    if ("listeClusters" %in% names(argumentsSupp))
      listeClusters <- argumentsSupp[["listeClusters"]]
    else
      listeClusters <- split(seq_len(n), partition)

    distances <- apply(distancesACalculer, 1L,  function(c)
    {
      centroid_linkage(pb$data[listeClusters[[as.numeric(c[1L])]], ],
                       pb$data[listeClusters[[as.numeric(c[2L])]], ],
                       pb$d)
    })

    return(distances)
  }

  else if (is_medoid_linkage(linkage))
  {
    argumentsSupp <- list(...)

    if ("medoides" %in% names(argumentsSupp))
      medoides <- argumentsSupp[["medoides"]]
    else
      medoides <- medoids_partition(pb, partition)

    if (is.null(names(medoides)))
    {
      distancesACalculer[, 1L] <- medoides[distancesACalculer[, 1L]]
      distancesACalculer[, 2L] <- medoides[distancesACalculer[, 2L]]
    }
    else
    {
      carCol1 <- as.character(distancesACalculer[, 1L])
      carCol2 <- as.character(distancesACalculer[, 2L])

      distancesACalculer[, 1L] <- medoides[carCol1]
      distancesACalculer[, 2L] <- medoides[carCol2]
    }

    distances <- pb[distancesACalculer[, 1L:2L, drop = FALSE]]

    return(distances)
  }

  fonctionCalcul <- .calcul_distances_mode_2
  if (modeEval == 1L)
    fonctionCalcul <- .calcul_distances_mode_1
  else if (modeEval == 2L)
    fonctionCalcul <- .calcul_distances_mode_2
  else if (modeEval == 3L)
    fonctionCalcul <- .calcul_distances_mode_3
  else if (modeEval == 4L)
    fonctionCalcul <- .calcul_distances_mode_4

  fonctionCalcul(distances_mat = pb, partition = partition,
                 distancesACalculer = distancesACalculer,
                 linkage = linkage, ...)
}

.calcul_distances_mode_1 <- function(distances_mat, partition,
                                     distancesACalculer, linkage, ...)
{
  n <- length(partition)


  linkage <- get_linkage_function(linkage)

  argumentsSupp <- list(...)

  if ("listeClusters" %in% names(argumentsSupp))
    listeClusters <- argumentsSupp[["listeClusters"]]
  else
    listeClusters <- split(seq_len(n), partition)

  if ("nbElementsClusters" %in% names(argumentsSupp))
    nbElementsClusters <- argumentsSupp[["nbElementsClusters"]]
  else
    nbElementsClusters <- lengths(listeClusters)


  distances <- rep(NA_real_, nrow(distances_mat))

  # Case of linkage with 2 clusters made both of one element.
  masqueTaille <- nbElementsClusters[distancesACalculer[, 1L]] == 1L &
    nbElementsClusters[distancesACalculer[, 2L]] == 1L


  if (any(masqueTaille))
  {

    # nolint start: indentation_linter
    distances[masqueTaille] <-
      distances_mat[
        cbind(do.call("c",
                      listeClusters[distancesACalculer[masqueTaille, 1L]]),
              do.call("c",
                      listeClusters[distancesACalculer[masqueTaille, 2L]])),
      drop = TRUE]
    # nolint end

    masque <- !masqueTaille
    distancesACalculer <- distancesACalculer[masque, , drop = FALSE]
  }
  else
    masque <- rep(TRUE, length(distances))

  if (!any(masque))
    return(distances)


  # Calcul of the remaining distancess
  masqueOrdreInverse <- distancesACalculer[, 1L] > distancesACalculer[, 2L]
  if (any(masqueOrdreInverse))
    distancesACalculer[masqueOrdreInverse, c(2L, 1L)] <-
    distancesACalculer[masqueOrdreInverse, 1L:2L, drop = FALSE]

  ordre <- order(distancesACalculer[, 1L])
  distancesACalculer <- distancesACalculer[ordre, , drop = FALSE]

  clustersPrincipaux <- sort(unique(distancesACalculer[, 1L]))


  res <- lapply(clustersPrincipaux, function(c) {

    clustersCorrespondants <-
      distancesACalculer[distancesACalculer[, 1L] == c, 2L]
    listeClustersCorrespondants <- listeClusters[clustersCorrespondants]

    data <- distances_mat[listeClusters[[c]],
                          do.call("c", listeClustersCorrespondants)]

    res <- split(data,
                 matrix(rep(seq_along(clustersCorrespondants),
                            times = lengths(listeClustersCorrespondants),
                            nrow = 1L)))

    vapply(res, linkage, numeric(1L))
  })

  res <- do.call("c", res)
  distances[masque][ordre] <- res

  distances
}

.calcul_distances_mode_2 <- function(distances_mat, partition,
                                     distancesACalculer, linkage, ...)
{
  # Evaluation by taking directly all distances under the form of a
  # then by splitting the list in submatrices
  if (!is_pbCon(distances_mat))
    return(.calcul_distances_mode_4(distances_mat, partition,
                                    distancesACalculer, linkage, ...))

  linkage <- get_linkage_function(linkage)
  n <- length(partition)

  argumentsSupp <- list(...)

  if ("listeClusters" %in% names(argumentsSupp))
    listeClusters <- argumentsSupp[["listeClusters"]]
  else
    listeClusters <- split(seq_len(n), partition)

  distancesACalculer[1L, 1L] <- as.character(distancesACalculer[1L, 1L])


  # For each cluster pair we get all the following lines
  listeIndices <- apply(distancesACalculer, 1L, function(c)
    list(i = listeClusters[[c[1L]]], j = listeClusters[[c[2L]]]))

  # Final conversion to submatrices
  sousMatrices <- distances_mat[listeIndices, drop = TRUE]

  vapply(sousMatrices, linkage, numeric(1L))
}

.calcul_distances_mode_3 <- function(distances_mat, partition,
                                     distancesACalculer, linkage, ...)
{
  # Calcul by getting all needed distances in one shot, stocking them in
  # a vector. This vector is them splitted before being send to the linkage
  # distance function (this impose that the linkage distance does not need
  # a matrix format).

  linkage <- get_linkage_function(linkage)

  argumentsSupp <- list(...)

  if ("listeClusters" %in% names(argumentsSupp))
    listeClusters <- argumentsSupp[["listeClusters"]]
  else
    listeClusters <- split(seq_len(n), partition)

  # Creation of a 2 column matrix indicating position of different elements
  # of distances to get.
  indexs <- apply(distancesACalculer, 1L, function(c)
  {
    c <- as.numeric(c)
    cluster1Points <- listeClusters[[c[1L]]]
    cluster2Points <- listeClusters[[c[2L]]]

    cbind(rep.int(cluster1Points, times = length(cluster2Points)),
          rep(cluster2Points, each = length(cluster1Points)))

  }, simplify = FALSE)

  tailleVecs_vec <- vapply(indexs, nrow, integer(1L))
  tailleVecsCum <- c(0L, cumsum(tailleVecs_vec)) + 1L

  indexs <- do.call("rbind", indexs)
  data_vec <- distances_mat[indexs, ]

  vapply(seq_len(nrow(distancesACalculer)),
         function(k) linkage(data_vec[tailleVecsCum[k] +
                                        0L:(tailleVecs_vec[k] - 1L)]),
         numeric(1L))
}

.calcul_distances_mode_4 <- function(distances_mat, partition,
                                     distancesACalculer, linkage, ...)
{
  # Classic calcul: recuperation of a submatrix per linkage


  n <- length(partition)

  linkage <- get_linkage_function(linkage)

  argumentsSupp <- list(...)

  if ("listeClusters" %in% names(argumentsSupp))
    listeClusters <- argumentsSupp[["listeClusters"]]

  else
    listeClusters <- split(seq_len(n), partition)

  distancesACalculer[, 1L] <- as.character(distancesACalculer[, 1L])
  distancesACalculer[, 2L] <- as.character(distancesACalculer[, 2L])
  sousMatricesDist <- apply(distancesACalculer, 1L, function(c)
    distances_mat[listeClusters[[c[1L]]],
                  listeClusters[[c[2L]]]], simplify = FALSE)

  vapply(sousMatricesDist, linkage, numeric(1L))
}


calculDistancesInter_Cpp <- function(x, partition,
                                     indexs_mat, comp_str = "min")
{
  comp_str <- get_Cppmode_linkage(comp_str)

  if (inherits(x, "DistMat"))
    x <- x$distances

  if (inherits(x, "SymMMat"))
  {
    return(distanceInterSymMMat(x$values, partition,
                                indexs_mat, comp_str))
  }

  else if (inherits(x, "SymVMat"))
  {
    aDefautDiag <-  !is.null(x$defaultDiag)
    defaultDiag <- ifelse(aDefautDiag, x$defaultDiag, 0.0)

    return(distanceInterSymVMat(x$values, partition,
                                indexs_mat, nrow(x),
                                aDefautDiag, defaultDiag, comp_str))
  }
  else if (is.matrix(x))
  {
    return(distanceInterSymMMat(x, partition,
                                indexs_mat, comp_str))
  }
  else
  {
    stop("Actually this function doesn't work for other type of data")
  }
}


# Linkage distance update -------------------------------------------------

#' @importFrom checkmate assertInt assertNumber assertIntegerish
lance_williams_update <- function(linkDistances,
                                  cluster1, cluster2,
                                  alpha1, alpha2,
                                  beta,
                                  gamma,
                                  coeffs = NULL,
                                  usefulClusters = seq_len(nrow(linkDistances)),
                                  ...)
{
  assertInt(cluster1, lower = 1L, upper = nrow(linkDistances))
  assertInt(cluster2, lower = 1L, upper = nrow(linkDistances))
  if (cluster1 == cluster2)
  {
    stop("`cluster1` and `cluster2` must be different")
  }

  if (!is.null(coeffs))
  {
    alpha1 <- coeffs[["alpha1"]]
    alpha2 <- coeffs[["alpha2"]]
    beta   <- coeffs[["beta"]]
    gamma <- coeffs[["gamma"]]
  }

  assertNumber(alpha1, finite = TRUE)
  assertNumber(alpha2, finite = TRUE)
  assertNumber(beta, finite = TRUE)
  assertNumber(gamma, finite = TRUE)

  if (!is.null(match.call()$usefulClusters))
  {
    assertIntegerish(usefulClusters, lower = 1L, upper = nrow(distances),
                     unique = TRUE, any.missing = FALSE, all.missing = FALSE)
  }

  linkDistances <- .lance_williams_update(linkDistances,
                                          cluster1, cluster2,
                                          alpha1, alpha2,
                                          beta,
                                          gamma,
                                          usefulClusters =
                                            seq_len(nrow(linkDistances)),
                                          ...)

  if (any(linkDistances[cluster1, usefulClusters] < 0.0, na.rm = TRUE))
  {
    stop("The values used for Lance-Williams update gave at least one negative linkage distance") # nolint: line_length_linter
  }

  linkDistances
}

.lance_williams_update <- function(linkDistances,
                                   cluster1, cluster2,
                                   alpha1, alpha2,
                                   beta,
                                   gamma,
                                   usefulClusters = seq_len(nrow(distances)),
                                   ...)

{
  # Particular cases
  # -- single and complete linkages
  if (alpha1 == 0.5 && alpha2 == 0.5 && beta == 0.0 && gamma %in% c(-0.5, 0.5))
  {
    linkage <- ifelse(gamma == -0.5, "single", "complete")
    return(.update_linkage(linkDistances,
                           cluster1, cluster2,
                           linkage,
                           usefulClusters, ...))
  }

  k <- nrow(linkDistances)

  # Check for clusters that don't exist anymore (with NaN values)
  usefulClusters <-
    usefulClusters[!(usefulClusters %in% c(cluster1, cluster2))]
  usefulClusters <-
    usefulClusters[!is.nan(linkDistances[cluster1, usefulClusters])]

  # Apply the Lance-Williams update formula
  update <-
    alpha1 * linkDistances[cluster1, usefulClusters] +
    alpha2 * linkDistances[cluster2, usefulClusters] +
    beta   * linkDistances[cluster1, cluster2] +
    gamma  * abs(linkDistances[cluster1, usefulClusters] -
                   linkDistances[cluster2, usefulClusters])


  linkDistances[cluster1, usefulClusters] <- update

  if (!inherits(linkDistances, c("DistMat", "AbstractSymMat")))
    linkDistances[usefulClusters, cluster1] <- update

  # Fix all distances with cluster 2 to NaN
  idsNoCluster2 <- seq_len(k)
  idsNoCluster2 <- idsNoCluster2[idsNoCluster2 != cluster2]

  linkDistances[cluster2, idsNoCluster2] <- NaN

  if (!inherits(linkDistances, c("DistMat", "AbstractSymMat")))
  {
    linkDistances[idsNoCluster2, cluster2] <-
      linkDistances[cluster2, idsNoCluster2]
  }

  linkDistances

}
update_single_linkage <- function(linkDistances,
                                  cluster1, cluster2,
                                  usefulClusters =
                                    seq_len(nrow(linkDistances)),
                                  elemDistances = NULL,
                                  partitionBefore = NULL,
                                  ...)
{
  update_linkage(linkDistances, cluster1, cluster2,
                 "single", usefulClusters,
                 elemDistances,
                 partitionBefore, ...)
}


.update_single_linkage <- function(linkDistances,
                                   cluster1, cluster2,
                                   usefulClusters =
                                     seq_len(nrow(linkDistances)),
                                   elemDistances = NULL,
                                   partitionBefore = NULL,
                                   ...)
{
  .update_linkage(linkDistances, cluster1, cluster2,
                  "single", usefulClusters, elemDistances,
                  partitionBefore, ...)
}

update_complete_linkage <- function(linkDistances,
                                    cluster1, cluster2,
                                    usefulClusters =
                                      seq_len(nrow(linkDistances)),
                                    elemDistances = NULL,
                                    partitionBefore = NULL,
                                    ...)
{
  update_linkage(linkDistances, cluster1, cluster2,
                 "complete", usefulClusters,
                 elemDistances, partitionBefore, ...)
}


.update_complete_linkage <- function(linkDistances,
                                     cluster1, cluster2,
                                     usefulClusters =
                                       seq_len(nrow(linkDistances)),
                                     elemDistances = NULL,
                                     partitionBefore = NULL,
                                     ...)
{
  .update_linkage(linkDistances, cluster1, cluster2,
                  "complete", usefulClusters,
                  elemDistances, partitionBefore, ...)
}

update_average_linkage <- function(linkDistances,
                                   cluster1, cluster2,
                                   usefulClusters =
                                     seq_len(nrow(linkDistances)),
                                   elemDistances = NULL,
                                   partitionBefore = NULL,
                                   nbElementsClustersBefore = NULL,
                                   ...)
{
  update_linkage(linkDistances, cluster1, cluster2,
                 "average", usefulClusters,
                 elemDistances, partitionBefore,
                 nbElementsClustersBefore, ...)
}
.update_average_linkage <- function(linkDistances,
                                    cluster1, cluster2,
                                    usefulClusters =
                                      seq_len(nrow(linkDistances)),
                                    elemDistances = NULL,
                                    partitionBefore = NULL,
                                    nbElementsClustersBefore = NULL,
                                    ...)
{

  .update_linkage(linkDistances, cluster1, cluster2,
                  "average", usefulClusters,
                  elemDistances, partitionBefore,
                  nbElementsClustersBefore, ...)

}

.update_linkage <- function(linkDistances, # nolint: cyclocomp_linter
                            cluster1, cluster2,
                            linkage = "single",
                            usefulClusters =
                              seq_len(nrow(linkDistances)),
                            elemDistances = NULL,
                            partitionBefore = NULL,
                            nbElementsClustersBefore = NULL,
                            ...)
{
  k <- nrow(linkDistances)

  # Check for clusters that don't exist anymore (with NaN values)
  usefulClusters <-
    usefulClusters[!(usefulClusters %in% c(cluster1, cluster2))]
  usefulClusters <-
    usefulClusters[!is.nan(linkDistances[cluster1, usefulClusters])]

  noSymmetricClass <- !inherits(linkDistances, c("DistMat", "AbstractSymMat"))

  # Check if any missing values in the distances including cluster1
  # or cluster2
  if (length(usefulClusters) > 0L)
  {
    # `NA` values means that new contiguities are possibles with # nolint: commented_code_linter
    # cluster1+cluster2 # nolint: commented_code_linter
    naDistcluster1 <-
      usefulClusters[is.na(linkDistances[cluster1, usefulClusters]) &
                   !is.nan(linkDistances[cluster1, usefulClusters])] # nolint: indentation_linter

    # Calculate NA values for cluster2 is only useful if there is an
    # update function for distances
    if (is_single_or_complete_linkage(linkage) ||
        (is_average_linkage(linkage) &&
         (!is.null(nbElementsClustersBefore) || !is.null(partitionBefore))))
    {
      naDistcluster2 <-
        usefulClusters[is.na(linkDistances[cluster2, usefulClusters]) &
                     !is.nan(linkDistances[cluster2, usefulClusters])] # nolint: indentation_linter
    }
    else
      naDistcluster2 <- NULL


    # if some NA have been found
    if (length(naDistcluster1) > 0L || length(naDistcluster2) > 0L)
    {
      # Case when the distances between elements have not
      # been given
      if (is.null(elemDistances) || is.null(partitionBefore))
      {
        usefulClusters <- setdiff(usefulClusters,
                                  c(naDistcluster1, naDistcluster2))
      }

      # When distances between elements are given missing values
      # can be calculated
      else
      {
        distancesToCalculate <- matrix(nrow = 0L, ncol = 2L)
        if (length(naDistcluster1) > 0L)
          distancesToCalculate <- rbind(distancesToCalculate,
                                        cbind(cluster1, naDistcluster1))

        if (length(naDistcluster2) > 0L)
          distancesToCalculate <- rbind(distancesToCalculate,
                                        cbind(cluster2, naDistcluster2))


        linkDistances[distancesToCalculate] <-
          .calcul_distances_inter(elemDistances, partitionBefore,
                                  distancesToCalculate, linkage)

        if (noSymmetricClass)
        {
          linkDistances[distancesToCalculate[, c(2L, 1L), drop = FALSE]] <-
            linkDistances[distancesToCalculate]
        }
      }
    }
  }


  if (length(usefulClusters) > 0L)
  {
    if (is_single_or_complete_linkage(linkage))
    {
      if (is_single_linkage(linkage))
        comp <- get(">")

      else
        comp  <- get("<")

      # Find which linkage distances must be updated
      masque <- comp(linkDistances[cluster1, usefulClusters],
                     linkDistances[cluster2, usefulClusters])

      # Update of the linkage distances
      if (any(masque))
      {

        update <- linkDistances[cluster2, usefulClusters[masque]]
        linkDistances[cluster1, usefulClusters[masque]] <- update


        if (noSymmetricClass)
        {
          linkDistances[usefulClusters[masque], cluster1] <- update
        }
      }
    }
    else if (is_average_linkage(linkage) &&
             (!is.null(nbElementsClustersBefore) || !is.null(partitionBefore)))
    {

      if (is.null(nbElementsClustersBefore))
      {
        cardCluster1 <- sum(partitionBefore == cluster1)
        cardCluster2 <- sum(partitionBefore == cluster2)
      }
      else
      {
        cardCluster1 <- nbElementsClustersBefore[cluster1]
        cardCluster2 <- nbElementsClustersBefore[cluster2]
      }

      if (cardCluster1 == 0L || cardCluster2 == 0L)
      {
        stop("`cluster1` and `cluster2` must have elements")
      }

      update <- cardCluster1 * linkDistances[cluster1, usefulClusters] +
        cardCluster2 * linkDistances[cluster2, usefulClusters]
      update <- update / (cardCluster1 + cardCluster2)

      linkDistances[cluster1, usefulClusters] <- update

      if (noSymmetricClass)
      {
        linkDistances[usefulClusters, cluster1] <-
          linkDistances[cluster1, usefulClusters]
      }

    }
    else
    {
      indexes <- cbind(cluster1, usefulClusters)
      newPartition <- partitionBefore
      newPartition[newPartition == cluster2] <- cluster1

      newDistances <- .calcul_distances_inter(linkDistances,
                                              newPartition,
                                              indexes, ...)

      linkDistances[cluster1, usefulClusters] <- newDistances

      if (noSymmetricClass)
        linkDistances[usefulClusters, cluster1] <- newDistances
    }
  }

  # Fix all distances with cluster 2 to NaN
  idsNoCluster2 <- seq_len(k)[-cluster2]

  linkDistances[cluster2, idsNoCluster2] <- NaN

  if (noSymmetricClass)
  {
    linkDistances[idsNoCluster2, cluster2] <-
      linkDistances[cluster2, idsNoCluster2]
  }

  linkDistances
}

#' @importFrom checkmate assertInt assertIntegerish
update_linkage <- function(linkDistances, # nolint: cyclocomp_linter
                           cluster1, cluster2,
                           linkage = "single",
                           usefulClusters =
                             seq_len(nrow(linkDistances)),
                           elemDistances = NULL,
                           partitionBefore = NULL,
                           nbElementsClustersBefore = NULL,
                           ...)
{
  assertInt(cluster1, lower = 1L, upper = nrow(linkDistances))
  assertInt(cluster2, lower = 1L, upper = nrow(linkDistances))
  if (cluster1 == cluster2)
  {
    stop("`cluster1` and `cluster2` must be different")
  }

  if (!is.null(elemDistances) && !is.null(partitionBefore))
    assertIntegerish(partitionBefore, len = nrow(elemDistances))

  if (!is.null(match.call()$usefulClusters))
  {
    assertIntegerish(usefulClusters, lower = 1L, upper = nrow(linkDistances),
                     unique = TRUE, any.missing = FALSE, all.missing = FALSE)
  }

  if (!is_single_or_complete_linkage(linkage) && !is_average_linkage(linkage))
  {
    stop("`linkage` must be single, complete or average")
  }

  if (is_average_linkage(linkage) &&
      is.null(partitionBefore) &&
      is.null(nbElementsClustersBefore))
  {
    stop("at least one of `partitionBefore` and `nbElementsClustersBefore` must be specified") # nolint: line_length_linter
  }

  if (!is.null(nbElementsClustersBefore))
    assertIntegerish(nbElementsClustersBefore, len = nrow(linkDistances))

  eval(body(.update_linkage), envir = environment())
}
