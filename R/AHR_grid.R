#' @include AHR_single.R
#' @include criteria.R

#' @title Agglomerative Hierarchical Regionalisation (AHR)
#' @description This function performs Agglomerative Hierarchical Clustering
#' (AHC)  on a given problem of regionalisation
#' (contiguity/connectivity constraint) with
#' optional size constraints on the clusters. Multiple iterations of AHR can be
#' made if multiple values are given for some
#' parameters.
#' @inheritParams arguments_problem
#' @param contiguity A contiguity matrix or an `igraph` contiguity graph. If not
#' provided, the problem is considered completely contiguous (all elements are
#' neighbors of each other). In that case solving time might be long for each
#' iteration.
#' @param M Maximum size constraint. Must be positive, superior or equal to `m`
#' and large enough for the problem to be feasible. Default is `Inf`
#' (no constraint). Multiple values can be given if different values for the
#' maximum size constraint must be tried.
#' (positive real vector)
#' @param nbTries The number of initialisations to test. The first one will
#' start from a single-element partition, and the remaining ones will
#' be generated randomly respecting connectivity and maximum cluster size
#' constraints. The number of random fusions to perform for each randomly
#' generated partition is determined by `propFusionsInit`.
#' (strictly positive integer)
#' @param propFusionsInit Proportion of fusions (relatively to the number
#' of elements) to be randomly performed for the initialisations. Ignored if
#'   `nbTries = 1`. Default is `1%`. (value in ]0,1[)
#' @param linkages A vector of distance measures to be tested. Can
#' include both function names and actual distance functions. These measures
#' determine the linkage criterion for the hierarchical clustering.
#' If using custom  distance functions, they must take the pairwise distance
#' matrix as an argument and return the linkage distance. Default implemented
#' linkages can be seen with [available_linkages()].
#' @param fusionConstraints Type of constraints to add on the fusions.
#' When a new fusion must be done the algorithm will choose the best fusion
#' in term of the linkage distance respecting the constraint. If no fusion check
#' the constraint a partial relaxation is realized. This kind of constraint can
#' be useful when the problem have a minimum size constraint. Multiple values
#' can be given.
#' The implemented constraints are the following:
#' * `FALSE` (or `NA`): no constraint added
#' * "`singletonMin`": fusion between a cluster of one element and a cluster
#' which doesn't verify the minimum size constraint.
#' * "`singletonAny`" : fusion between a cluster of one element and any other
#' cluster
#' * "`minMin`" : fusion between clusters which do not verify the minimum
#' size constraint
#' * "`minAny`" : fusion between a cluster which do not verify the minimum
#' size constraint and any other one.
#' (vector)
#' @param fusionConstraintModes The way the fusion constraint, if any has been
#' supplied in `fusionConstraints` parameters (ignored otherwise)
#' should be applied. Actually two modes are available:
#' * "`deterministic`": the fusion constraint is applied constantly until it is
#' unfeasible.
#' * "`random`": the fusion constraint is applied randomly with a probability
#' increasing with the number of iterations.
#'
#' Default to "deterministic". (character vector)
#' @param splitConnectedComponents A flag indicating if AHR should be done
#' independently on connected components. This can consequently improve
#' calculation time if the problem is not connected and will not impact the
#' results of the AHR. If `TRUE` the grid will be realized on the different
#' connected components and the return will be a list with one element
#' per component. Ignored if the problem is connected.
#' Default to `FALSE`. (flag)
#' @param criteria A vector of criteria for cluster evaluation to be
#' calculated at the end of the hierarchical clustering. Optional.
#' Use [available_criteria()] to see what criteria are
#' available. (character vector, possibly `NULL`)
#' @param evalLinkages a vector of linkages that will be used by each criterion
#' in `criteria` needing this parameter. By default it is equal to `linkages`.
#' @param minNbClusters Minimum number of clusters allowed. Default is `2`.
#'   Setting a higher value for `minNbClusters` can reduce computation time.
#'   (strictly positive integer)
#' @param maxNbClusters Maximum number of clusters allowed. At the end of
#' each iteration if solutions with a number of cluster inferior to this value
#' exist, those with a superior number of clusters will be removed. Must be
#' superior or equal to `minNbClusters`. Can reduce computation time.
#' Default to `Inf`. (strictly positive integer)
#' @inheritParams parallel_arguments
#' @inheritParams verbose_argument
#' @returns Depending on `splitConnectedComponents` and if the problem is
#' connected or not.
#' If `splitConnectedComponents = FALSE` or the problem is connected, a list
#' with 3 elements:
#' * `results`: A [tibble][tibble::tibble-package] containing partitions given
#' by AHR with some information like  what constraint are respected.
#' If at least one feasible solution have been found only those will be return,
#' with calculated criteria if any given in `criteria`. If at least one
#' criterion is given results will be order by quality regarding the
#' first criterion.
#' Otherwise every solution will be returned, ordered by their size constraint
#' score.
#' * `grid`: the different parameters used for each AHR in a
#' [tibble][tibble::tibble-package].
#' * `initialPartitions`: the different partitions that have been used for
#' the first iteration of the AHRs. Stored in a
#' [tibble][tibble::tibble-package].
#'
#' Otherwise a list with one element per connected component. Each of the
#' elements are a list corresponding to the return of this function for
#' the connected component.
#' @name AHR
#' @details `r badge('experimental')`
#' @example inst/examples/AHR.R
#' @export
#' @seealso [merge_cc_partitions()]
#' @seealso [stats::hclust()] if you don't have any constraint.
#' @importFrom parallel detectCores
#' @importFrom checkmate assertDouble assertFlag
AHR <- function(distances = NULL, contiguity = NULL, sizes = NULL,
                d = NULL, data = NULL,
                m = 0.0, M = Inf,
                standardQuant = FALSE, binarQual = FALSE,
                nbTries = 5L,  propFusionsInit = 0.01,
                linkages = "saut_min",
                fusionConstraints = NA,
                fusionConstraintModes = available_fusion_modes(),
                splitConnectedComponents = FALSE,
                criteria = NULL,
                evalLinkages = linkages,
                minNbClusters = 2L,
                maxNbClusters = Inf,
                parallele = TRUE,
                nbCores = detectCores() - 1L,
                verbose = TRUE)
{

  # Checking arguments
  assertNumericVector(m, lower = 0.0, any.missing = FALSE, all.missing = FALSE)
  assertNumericVector(M, lower = 0.0, any.missing = FALSE, all.missing = FALSE)
  assertFlag(parallele)
  assertFlag(verbose)
  assertFlag(splitConnectedComponents)

  pb <- constructor_pbCon(distances = distances,
                          contiguity = contiguity,
                          sizes = sizes,
                          d = d, data = data,
                          m = max(m), M = min(M),
                          standardQuant = standardQuant,
                          binarQual = binarQual)

  options(c3t_verbose = verbose) # nolint


  if (parallele)
  {
    assertCount(nbCores)
    if (nbCores > 1L)
    {
      c3t_create_clusters(nbCores)
      on.exit(c3t_stop_clusters(), add = TRUE, after = TRUE)
    }
  }

  AHR_pb(pb, nbTries = nbTries, criteria = criteria,
         linkages = linkages,
         evalLinkages = evalLinkages,
         minNbClusters = minNbClusters,
         fusionConstraints = fusionConstraints,
         fusionConstraintModes = fusionConstraintModes,
         splitConnectedComponents = splitConnectedComponents,
         propFusionsInit = propFusionsInit)
}


#' @importFrom tibble as_tibble add_column tibble
#' @importFrom dplyr bind_rows arrange
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom parallel clusterEvalQ clusterExport
#' @importFrom parallel parLapply
#' @importFrom stats hclust
#' @importFrom cli cli_alert
.AHR_grid <- function(pb, grid, # nolint: cyclocomp_linter
                      removeUnfeasible = TRUE,
                      removeTrivialPart = TRUE,
                      maxNbClusters = Inf)
{

  idExp <- NULL # Only used for CHECK warnings
  grid <- arrange(grid, idExp) # nolint
  verbose <- getOption("c3t_verbose", default = TRUE)

  calculRAHgrid <- function(hyperParam)
  {
    set.seed(seed)

    linkage <- hyperParam$linkage
    partitionInit <- hyperParam$partitionInit
    fusionConstraint <- hyperParam$fusionConstraint
    fusionConstraintMode <- hyperParam$fusionConstraintMode
    minNbClusters <- hyperParam$minNbClusters
    idExp <- hyperParam$idExp
    pb$m <- hyperParam$m
    pb$M <- M <- hyperParam$M
    contiguity <- hyperParam$contiguity

    if (M == Inf && !contiguity)
    {
      method <- get_hclust_equivalent(linkage)
      resHclust <- hclust(as.dist(pb), method = method)
      nbClusters <-
        2L:min(maxNbClusters, nrow(pb) - 1L * removeTrivialPart)
      resHclust <- lapply(nbClusters, function(k) cutree(resHclust, k = k))
      resRAH <- tibble(partition  =  resHclust,
                       nbClusters = nbClusters,
                       checkContiguityConst = TRUE,
                       checkMinSizeConst  = TRUE,
                       checkMaxSizeConst  = TRUE,
                       scoreMinSizeConst = 0.0,
                       scoreMaxSizeConst = 0.0,
                       scoreSizeConsts    = 0.0)
    }

    else
    {
      if (!contiguity)
      {
        oldContiguity <- pb$getContiguity()
        pb$setCompleteContiguity()
      }

      resRAH <- AHR_single(pb, linkage = linkage,
                           minNbClusters = minNbClusters,
                           fusionConstraint = fusionConstraint,
                           fusionConstraintMode = fusionConstraintMode,
                           partitionInit = partitionInit[[1L]])

      resRAH <- as_tibble(resRAH)

      if (!contiguity)
      {
        pb$setOldContiguity(oldContiguity)
        rm(oldContiguity)
      }
    }

    if (maxNbClusters < Inf)
    {
      mask <- resRAH$nbClusters <= maxNbClusters
      if (any(mask) && !all(mask))
        resRAH <- resRAH[mask, ]
    }

    if (removeUnfeasible)
    {
      mask <-
        resRAH$checkContiguityConst &
        resRAH$checkMinSizeConst &
        resRAH$checkMaxSizeConst

      if (any(mask))
        resRAH <- resRAH[mask, ]
    }

    if (removeTrivialPart)
    {
      mask <- resRAH$nbClusters > 1L & resRAH$nbClusters < nrow(pb)
      if (any(mask) && !all(mask))
        resRAH <- resRAH[mask, ]
    }

    resRAH <- add_column(resRAH, idExp = idExp, .before = 1L)

    resRAH
  }

  listIdsExp <- grid$idExp
  grid <- split(grid, seq_len(nrow(grid)))
  listErrors <- list()
  seed <- .Random.seed

  results <- c3t_lapply(grid, calculRAHgrid,
                        export = c("seed",
                                   "listErrors",
                                   "removeUnfeasible",
                                   "removeTrivialPart",
                                   "maxNbClusters"),
                        exportIfAbsent = c("pb",
                                           "calculRAHgrid"),
                        elementsToRemove =
                          c("removeUnfeasible",
                            "removeTrivialPart",
                            "maxNbClusters"))

  if (length(listErrors) > 0L)
  {
    warning("Some errors happened during AHR calculations:")
    print(listErrors)
    listeIdsErreurs <- do.call("c", lapply(listErrors, "[[", "idExp"))

    # Remove of the results that sent an error
    if (length(listeIdsErreurs) == length(grid))
    {
      stop("All experiences sent an error")
    }
    else
    {
      mask <- listIdsExp %in% listeIdsErreurs
      results <- results[!mask]
    }
  }

  tailleDFs <- vapply(results, nrow, integer(1L))

  if (verbose)
  {
    alert <-
      gettext("{sum(tailleDFs)} non-trivial regionalisation{?s} obtained")
    cli_alert(alert)
  }


  maskResultatsVides <- tailleDFs == 0L
  if (any(maskResultatsVides))
    results <- results[!maskResultatsVides, ]

  if (removeUnfeasible)
  {
    posFeasible <-
      vapply(results,
             function(resRAH) all(resRAH[1L, c("checkContiguityConst",
                                               "checkMinSizeConst",
                                               "checkMaxSizeConst")]),
             logical(1L))

    if (any(posFeasible))
    {
      results <- results[posFeasible]

      if (verbose)
      {
        success <-
          gettext("{sum(sapply(results, nrow))} feasable partitions obtained")
        cli_alert_success(success)
      }
    }

    else if (verbose)
    {
      danger <- gettext("No feasible solution obtained")
      cli_alert_danger(danger)
    }

  }

  bind_rows(results)
}

#' @importFrom dplyr cross_join distinct
#' @importFrom tibble tibble add_column
#' @importFrom checkmate assertInt assertNumber assertCount assertVector
#' @importFrom tidyr expand_grid
creation_grid_RAH <- function(pb, linkages, nbTries, minNbClusters, # nolint: cyclocomp_linter
                              sizeContGrid,
                              fusionConstraints, fusionConstraintModes,
                              propFusionsInit)

{
  # Checking arguments
  assertInt(minNbClusters, upper = nrow(pb))
  assertCount(nbTries, positive = TRUE)

  if (nbTries > 1L)
  {
    assertProportions(propFusionsInit, len = 1L,
                      zeros.ok = FALSE, ones.ok = FALSE)

  }
  assertVector(linkages, TRUE)
  assertLinkages(linkages)
  linkages <- corresponding_linkage(linkages)

  if (length(linkages) == 0L)
  {
    stop("No linkage has been given")
  }

  n <- pb$n()
  # We generate the initial mergers. The first proposed partition will
  # systematically be the trivial partition (each element is its own region).
  # The other partitions (if any) are generated randomly, based on mergers
  # starting from the trivial partition and respecting the contiguity of the
  # problem. For each randomly generated partition, a certain number of mergers
  # (if possible) are performed, corresponding to `propFusionsInit * n`.
  # This saves computation time and allows testing different initial states.
  nbFusions <- max(ceiling(propFusionsInit * n), 1.0)

  modes <- c("unitary", rep("random", nbTries - 1L))
  initialPartitions <- gen_initial_partitions(pb,
                                              modes = modes,
                                              nbFusions = nbFusions)

  nbFusions <- c(0L, rep(nbFusions, length(initialPartitions) - 1L))

  initialPartitions <- tibble(idpartitionInit =
                                seq_len(length(initialPartitions)),
                              partitionInit = initialPartitions,
                              mode = modes,
                              nbFusions = nbFusions)

  # Creation of the tuning grid
  grid <- expand_grid(linkages, minNbClusters,
                      initialPartitions$idpartitionInit)

  fusionConstraintGrid <- fusion_constraint_mode_grid(fusionConstraints,
                                                      fusionConstraintModes,
                                                      pb)
  grid <- cross_join(grid, fusionConstraintGrid)

  colnames(grid) <- c("linkage", "minNbClusters", "idpartitionInit",
                      "fusionConstraint", "fusionConstraintMode")

  # Looking for redundant size constraints
  sizeContGrid$m <- pb$simplify_min_constraint(sizeContGrid$m)
  sizeContGrid$M <- pb$simplify_max_constraint(sizeContGrid$M)
  if (pb$isConnected())
    sizeContGrid$contiguity <- TRUE

  sizeContGrid   <- distinct(sizeContGrid)

  # Size and contiguity constraints are added to the main grid
  grid <- cross_join(grid, sizeContGrid)

  # Remove of potential redundancies
  maskAucuneContrainte <- grid$M == Inf & !grid$contiguity
  if (any(maskAucuneContrainte))
  {
    if (length(initialPartitions) > 1L)
    {
      maskPartitionNonUnitaire <-
        maskAucuneContrainte & grid$idpartitionInit != 1L
      if (any(maskPartitionNonUnitaire))
      {
        grid <- grid[!maskPartitionNonUnitaire, ]
        maskAucuneContrainte <-  grid$M == Inf & !grid$contiguity
      }
    }
    grid[maskAucuneContrainte, "fusionConstraint"] <- "No"

    linkagehclust <-
      get_hclust_equivalent(grid$linkage[maskAucuneContrainte])
    calculImpossible <- is.na(linkagehclust)
    if (any(calculImpossible))
      grid <- grid[-which(maskAucuneContrainte)[calculImpossible], ]
  }

  maskContrainteFusionInutile <- grid$m == 0.0 &
    grid$fusionConstraint %in% c("singletonMin", "minAny", "minMin")

  if (any(maskContrainteFusionInutile))
    grid[maskContrainteFusionInutile, "fusionConstraint"] <- "No"


  grid <- distinct(grid)

  grid <- add_column(grid,
                     partitionInit =
                       initialPartitions$partitionInit[grid$idpartitionInit],
                     .after = "idpartitionInit")

  grid <- add_column(grid,
                     idExp = seq_len(nrow(grid)),
                     .before = 1L)

  list(grid = grid,
       initialPartitions = initialPartitions)
}

#' @param pb the region partitioning problem with contiguity and / or size
#' constraints to be classified. (`pbCon` object)
#' @param sizeContGrid a data.frame containing 3 columns: `m`, `M`, and
#' `contiguity` to specify the constraint applications. `m` and `M` are positive
#' reals, and `contiguity` is a flag indicating whether to
#' consider the contiguity constraint in the calculation.
#' @importFrom dplyr filter row_number select cross_join
#' @importFrom tibble tibble
#' @importFrom tidyr expand_grid
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom parallel clusterEvalQ clusterExport
#' @importFrom parallel parSapply
#' @importFrom rlang :=
#' @importFrom cli cli_alert
#' @keywords internal
#' @importFrom checkmate assertFlag assertCount assertInt
AHR_pb <- function(pb, nbTries = 5L, criteria = NULL, # nolint: cyclocomp_linter
                   linkages = available_linkages(),
                   evalLinkages = linkages,
                   minNbClusters = 2L,
                   fusionConstraints = available_fusion_constraints(),
                   fusionConstraintModes = available_fusion_modes(),
                   propFusionsInit = 0.01, splitConnectedComponents = TRUE,
                   maxNbClusters = Inf,
                   sizeContGrid = tibble(m = pb$m, M = pb$M, contiguity = TRUE))
{

  assertFlag(splitConnectedComponents)
  verbose <- isTRUE(getOption("c3t_verbose", default = TRUE))

  if (is.infinite(maxNbClusters))
    maxNbClusters <- pb$n()

  else
    assertCount(maxNbClusters, positive = TRUE)

  assertInt(minNbClusters, lower = 1L, upper = maxNbClusters)

  t0 <- Sys.time() # nolint
  if (verbose)
  {
    startingTime <- gettext("Starting time: {t0}")
    cli_alert_info(startingTime)
  }

  n <- pb$n()
  # Checking arguments

  if (length(criteria) > 0L)
    assertCompatibleCriteria(criteria, pb, linkage = evalLinkages)

  if (verbose)
  {
    executionTimeStr <- gettext("Execution time: {Sys.time() - t0}")
    on.exit({cli_alert(executionTimeStr)}, add = TRUE)
  }

  splitConnectedComponents <- splitConnectedComponents &&
    !pb$isConnected()

  if (splitConnectedComponents)
  {
    listeSPb <- pb$split_connected_spb()

    if (verbose)
    {
      info <-
        gettext("the problem has ben divided in {length(listeSPb)} connected components") # nolint: line_length_linter
      cli_alert_info(info)
    }

    resRAHSPb <- list()
    for (spb in listeSPb)
    {
      c3t_clusters_rm("pb")
      res <- AHR_pb(spb, nbTries, criteria,
                    linkages, evalLinkages,
                    minNbClusters,
                    fusionConstraints,
                    fusionConstraintModes,
                    propFusionsInit,
                    splitConnectedComponents = FALSE,
                    maxNbClusters,
                    sizeContGrid)
      resRAHSPb[[length(resRAHSPb) + 1L]] <- res
    }
    c3t_clusters_rm("pb")

    return(resRAHSPb)
  }

  data <- creation_grid_RAH(pb, linkages, nbTries, minNbClusters,
                            sizeContGrid,
                            fusionConstraints, fusionConstraintModes,
                            propFusionsInit)

  rm(linkages, nbTries, minNbClusters, sizeContGrid,
     fusionConstraints, fusionConstraintModes,
     propFusionsInit)

  if (verbose)
  {
    info <- gettext("{nrow(data$grid)} AHC to evaluate")
    cli_alert_info(info)
  }

  # Calculation of the AHC trees. For each tree if some partitions respect
  # the minimum size constraint only those are conserved. Otherwise all are
  # kept.
  results <- .AHR_grid(pb = pb, grid = data$grid,
                       removeUnfeasible = TRUE,
                       removeTrivialPart = TRUE,
                       maxNbClusters = maxNbClusters)

  partitionInit <- NULL # Only used for remove check warning
  data$grid <- select(data$grid, -partitionInit) # nolint

  mask <- results$nbClusters > 1L & results$nbClusters < n
  results <- results[mask, ]
  if (nrow(results) == 0L)
  {
    warning("Imporssible to generate non-trivial regionalisations")
    return(NULL)
  }

  results <- add_column(results, method = "AHR", .before = 1L)

  # If some partitions check all constraints only those are conserved.
  mask <- results$checkContiguityConst &
           results$checkMinSizeConst  &
           results$checkMaxSizeConst

  # If none check all constraints we order results by constraints size score
  # and return it.
  if (!any(mask))
  {
    results <- results[order(results$scoreSizeConsts), ]
    data$results <- results
    return(data)
  }

  results <- results[mask, ]

  # Remove of potential redundant partitions
  tableau <- table(results$nbClusters)
  lignesASupprimmer <- integer(0L)

  for (k in as.numeric(names(tableau)[tableau > 1L]))
  {
    lignesKclusters <- results$nbClusters == k
    partitionsKClusters <- results$partition[lignesKclusters]
    partitionsDupliquees <- duplicated(partitionsKClusters)
    if (any(partitionsDupliquees))
    {
      lignesASupprimmer <- c(lignesASupprimmer,
                             which(lignesKclusters)[partitionsDupliquees])
    }
  }

  if (length(lignesASupprimmer) > 0L)
  {
    results <- filter(results, !row_number() %in% lignesASupprimmer)

    if (verbose)
    {
      alert <-
        gettext("{length(lignesASupprimmer)} redundancie{?s} {?has/have} been removed.") # nolint: line_length_linter
      cli_alert(alert)
    }
  }

  if (length(criteria) > 0L)
  {
    criteria <- criteria_linkages_grid(criteria = criteria,
                                       linkages = evalLinkages)
  }


  for (i in seq_len(nrow(criteria)))
  {
    criterion <- criteria$criterion[[i]]
    linkage <- criteria$linkage[[i]]
    name <- criteria$name[[i]]

    if (verbose)
    {
      alert <- gettext("Calculation of the {criterion} criterion")
      cli_alert(alert)
    }


    critValues <- pb$quality_partition(results$partition,
                                       criterion, linkage,
                                       connected = FALSE,
                                       valueOnly = TRUE)

    results <- add_column(results,
                          !!name := critValues)

    # Results are order by the first criterion,
    # from the best to the worst
    if (i == 1L)
    {
      orderCriterion <- optimal_partitions(critValues,
                                           criterion,
                                           byName = FALSE)
      results <- results[orderCriterion, ]
      rm(orderCriterion)
    }
  }

  data$results <- results
  data

}

#' @param verbose Logical indicating whether to display progress messages.
#' Default is TRUE.
#' @name verbose_argument
NULL
