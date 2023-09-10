#' @param pb problem with connectivity and/or size constraint
#' @param regionalisation feasible regionalisation to optimize.
#' @param maxIt maximum number of allowed iterations. Default is `Inf`.
#' (strictly positive integer)
#' @param linkage linkage distance to use if `criterion` in {"`Dunn`", "`AHC`"}
#' @returns A list containing the following elements:
#' * `statut`: state of improvement. Indicates whether an improvement
#' could be made or not.
#' * `regionalisation`: the new regionalisation. Identical to the input
#' argument if no improvement could be made.
#' * `it`: number of improving iterations performed
#' * `criterion`: name of the criterion used
#' * `initialCriterion`: initial value of the criterion
#' * `finalCriterion`: final value of the criterion, at least as good as at the
#' beginning of the algorithm.
#' @name enhancer
NULL

#' Improvement of an existing feasible solution
#'
#'
#' @description Algorithm to improve (according to a certain criterion) a
#' solution that is feasible for a certain classification problem with
#' connectivity and size constraints. `r badge('experimental')`
#' @inheritParams arguments_problem
#' @inheritParams enhancer
#' @inheritParams enhancer_grid_creation
#' @inheritParams grid_enhancer_pb
#' @inheritParams verbose_argument
#' @inheritParams parallel_arguments
#' @importFrom parallel detectCores
#' @returns a [tibble][tibble::tibble-package] with one row per try.
#' For each row the following variables:
#' * `criterion`: name of the criterion used for improvement.
#' * `linkage`: type of linkage distance used (`NA` if this argument is
#' irrelevant for the actual criterion).
#' * `sampleSize`: size of the sample for the calculation of the criterion
#' (`NA` if irrelevant).
#' * `statut`: state of improvement. Indicates whether an improvement
#' could be made or not.
#' * `iterations`: number of improving iterations performed.
#' * `regionalisationOpti`: the new regionalisation. Identical to the input
#' argument if no improvement could be made.
#' * one column per criterion indicated in  `critereEvaluation`. If some of
#' those criteria use a linkage distance, there will be one column per linkage
#' distance given in `linkage` and per criterion.
#  @aliases enhance enhancer amelioration ameliorant
#' @export
#' @references Marc Christine and Michel Isnard. "Un algorithme de regroupement
#' d'unités statistiques selon certains critères de similitudes"
#' Insee Méthodes, 2000, p. 50`
#' @importFrom checkmate assertFlag
enhance_feasible <-
  function(regionalisation,
           distances = NULL,
           contiguity = NULL,
           sizes = NULL,
           d = NULL, data = NULL,
           m = 0.0, M = Inf,
           standardQuant = FALSE, binarQual = FALSE,
           enhanceCriteria = c("AHC", "Silhouette", "Dunn"),
           linkages = "saut max",
           evaluationCriteria = enhanceCriteria,
           maxIt = Inf,
           parallele = TRUE, nbCores = detectCores() - 1L,
           verbose = TRUE)
  {

    assertFlag(verbose)
    options(c3t_verbose = verbose) # nolint: undesirable_function_linter

    pb <- constructor_pbCon(distances = distances,
                            contiguity = contiguity,
                            sizes = sizes,
                            d = d, data = data,
                            m = m, M = M,
                            standardQuant = standardQuant,
                            binarQual = binarQual)

    grille <- enhancer_grid_creation(pb,
                                     enhanceCriteria,
                                     linkages)

    parallele <- parallele && nrow(grille) > 1L

    if (parallele)
    {
      nbCores <- min(nbCores, nrow(grille))
      c3t_create_clusters(nbCores)
      on.exit(c3t_stop_clusters(), add = TRUE, after = TRUE)
    }

    resAmelioration <- grid_enhancer_pb(pb, regionalisation,
                                        grille, evaluationCriteria,
                                        maxIt)

    options(c3t_verbose = FALSE) # nolint: undesirable_function_linter

    resAmelioration
  }


#' @importFrom checkmate assertString assertCount
enhancer_pb <- function(pb, regionalisation, criterion = "Dunn", # nolint: cyclocomp_linter
                        maxIt = Inf,
                        sampleSize = ceiling(0.8 * pb$n()),
                        linkage = "single")
{
  # Verification of arguments
  assertString(criterion)

  if (toupper(criterion) == "AHC")
    globalCriterion <- FALSE

  else
  {
    assertCompatibleCriterion(criterion, pb, linkage = linkage)
    globalCriterion <- TRUE
    criterion <- corresponding_criterion(criterion)
    updatableCriterion <- updatable_criterion(criterion)
    if (updatableCriterion)
      nameDataCrit <- get_list_name_criterion(criterion)


    ## ASSERT LINKAGE
  }

  if (!is.infinite(maxIt))
    assertCount(maxIt, positive = TRUE)

  # Verifying that the given partition is a feasible solution
  if (!pb$isFeasibleSolution(regionalisation))
  {
    stop("Given solution is not a feasible solution")
  }

  regionalisation <- standardize_partition(regionalisation)

  # Retrieving the sizes of each region
  regionsSizes <- .clusters_sizes(regionalisation, pb$sizes)
  regionsCardinals <- table(regionalisation)

  t0 <- Sys.time()
  statut <- NULL

  nbRegions <- .nbRegions(regionalisation)
  n <- pb$n()

  # Drawing an arbitrary permutation of all elements
  listeElements_vec <- sample(seq_len(n))

  dataElementsTransferables <- NULL
  oldRegionalisation <- regionalisation

  calcul_elementsTransferables <- function()
  {

    # Calculating or updating the list of transferable elements in
    # the current state
    if (it == 0L)
    {
      assign("dataElementsTransferables",
             transferable_elements_pb(regionalisation, pb = pb),
             envir = parent.frame(1L))
    }
    else
    {
      assign("dataElementsTransferables",
             maj_elements_transferables(oldRegionalisation,
                                        regionalisation,
                                        dataElementsTransferables,
                                        donor,
                                        receiver,
                                        pb = pb),
             envir = parent.frame(1L))
    }

    assign("dataElementsTransferables",
           remove_transferables(pb,
                                regionalisation,
                                dataElementsTransferables,
                                regionsSizes,
                                rmFromSingleton = TRUE,
                                rmCreationMin = TRUE,
                                rmCreationMax = TRUE,
                                rmNormToNorm = FALSE,
                                smallClusters = integer(0L),
                                bigClusters = integer(0L),
                                nbElemsClusters = regionsCardinals),
           envir = parent.frame(1L))


    dataElementsTransferables$transferable
  }

  it <- 0L
  elementsTransferables <- calcul_elementsTransferables()

  if (globalCriterion)
  {
    if (updatableCriterion)
    {
      dataCrit <- pb$quality_partition(regionalisation,
                                       criterion,
                                       linkage,
                                       connected = FALSE,
                                       valueOnly = FALSE)

       initialCriterion <- actualCrit <-
        finalCriterion <- dataCrit[[nameDataCrit]]

    }
    else
    {
      initialCriterion <-
        actualCrit <-
        finalCriterion <-
        pb$quality_partition(regionalisation,
                             criterion,
                             linkage,
                             connected = FALSE,
                             valueOnly = TRUE)
    }
  }
  else
    initialCriterion <- actualCrit <- finalCriterion <- NA_real_


  # We continue transfers as long as it's possible and beneficial
  k <- 1L
  subList <- listeElements_vec[listeElements_vec %in% elementsTransferables]

  tryCatch(
    interrupt = # nolint: indentation_linter
      function(cond)
      {
        message("interruption in the improvement")
        assign("statut", "interruption", envir = parent.frame(1L))
      },
    expr = # nolint: indentation_linter
      {
        while (k <= length(subList) && it <= maxIt)
        {
          # We take the k-th element from the list of randomly generated
          # elements, which is present in the list of transferable elements.
          x1 <- subList[k]

          # If it appears multiple times, we randomly shuffle the transfer
          # possibilities.
          lignes <- which(dataElementsTransferables$frontier[, "x1"] == x1)
          if (length(lignes) > 1L)
            lignes <- sample(lignes)

          transfert <- FALSE
          donor <- dataElementsTransferables$frontier[lignes[1L], "cluster1"]

          if (!globalCriterion)
          {
            partitionSeparation <- regionalisation
            partitionSeparation[x1] <- nbRegions + 1L

            distanceDonor <-
              .calcul_distances_inter(pb,
                                      partitionSeparation,
                                      matrix(c(donor,
                                               nbRegions + 1L),
                                             nrow = 1L, ncol = 2L),
                                      linkage)

          }

          numLigne <- 1L
          while (!transfert && numLigne <= length(lignes))
          {
            l <- lignes[numLigne]
            receiver <- dataElementsTransferables$frontier[l, "cluster2"]


            nouvelleRegionalisation <- regionalisation
            nouvelleRegionalisation[x1] <- receiver

            if (globalCriterion)
            {
              if (updatableCriterion)
              {
                newDataCrit <-
                  pb$update_quality_partition(regionalisation,
                                              criterion,
                                              dataCrit,
                                              donor,
                                              receiver,
                                              x1,
                                              linkage)

                newCrit <- newDataCrit[[nameDataCrit]]

              }
              else
              {
                newCrit <- pb$quality_partition(nouvelleRegionalisation,
                                                criterion,
                                                linkage,
                                                connected = FALSE,
                                                valueOnly = TRUE)
              }

              if (optimal_partitions(c(actualCrit, newCrit),
                                     criterion)[1L] == 2L)
              {
                transfert <- TRUE
                actualCrit <- newCrit
                if (updatableCriterion)
                  dataCrit <- newDataCrit
              }
            }


            # Agglomerative Hierarchical Clustering (AHC) Criterion:
            # we detach the point from its cluster and check if its linkage
            # distance is smaller with the potential receiving cluster
            # than with its current cluster.
            else
            {
              distanceReceiver <-
                .calcul_distances_inter(pb,
                                        partitionSeparation,
                                        matrix(c(receiver,
                                                 nbRegions + 1L),
                                               nrow = 1L, ncol = 2L),
                                        linkage)

              if (distanceReceiver < distanceDonor)
                transfert <- TRUE
            }

            numLigne <- numLigne + 1L
          }

          # If transfert happened we update regionalisation data
          if (transfert)
          {
            oldRegionalisation <- regionalisation

            regionalisation <- nouvelleRegionalisation

            regionsSizes[donor] <-
              regionsSizes[donor] - pb$sizes[x1]

            regionsCardinals[donor] <-
              regionsCardinals[donor] - 1L

            regionsSizes[receiver] <-
              regionsSizes[receiver] + pb$sizes[x1]

            regionsCardinals[receiver] <-
              regionsCardinals[receiver] + 1L

            # Update of the transferable elements
            elementsTransferables <- calcul_elementsTransferables()

            if (length(elementsTransferables) > 0L)
            {
              subList <- listeElements_vec[listeElements_vec %in%
                                               elementsTransferables]
            }
            else
              subList <- NULL

            k <- 1L
            it <- it + 1L
          }

          else
            k <- k + 1L
        }
      })


  if (is.null(statut) && it == 0L)
    statut <- "non_ameliorable"

  else if (is.null(statut))
    statut <- "amelioration"


  if (it > 0L)
  {
    regionalisation <- standardize_partition(regionalisation)

    finalCriterion <- actualCrit
  }

  list(statut = statut, regionalisation = regionalisation, it = it,
       criterion = criterion, initialCriterion = initialCriterion,
       finalCriterion = finalCriterion,
       tempsCalcul = Sys.time() - t0)
}

#' @description Create a grid composed of the different parameters used by an
#' iteration of the enhancer. Depending of the `pb` object the grid will be
#' simplified. Also some criteria might not be available.
#' @param enhanceCriteria A vector of criteria used for the enhancement of the
#' actual feasible solution.
#' Currently available choices are those in [available_criteria()],
#' plus "AHC" (depends of the `linkages` parameter). Compared to others AHC
#' doesn't improve a global criterion but do this
#' locally, hoping to reduce computing time. Regarding to this criterion
#' a feasible solution, built by move a unique element from a cluster to another
#' is better if the element is closer to the other cluster than it's actual
#' (depending of some linkage).
#' @param linkages Vector of linkage distances used when a criterion ("`Dunn`",
#' "`AHC`") needs it.
#' @keywords internal
#' @importFrom tidyr expand_grid
#' @importFrom dplyr distinct
enhancer_grid_creation <- function(
    pb,
    enhanceCriteria = c(available_criteria(), "AHC"),
    linkages = "saut_min")
{
  linkages <- corresponding_linkage(linkages, simplify = TRUE)

  enhanceCriteria <- toupper(unique(enhanceCriteria))

  if ("AHC" %in% enhanceCriteria && length(enhanceCriteria) > 1L)
  {
    assertCompatibleCriteria(enhanceCriteria[enhanceCriteria != "AHC"], pb,
                             linkage = linkages)
    enhanceCriteria[enhanceCriteria != "AHC"] <-
      corresponding_criterion(enhanceCriteria[enhanceCriteria != "AHC"])
  }
  else if (!("AHC" %in% enhanceCriteria))
  {
    assertCompatibleCriteria(enhanceCriteria, pb, linkage = linkages)
    enhanceCriteria <- corresponding_criterion(enhanceCriteria)
  }

  grille <- expand_grid(criterion = enhanceCriteria,
                        linkage = linkages)


  noLinkageCriteria <- CRITERIA$criterion[!CRITERIA$needsLinkage]
  if (length(noLinkageCriteria) > 0L &&
      any(enhanceCriteria %in% noLinkageCriteria))
  {
    grille[grille$criterion %in% noLinkageCriteria, "linkage"] <- NA_character_
  }

  distinct(grille)

}

#' @description Function applying on a partition the enhancer algorithm
#' for each set of parameters given in a grid.
#' @param evaluationCriteria criteria used for comparison after enhancement.
#' They are evaluated on each feasible solution given by each criterion used
#' for enhancement. Must be a vector composed of the available criteria in c3t.
#' For the Dunn index there will be one criterion per linkage given.
#' See [available_criteria()].
#' @keywords internal
#' @importFrom tibble as_tibble add_column tibble
#' @importFrom dplyr  bind_cols select pull
#' @importFrom purrr partial
#' @importFrom cli cli_alert
grid_enhancer_pb <- function(pb, # nolint: cyclocomp_linter
                             regionalisation,
                             grille,
                             evaluationCriteria,
                             maxIt = 1500L,
                             ...)
{
  # Checking arguments
  if (!pb$isFeasibleSolution(regionalisation))
  {
    stop("`regionalisation` is not a feasible solution")
  }

  verbose <- isTRUE(getOption("c3t_verbose", default = TRUE))

  if (length(evaluationCriteria) > 0L)
  {
    evaluationCriteria <- unique(toupper(evaluationCriteria))
    evaluationCriteria <- setdiff(evaluationCriteria, "AHC")
  }


  if (length(evaluationCriteria) > 0L)
  {
    linkages <- unique(grille$linkage[!is.na(grille$linkage)])
    assertCompatibleCriteria(evaluationCriteria, pb, linkage = linkages)
    evaluationCriteria <- unique(corresponding_criterion(evaluationCriteria))

    evaluationCriteriaGrid <-
      criteria_linkages_grid(criteria = evaluationCriteria,
                             linkages = unique(grille$linkage))

    evaluationCriteriaGrid$initialValue <- NA_real_

  }

  amelioration_fun <-
    function(params) enhancer_pb(pb = pb,
                                 regionalisation = regionalisation,
                                 maxIt = maxIt,
                                 criterion = params["criterion"],
                                 linkage = params["linkage"])


  if (verbose)
  {
    alert <- gettext("Evaluation of the {nrow(grille)} enhancement{?s}")
    cli_alert(alert)
  }

  resAmeliorant <- c3t_apply(grille, 1L, amelioration_fun,
                             export = "regionalisation",
                             exportIfAbsent = "pb")

  linkage <- NULL # ony used for CHECk warnings
  if (all(is.na(grille$linkage)))
    grille <- select(grille, -linkage)


  results <- tibble(statut = vapply(resAmeliorant, "[[",
                                    character(1L), "statut"),
                    iterations = vapply(resAmeliorant, "[[",
                                        integer(1L), "it"),
                    regionalisationOpti = lapply(resAmeliorant, "[[",
                                                 "regionalisation"))

  if (verbose)
  {
    tempsCalcul <- lapply(resAmeliorant, "[[", "tempsCalcul")
    tempsCalcul <- vapply(tempsCalcul, as.double, double(1L), units = "mins")
    results <- add_column(results, tempsCalcul_mins = tempsCalcul)
  }

  rm(resAmeliorant)

  if (length(evaluationCriteria) == 0L)
    return(list(results = bind_cols(grille, results)))


  valeurfinalCriterion <-
    matrix(NA_real_, nrow = nrow(grille), ncol = nrow(evaluationCriteriaGrid))
  colnames(valeurfinalCriterion) <- evaluationCriteriaGrid$name
  results <- bind_cols(grille, results, valeurfinalCriterion)

  rm(grille, valeurfinalCriterion)


  if (verbose)
  {
    alert <-
      gettext("Calculation of {length(evaluationCriteria)} evaluation criteria on the initial partition") # nolint: line_length_linter
    cli_alert(alert)
  }



  # Calculation of the non-evaluated criteria
  # -- Initial partition
  for (l in which(is.na(evaluationCriteriaGrid$initialValue)))
  {
    criterion <- evaluationCriteriaGrid$criterion[[l]]
    linkage <- evaluationCriteriaGrid$linkage[[l]]
    evaluationCriteriaGrid$initialValue[[l]] <-
      pb$quality_partition(regionalisation, criterion, linkage)
  }

  if (verbose)
  {
    alert <-
      gettext("Calcul of {length(evaluationCriteria)} evaluation criteria on the {nrow(results)} enhanced partition{?s}") # nolint: line_length_linter
    cli_alert(alert)
  }

  # -- Enhanced partitions
  for (critere in evaluationCriteria)
  {
    if (anyNA(results[, critere]))
    {
      masqueNA <- is.na(results[, critere])
      results[masqueNA, critere] <-
        pb$quality_partition(results$regionalisationOpti[masqueNA], critere)
    }
  }

  # Add the delta improvment for each evaluation criterion
  for (critere in evaluationCriteria)
  {

    nomColonne <- paste0("dt_", critere)
    posColonne <- which(colnames(results) == critere)
    initialValue <-
      evaluationCriteriaGrid$initialValue[
        evaluationCriteriaGrid$criterion == critere] # nolint: indentation_linter
    results  <- add_column(results,
                           !!nomColonne := # nolint
                             pull(results[, critere]) - initialValue,
                           .after = posColonne)
  }

  results <- select(results, -!!evaluationCriteriaGrid$name)

  initialValues <- evaluationCriteriaGrid$initialValue
  names(initialValues) <- evaluationCriteriaGrid$name

  list(initialValues = initialValues,
       results = results)
}
