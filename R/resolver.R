#' Resolve an unfeasible regionalisation
#'
#' This algorithm tries to transform an unfeasible solution that check
#' connectivity constraint but not one or more size constraints. Even if
#' the problem is feasible the algorithm doesn't ensure to succeed.
#' `r badge('experimental')`
#'
#' @inheritParams arguments_problem
#' @inheritParams verbose_argument
#' @param regionalisation A partition checking the contiguity constraint but
#' not checking at least the min or the max size constraint.
#' @param maxItTransfers positive integer indicating the maximum number
#' of elements can be transferred.
#' @param ... used for development.
#' @example inst/examples/resolve_unfeasible.R
#' @export
#' @importFrom checkmate assertFlag
resolve_unfeasible <- function(distances = NULL, contiguity = NULL,
                               sizes = NULL,
                               d = NULL, data = NULL, m = 0.0, M = Inf,
                               standardQuant = FALSE, binarQual = FALSE,
                               regionalisation, # nolint : function_argument_linter
                               maxItTransfers = Inf,
                               verbose = FALSE,
                               ...)
{
  assertFlag(verbose)
  options(c3t_verbose = verbose) # nolint

  pb <- constructor_pbCon(distances = distances,
                          contiguity = contiguity,
                          sizes = sizes,
                          d = d, data = data,
                          m = m, M = M,
                          standardQuant = standardQuant,
                          binarQual = binarQual)


  resSolve <- resolve_regionalisation_pb(pb, regionalisation,
                                         allowTransferts = TRUE,
                                         authorizedFusions = NA_character_,
                                         maxItTransfers = maxItTransfers,
                                         maxItFusions = 0L,
                                         ...)


  resSolve[!names(resSolve) %in% c("relaxationMin", "relaxationMax")]

}


#' @importFrom cli cli_alert cli_alert_danger
#' @importFrom cli cli_alert_warning cli_alert_success
#' @importFrom checkmate assertFlag assertNumber assertString
resolve_regionalisation_pb <- # nolint: cyclocomp_linter
  function(pb, regionalisation,
           relaxationMin = 1.0, relaxationMax = 1.0,
           authorizedFusions = NA_character_,
           allowSingletonTransfert = TRUE,
           maxItTransfers = Inf,
           maxItFusions = Inf,
           allowTransferts = TRUE)
{

  verbose <- getOption("c3t_verbose", default = FALSE)

  # Checking arguments
  validObject(pb)
  assertNumber(relaxationMin, lower = 0.0, upper = 1.0)
  assertNumber(relaxationMax, lower = 1.0, finite = FALSE)

  assertString(authorizedFusions, na.ok = TRUE)

  assertCountOrInfinite(maxItFusions)

  if (maxItFusions == 0L)
    allowFusions <- FALSE
  else
    allowFusions <- !is.na(authorizedFusions)

  assertFlag(allowTransferts)



  if (allowTransferts)
  {
    assertCountOrInfinite(maxItTransfers)
    if (maxItTransfers == 0L)
      allowTransferts <- FALSE
  }

  if (allowTransferts)
  {
    if (allowFusions)
      allowSingletonTransfert <- TRUE
    else
      assertFlag(allowSingletonTransfert)
  }
  else
    allowSingletonTransfert <- FALSE

  if (!allowTransferts && !allowFusions)
  {
    stop("Transfert or fusion must be allowed")
  }


  # Given partition must check contiguity partition
  if (!is_regionalisation(regionalisation, pb$contiguity))
  {
    stop("`regionalisation` doesn't check connectivity constraint")
  }

  regionalisation <- standardize_partition(regionalisation)

  # Application of size constraint relaxation
  m <- pb$m
  if (pb$hasMinConstraint())
    pb$m <- relaxationMin * m

  M <- pb$M
  if (pb$hasMaxConstraint())
    pb$M <- relaxationMax * M

  # Regions sizes
  initialRegionsSizes <- regionsSizes <-
    .clusters_sizes(regionalisation, pb$sizes)

  regionsCardinals <- table(regionalisation)
  nbClustersInit <- length(regionsCardinals)

  # Regions are sorted depending of their respect of the size constraints
  if (pb$hasMinConstraint())
    smallRegions <- initialSmallRegions <- which(regionsSizes < pb$m)

  else
    smallRegions <- initialSmallRegions <- integer(0L)


  if (length(initialSmallRegions) == 0L)
    initialSmallRegions <- NULL

  if (pb$hasMaxConstraint())
    bigRegions <- initialBigRegions <- which(regionsSizes > pb$M)

  else
    bigRegions <- initialBigRegions <- integer(0L)

  if (length(initialBigRegions) == 0L)
    initialBigRegions <- NULL

  status <- NULL

  # If all regions check all size constraints then there is nothing to
  # resolve
  if (length(smallRegions) + length(bigRegions) == 0L)
    status <- "already_feasible"

  resolvedRegions <- integer(0L)
  n <- pb$n()

  # We randomly take a permutation of the elements set
  listeElements_vec <- sample(seq_len(n))

  dataElementsTransferables <- NULL
  oldRegionalisation <- regionalisation

  itTransfers <- 0L
  calcul_elementsTransferables <- function()
  {

    # Calculation of update of the transferable elements set
    if (itTransfers == 0L)
    {
      assign("dataElementsTransferables",
             transferable_elements_pb(regionalisation, pb = pb),
             envir = parent.frame())
    }

    else
    {
      assign("dataElementsTransferables",
             maj_elements_transferables(oldRegionalisation,
                                        regionalisation,
                                        dataElementsTransferables,
                                        donor, receiver,
                                        pb = pb),
             envir = parent.frame())
    }

    # We only consider transitions that involve:
    # - a region r1 that is too small (which would receive) and a region r2
    #   of valid or too large size (which would give the element), such that
    #   if r2 gives this element, the region wouldn't be too small, and r1
    #   wouldn't be too large.
    # - a region r1 that is too large (which would give) and a region r2 of
    #   valid size (which would receive), such that if r1 gives this element
    #   to r2, r2 wouldn't become too large.

    # Specifically, regions that are too small cannot give, and those that are
    # too large cannot receive.

    assign("dataElementsTransferables",
           remove_transferables(pb,
                                regionalisation,
                                dataElementsTransferables,
                                regionsSizes,
                                rmFromSingleton = !allowSingletonTransfert,
                                rmFromMin = TRUE,
                                rmToMin = FALSE,
                                rmFromMax = FALSE,
                                rmToMax = TRUE,
                                rmCreationMin = TRUE,
                                rmCreationMax = TRUE,
                                rmNormToNorm = TRUE,
                                smallClusters = smallRegions,
                                bigClusters = bigRegions,
                                nbElemsClusters = regionsCardinals),
           envir = parent.frame(1L))


    dataElementsTransferables$transferable
  }

  elementsTransferables <- calcul_elementsTransferables()

  itFusions <- 0L
  # Gestion de la fusion des clusters
  if (allowFusions && length(smallRegions) > 0L)
  {
    # nocov start
    if (verbose)
    {
      alert <- gettext("Clusters fusion")
      cli_alert(alert)

    }
    # nocov end

    fusion <- FALSE
    clustersVoisins <-
      unique(dataElementsTransferables$frontier[, c("cluster1", "cluster2")],
             MARGIN = 1L)

    clustersVoisins <-
      clustersVoisins[clustersVoisins[, 1L] < clustersVoisins[, 2L], ,
                      drop = FALSE]

    # Suppression des clusters trop gros
    if (nrow(clustersVoisins) > 0L)
    {
      suppressionMax <- clustersVoisins[, "cluster1"] %in% bigRegions |
        clustersVoisins[, "cluster2"] %in% bigRegions
      clustersVoisins <- clustersVoisins[!suppressionMax, , drop = FALSE]
    }

    while (itFusions <= maxItFusions && nrow(clustersVoisins) > 0L &&
          length(smallRegions) > 0L)
    {
      fusion <- FALSE

      auMoinsUnMin <- clustersVoisins[, "cluster1"] %in% smallRegions |
        clustersVoisins[, "cluster2"] %in% smallRegions
      clustersVoisins <- clustersVoisins[auMoinsUnMin, , drop = FALSE]

      if (authorizedFusions != "all")
      {
        suppressionNonMins <-
          clustersVoisins[, "cluster1"] %in% smallRegions &
          clustersVoisins[, "cluster2"] %in% smallRegions
        clustersVoisins <- clustersVoisins[!suppressionNonMins, , drop = FALSE]
      }

      # Checking size after fusion
      if (nrow(clustersVoisins) > 0L && pb$hasMaxConstraint())
      {
        masqueTaille <- regionsSizes[clustersVoisins[, "cluster1"]] +
                        regionsSizes[clustersVoisins[, "cluster2"]] <= pb$M
        clustersVoisins <- clustersVoisins[masqueTaille, , drop = FALSE]
      }

      if (nrow(clustersVoisins) > 0L)
      {
        itFusions <- itFusions + 1L
        # nocov start
        if (verbose && (itFusions == 1L || itFusions %% 10L == 0L))
        {
          alert <- gettext("Fusion {itFusions} (max: {maxItFusions})")
          cli_alert(alert)
        }
        # nocov end



        # Fusion between two mins have priority
        couplesMins <-
          which(clustersVoisins[, "cluster1"] %in% smallRegions &
                clustersVoisins[, "cluster2"] %in% smallRegions)

        if (length(couplesMins) > 0L)
          k <- sample(couplesMins, size = 1L)
        else
          k <- sample(seq_len(nrow(clustersVoisins)), size = 1L)

        cluster1 <- clustersVoisins[k, "cluster1"]
        cluster2 <- clustersVoisins[k, "cluster2"]
        regionalisation[regionalisation == cluster2] <- cluster1

        # Contiguities update
        clustersVoisins <- update_contiguity_list(clustersVoisins,
                                                  cluster1,
                                                  cluster2)
        # Clusters parameters update
        regionsSizes[cluster1] <- regionsSizes[cluster1] +
                                  regionsSizes[cluster2]

        regionsSizes[cluster2] <- NaN

        if (cluster2 %in% smallRegions)
        {
          smallRegions <- smallRegions[smallRegions != cluster2]
          resolvedRegions <- c(smallRegions, cluster2)
        }

        # If the new cluster now checks the minimum size constraints
        if (cluster1 %in% smallRegions &&
            regionsSizes[cluster1] >= pb$m)
        {
          smallRegions <-
            smallRegions[smallRegions != cluster1]

          resolvedRegions <- c(smallRegions, cluster1)
        }

        regionsCardinals[cluster1] <- regionsCardinals[cluster1] +
                                      regionsCardinals[cluster2]

        regionsCardinals[cluster2] <- NaN
      }
    }
    if (fusion)
      elementsTransferables <- calcul_elementsTransferables()
  }

  # nocov start
  if (verbose && itFusions > 0L)
  {
    alert <- gettext("Number of iterations (fusions) : {itFusions}")
    cli_alert(alert)
  }
  # nocov end

  # nocov start
  # Transferts are made while there is at least one region which doesn't
  # check size constraints and there exists transferable elements
  if (allowTransferts && verbose)
  {
    ALERT <- gettext("Transfert of elements one-by-one")
    cli_alert(ALERT)
  }
  # nocov end

  while (allowTransferts &&
        length(smallRegions) + length(bigRegions) > 0L &&
        length(elementsTransferables) > 0L)
  {
    itTransfers <- itTransfers + 1L

    # We take the first element in the permuted set that is transferable
    x1 <- listeElements_vec[listeElements_vec %in% elementsTransferables][1L]

    # If it can be transferred with multiple elements we chose the transfert
    # randomly
    line <- which(dataElementsTransferables$frontier[, "x1"] == x1)
    if (length(line) > 1L)
      line <- sample(line, size = 1L)

    donor <- regionalisation[x1]
    receiver <- unname(dataElementsTransferables$frontier[line, "cluster2"])


    # Update of the partition
    oldRegionalisation <- regionalisation
    regionalisation[x1] <- receiver

    # Regions sizes and their status regarding size constraints are updated
    regionsSizes[donor] <- regionsSizes[donor] - pb$sizes[x1]
    if (donor %in% bigRegions &&
        regionsSizes[donor] <= pb$M)
    {
      bigRegions <- setdiff(bigRegions, donor)
      resolvedRegions <- c(resolvedRegions, donor)
    }
    regionsCardinals[donor] <- regionsCardinals[donor] - 1L

    if (donor %in% smallRegions &&
        regionsCardinals[donor] == 0L)
    {
      smallRegions <- setdiff(smallRegions, donor)
      resolvedRegions <- c(resolvedRegions, donor)
    }

    regionsSizes[receiver] <- regionsSizes[receiver] + pb$sizes[x1]
    if (receiver %in% smallRegions &&
        regionsSizes[receiver] >= pb$m)
    {
      smallRegions <- setdiff(smallRegions, receiver)
      resolvedRegions <- c(resolvedRegions, receiver)
    }
    regionsCardinals[receiver] <- regionsCardinals[receiver] + 1L

    # Update of transferable elements set
    if (length(smallRegions) + length(bigRegions) > 0L)
      elementsTransferables <- calcul_elementsTransferables()


  }

  # nocov start
  if (verbose && itFusions > 0L)
  {
    alert <- gettext("Number of iterations (transferts) : {itTransfers}")
    cli_alert(alert)
  }
  # nocov end


  if (is.null(status))
  {
    if (itTransfers + itFusions == 0L)
      status <- "unresolvable"

    else if (length(smallRegions) + length(bigRegions) > 0L)
      status <- "partially_resolved"

    else
      status <- "fully_resolved"
  }

  # nocov start
  if (verbose)
  {
    message <-
      switch(status,
             unresolvable = gettext("unresolvable partition"),
             partially_resolved = gettext("partially resolved partition"),
             fully_resolved = gettext("fully resolved partition"))


    switch(status,
           unresolvable = cli_alert_danger(message),
           partially_resolved = cli_alert_warning(message),
           fully_resolved = cli_alert_success(message)
    )
  }
  # nocov end

  pb$m <- m
  pb$M <- M

  if (length(smallRegions) == 0L)
    smallRegions <- NULL

  if (length(bigRegions) == 0L)
    bigRegions <- NULL

  if (length(resolvedRegions) == 0L)
    resolvedRegions <- NULL

  list(status = status,
       itTransfers = itTransfers,
       itFusions = itFusions,
       initialSmallRegions = unname(initialSmallRegions),
       finalSmallRegions = unname(smallRegions),
       initialBigRegions = unname(initialBigRegions),
       finalBigRegions = unname(bigRegions),
       resolvedRegions = resolvedRegions,
       initialRegionsSizes = initialRegionsSizes,
       finalRegionsSizes = regionsSizes,
       regionalisation = standardize_partition(regionalisation),
       initialNbClusters = nbClustersInit,
       finalNbClusters = length(unique(regionalisation)),
       relaxationMin = relaxationMin,
       relaxationMax = relaxationMax)
}
