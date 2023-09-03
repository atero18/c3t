envirClustersc3t <- new.env()
envirClustersc3t$c3tClusters <- NULL
options(c3t_parallel = FALSE) # nolint
options(c3t_parallelizing = FALSE) # nolint
options(c3t_verbose_parallel = FALSE) # nolint


# Create and remove c3t clusters ------------------------------------------


c3t_get_clusters <- function()
{
  get("c3tClusters", envir = envirClustersc3t)
}

c3t_clusters_exist <- function()
{
  !is.null(c3t_get_clusters())
}

#' @importFrom rlang env_has
#' @importFrom cli cli_alert cli_alert_warning cli_alert_info cli_alert_danger
doInPar <- function()
{
  isTRUE(getOption("c3t_parallel", default = FALSE)) &&
    c3t_clusters_exist() &&
    isFALSE(getOption("c3t_parallelizing", default = TRUE))
}

#' @param parallele Logical indicating whether to use parallel processing.
#'  Default is TRUE.
#' @param nbCores Number of CPU cores to use for parallel processing
#' (sockets method). Default is one less than the detected number of cores.
#' @name parallel_arguments
NULL

#' @importFrom parallel detectCores makeCluster clusterEvalQ
#' @importFrom rlang env_has
c3t_create_clusters <- function(nbCores = detectCores() - 1L,
                                verbose = getOption("c3t_verbose",
                                                    default = FALSE))
{

  assertCount(nbCores, positive = TRUE)
  options(c3t_parallel = FALSE) # nolint: undesirable_function_linter

  if (c3t_clusters_exist())
  {
    # Case where clusters are currently in use
    if (isTRUE(getOption("c3t_parallelizing", default = TRUE)))
    {
      WARNING <- gettext("Clusters already created and in use")
      cli_alert_warning(WARNING)
      return(invisible(c3t_get_clusters()))
    }
    # If the clusters already exist and are of the desired number,
    # we clear them
    else if (c3t_nbCores() == nbCores)
    {
      INFO <- gettext("Existing clusters. Clearing their data")
      cli_alert_info(INFO)
      clusterEvalQ(c3t_get_clusters(), rm(list = ls()))
      clusterEvalQ(c3t_get_clusters(), library(c3t)) # nolint: undesirable_function_linter unused_imoprt_linter
      options(c3t_parallel = FALSE) # nolint: undesirable_function_linter
      return(invisible(c3t_get_clusters()))
    }
    else
    {
      WARNING <- gettext("Existing clusters. Complete reset")
      cli_alert_warning(WARNING)
      c3t_stop_clusters()
    }
  }

  options(c3t_parallelizing = FALSE) # nolint: undesirable_function_linter

  if (nbCores > detectCores())
  {
    warning <-
      gettext("The number of cores proposed ({nbCores}) is greater than the number of machine cores ({detectCores()}). Reducing") # nolint: line_length_linter
    cli_alert_warning(warning)
    nbCores <- detectCores()
  }

  if (nbCores == 1L)
  {
    WARNING <- gettext("The requested number of cores is 1. Unnecessary")
    cli_alert_warning(WARNING)
  }

  if (verbose)
  {
    alert <- gettext("Creating {nbCores} clusters")
    cli_alert(alert)
  }

  c3tClusters <- makeCluster(nbCores)
  assign("c3tClusters", c3tClusters, envir = envirClustersc3t)
  clusterEvalQ(c3tClusters, library(c3t)) # nolint: undesirable_function_linter
  clusterEvalQ(c3tClusters,
               options(c3t_parallelizing = TRUE)) # nolint: undesirable_function_linter
  options(c3t_parallel = TRUE) # nolint: undesirable_function_linter

  invisible(c3tClusters)
}


#' @importFrom parallel stopCluster
#' @importFrom cli cli_alert_info
c3t_stop_clusters <- function(forceStop = FALSE)
{
  assertFlag(forceStop)

  if (c3t_clusters_exist())
  {
    if (isTRUE(getOption("c3t_parallelizing", default = TRUE)))
    {
      if (!forceStop)
      {
        INFO <-
          gettext("Clusters are present but in use. Set `forceStop` to `TRUE` to force stopping") # nolint: line_length_linter
        cli_alert_info(INFO)
        return(invisible(FALSE))
      }
      else
      {
        INFO <- gettext("Clusters are present but in use. Forceful deletion")
        cli_alert_info(INFO)
      }
    }

    options(c3t_parallel = FALSE) # nolint: undesirable_function_linter
    options(c3t_parallelizing = FALSE) # nolint: undesirable_function_linter

    stopCluster(c3t_get_clusters())

    assign("c3tClusters", NULL, envir = envirClustersc3t)
  }
  invisible(TRUE)
}



c3t_nbCores <- function()
{
  ifelse(c3t_clusters_exist(), length(c3t_get_clusters()), 0L)
}



# Send and remove data on c3t clusters ------------------------------------


#' @importFrom parallel clusterEvalQ
#' @importFrom rlang inject
c3t_clusterEvalQ <- function(commande)
{
  invisible(inject(clusterEvalQ(c3t_get_clusters(), !!commande)))
}

#' @importFrom parallel clusterExport clusterEvalQ
c3t_clusterExport <- function(export  = character(0L),
                              exportIfAbsent = character(0L),
                              envirExport = parent.frame(n = 2L))
{
  c3tClusters <- c3t_get_clusters()
  if (length(export) > 0L)
    clusterExport(c3tClusters, export, envir = envirExport)



  if (length(exportIfAbsent) > 0L)
  {
    clusterExport(c3tClusters, "exportIfAbsent", envir = environment())
    presence <- clusterEvalQ(c3tClusters, vapply(exportIfAbsent,
                                                 exists,
                                                 logical(1L)))

    clusterEvalQ(c3tClusters, rm(exportIfAbsent))
    presence <- do.call("rbind", presence)
    elementsNonPresents <- !apply(presence, 2L, all)

    if (any(elementsNonPresents))
      clusterExport(c3tClusters, exportIfAbsent[elementsNonPresents],
                    envir = envirExport)
  }

  return(TRUE)
}

#' @importFrom parallel clusterExport clusterEvalQ
c3t_clusters_rm <- function(elementsToRemove = character())
{
  if (c3t_clusters_exist() && length(elementsToRemove) > 0L)
  {
    c3tClusters <- c3t_get_clusters()
    clusterExport(c3tClusters, "elementsToRemove", envir = environment())
    clusterEvalQ(c3tClusters, rm(list = elementsToRemove))
    clusterEvalQ(c3tClusters,  rm(elementsToRemove))
  }

  TRUE
}


# apply functions for c3t (sequential / parallel) -------------------------

#' @importFrom cli cli_alert_info cli_alert_danger
#' @importFrom checkmate assertFlag
c3t_APP <- function(fonctionApply, fonctionParApply,
                    parallele,
                    donnees, FUN,
                    export = character(0L),
                    exportIfAbsent = character(0L),
                    elementsToRemove = export,
                    envirExport)
{
  assertFlag(parallele)

  parallele <- parallele && !is.null(fonctionParApply)
  if (parallele)
  {
    if (getOption("c3t_verbose_parallel", default = FALSE))
    {
      info <- gettext("Parallel computation using {c3t_nbCores()} clusters")
      cli_alert_info(info)
    }

    options(c3t_parallelizing = TRUE) # nolint
    on.exit(options(c3t_parallelizing = FALSE), # nolint
            add = TRUE,
            after = TRUE)

    erreurRealisee <- FALSE
    res <- tryCatch(
      expr =
      {
        c3t_clusterExport(export,
                          exportIfAbsent,
                          envirExport = envirExport)
        fonctionParApply(donnees, FUN)
      },
      error = function(cond)
      {
        DANGER <-
          gettext("Parallelization did not work. Attempting sequential computation") # nolint: line_length_linter
        cli_alert_danger(DANGER)
        erreurRealisee <<- TRUE # nolint: undesirable_operator_linter
      }
    )

    if (!erreurRealisee)
    {
      c3t_clusters_rm(elementsToRemove)
      return(res)
    }
    else
      parallele <- FALSE
  }

  if (!parallele)
    return(fonctionApply(donnees, FUN))
}

#' @importFrom parallel parLapply
#' @importFrom purrr partial

c3t_lapply <- function(donnees, FUN, ...,
                       export = character(0L),
                       exportIfAbsent = character(0L),
                       elementsToRemove = export,
                       envirExport = parent.frame(1L))
{
  parallele <-  doInPar() && length(donnees) > 1L

  if (parallele)
  {
    c3t_fonct_parLapply <-
      function(donnees, fun) parLapply(cl = c3t_get_clusters(),
                                       X = donnees, fun = FUN, ...)
  }
  else
    c3t_fonct_parLapply <- NULL

  c3t_fonc_lapply <- function(donnees, FUN)  lapply(X = donnees, FUN = FUN, ...)

  c3t_APP(c3t_fonc_lapply, c3t_fonct_parLapply,
          parallele,
          donnees, FUN,
          export,
          exportIfAbsent,
          elementsToRemove,
          envirExport)
}

#' @importFrom parallel parSapply
#' @importFrom purrr partial
c3t_sapply <- function(donnees, FUN, ..., simplify = TRUE,
                       export = character(0L),
                       exportIfAbsent = character(0L),
                       elementsToRemove = export,
                       envirExport = parent.frame(1L))
{
  parallele <- doInPar() && length(donnees) > 1L

  if (parallele)
  {
    c3t_fonct_parSapply <-
      function(donnees, FUN) parSapply(cl = c3t_get_clusters(),
                                       X = donnees,
                                       FUN = FUN,
                                       simplify = simplify,
                                       ...) # nolint
  }
  else
    c3t_fonct_parSapply <- NULL

  c3t_fonct_sapply <-
    function(donnees, FUN) sapply(X = donnees, # nolint: undesirable_function_linter
                                  FUN = FUN,
                                  simplify = simplify, ...)

  c3t_APP(c3t_fonct_sapply,
          c3t_fonct_parSapply,
          parallele,
          donnees, FUN,
          export,
          exportIfAbsent,
          elementsToRemove,
          envirExport)
}

#' @importFrom parallel parApply
#' @importFrom purrr partial
c3t_apply <- function(donnees, MARGIN, FUN, ...,
                      export = character(0L),
                      exportIfAbsent = character(0L),
                      elementsToRemove = export,
                      envirExport = parent.frame(1L))
{
  parallele <- doInPar()

  if (parallele && MARGIN == 1L && ncol(donnees) == 1L)
    parallele <- FALSE

  if (parallele && MARGIN == 2L && nrow(donnees) == 1L)
    parallele <- FALSE


  if (parallele)
  {
    c3t_fonct_parApply <-
      function(donnees, FUN) parApply(c3t_get_clusters(),
                                      X = donnees,
                                      FUN = FUN,
                                      MARGIN = MARGIN,
                                      ...)
  }
  else
    c3t_fonct_parApply <- NULL

  c3t_fonct_Apply <-
    function(donnees, FUN) apply(X = donnees,
                                 FUN = FUN,
                                 MARGIN = MARGIN,
                                 ...)

  c3t_APP(c3t_fonct_Apply,
          c3t_fonct_parApply,
          parallele,
          donnees, FUN,
          export,
          exportIfAbsent,
          elementsToRemove,
          envirExport)
}
