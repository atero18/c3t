#' @param nbIndividuals Number of individuals to distribute
#' (strictly positive integer).
#' @param nbMinEmptyZones Number of zones where no individuals should
#' be present. Default is 0. `nbMinEmptyZones + nbMetropolises` must be smaller
#'  than `x_int * y_int`.
#' (positive integer)
#' @param nbMetropolises Number of zones that should be more populated.
#' Default is 0. `nbMinEmptyZones_int + nbMetropolises`_int must be smaller than
#' `x_int * y_int`. (positive integer)
#' @param propInMetropolises Proportion of individuals found in metropolises
#' (if any). (floating point in `[0,1]`)
#' @name gen_repartition
NULL

# Simulates a distribution of sizes across different elements / zones of a set.
# It is possible to provide specific characteristics to a number of elements.
#' @importFrom checkmate assertCount assertNumber
gen_repartition <- function(nbElements_int,
                            nbIndividuals_int,
                            nbMinEmptyZones_int = 0L,
                            nbMetropolises_int = 0L,
                            ppt_in_metropolises_real = 0.4)
{
  # Checking arguments
  assertCount(nbElements_int, positive = TRUE)
  assertCount(nbIndividuals_int)

  if (nbIndividuals_int == 0L)
  {
    return(data.frame(nbIndividuals = integer(nbElements_int),
                      emptyFixedZone = logical(nbElements_int),
                      metropolis = logical(nbElements_int)))
  }
  assertCount(nbMinEmptyZones_int)
  assertCount(nbMetropolises_int)
  if (nbMetropolises_int > 0L)
    assertNumber(ppt_in_metropolises_real, lower = 0.0, upper = 1.0)

  if (nbMinEmptyZones_int + nbMetropolises_int > nbElements_int)
  {
    stop("There are not enough available zones to consider all parameters")
  }

  df_repartition <- data.frame(nbIndividuals = integer(nbElements_int),
                               emptyFixedZone = logical(nbElements_int))
  specialZones_vec <- NULL
  # Creating empty zones
  if (nbMinEmptyZones_int > 0L)
  {
    emptyZones_vec <- sample(seq_len(nbElements_int),
                             nbMinEmptyZones_int, FALSE)
    specialZones_vec <- union(specialZones_vec, emptyZones_vec)
    df_repartition$emptyFixedZone[emptyZones_vec] <- TRUE
  }
  else
    emptyZones_vec <- NULL


  df_repartition$metropolis <- FALSE

  # Creating metropolises
  if (nbMetropolises_int > 0L)
  {
    nbIndividualsMetropolis <- ceiling(ppt_in_metropolises_real *
                                    nbIndividuals_int)
    nbIndividualsOutsideMetropolis <-
      nbIndividuals_int - nbIndividualsMetropolis
    metropolisZones_vec <-
      sample(setdiff(seq_len(nbElements_int), emptyZones_vec),
             nbMetropolises_int, FALSE)

    specialZones_vec <- union(specialZones_vec, metropolisZones_vec)

    df_repartition$metropolis[metropolisZones_vec] <- TRUE
    rep_inhabitants_metropolis <- sample(seq_len(nbMetropolises_int),
                                         nbIndividualsMetropolis,
                                         TRUE)

    df_repartition[metropolisZones_vec, "nbIndividuals"] <-
      as.vector(table(rep_inhabitants_metropolis))
  }
  else
  {
    nbIndividualsOutsideMetropolis <- nbIndividuals_int
    emptyZones_vec <- NULL
  }

  # Managing residents who have not yet been placed
  if (nbIndividualsOutsideMetropolis > 0L)
  {
    normalZones_vec <- setdiff(seq_len(nbElements_int), specialZones_vec)
    repInhOutMetropolis <- sample(normalZones_vec,
                                  nbIndividualsOutsideMetropolis,
                                  TRUE)
    repInhOutMetropolis <-
      table(repInhOutMetropolis)

    df_repartition[names(repInhOutMetropolis),
                   "nbIndividuals"] <-
      as.vector(repInhOutMetropolis)
  }

  return(df_repartition)
}

#' @importFrom checkmate assertCount assertChoice
one_person_per_case <- function(n, type = "double")
{
  assertCount(n)
  assertChoice(type, c("double", "integer"))
  value <- ifelse(type == "double", 1.0, 1L)
  rep(value, n)
}
