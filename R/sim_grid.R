#' @include sim_repartition.R
#' @include sim_context.R
#' @include sim_contiguity.R

#' Converts the coordinates (i,j) of a square in a grid into a position/rank
#' within the contiguity matrix.
#'
#' @param i The row index of the square.
#' @param j The column index of the square.
#' @param x_int The number of rows in the grid.
#' @param y_int The number of columns in the grid.
#'
#' @returns The position/rank of the square in the contiguity matrix.
#' @noRd
#' @importFrom checkmate assertCount assertIntegerish
xy_to_rank_grid <- function(i, j, x_int, y_int)
{
  # Argument verification
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)
  assertIntegerish(i, lower = 1L, upper = x_int)
  assertIntegerish(j, lower = 1L, upper = y_int)

  return((j - 1L) * x_int + i)
}

#' Converts a position/rank in a contiguity matrix of a grid into a pair
#' (i,j) of coordinates corresponding to the position of the element
#' in the grid.
#'
#' @param k_int The position/rank in the contiguity matrix.
#' @inheritParams xy_to_rank_grid
#'
#' @returns A list with elements 'i' and 'j', representing the row and
#' column coordinates respectively.
#' @noRd
#' @importFrom checkmate assertCount assertInt
rank_to_xy_grid <- function(k_int, x_int, y_int)
{
  # Argument verification
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)
  assertInt(k_int, lower = 1L, upper = x_int * y_int)

  j <- as.integer(ceiling(k_int / x_int))
  res <- k_int %% x_int
  i <- ifelse(res > 0L, res, x_int)

  return(list(i = i, j = j))
}


# nocov start
#' Grid Simulation
#'
#' Simulates a territory in the form of a grid with its contiguity,
#' a distribution of individuals across the cells, and a context for
#' different sizes.
#' @inheritParams grid_contiguity_matrix
#' @inheritParams gen_context
#' @inheritParams gen_repartition
#' @param seed indicates if the random seed should be fixed before generating
#' values. `NULL` (default) if the seed must not be fixed. Otherwise must
#' be an integer. (integer)
#' @returns a list including:
#' 1. `contiguity`: grid's contiguity matrix
#' 2. `distribution`: a dataframe specifying the distribution of sizes
#' across different cells. Contains:
#'    * `nbIndividuals`: number of individuals in the zone. The sum equals
#'    `nbElems_int`.
#'    * `emptyFixedZone`: `TRUE` if the zone has been fixed as empty,
#'    `FALSE` otherwise.
#'    * `metropolis`: `TRUE` if the zone is a metropolis, `FALSE` otherwise.
#' 3. `context`: a dataframe where each row corresponds to the context/data
#' related to an individual, and each column is a variable.
#' @name gen_grid
#' @keywords internal
#' @export
gen_grid <- function(x_int, y_int, contiguityType = "Queen",
                     nbIndividuals = x_int * y_int,
                     nbMinEmptyZones = 0L,
                     nbMetropolises = 0L,
                     propInMetropolises = 0.4,
                     quantitatives_mat = c(0.0, 1.0),
                     qualitatives_list = NULL,
                     nbQuantitatives = 0L,
                     seed = NULL)
{
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)
  assertInt(seed, null.ok = TRUE)

  if (!is.null(seed))
    set.seed(seed)

  nbCases <- x_int * y_int

  matCont <- grid_contiguity_matrix(x_int, y_int, contiguityType = "Queen")

  repartition_df <-
    gen_repartition(nbElements_int = nbCases,
                    nbIndividuals_int = nbIndividuals,
                    nbMinEmptyZones_int = nbMinEmptyZones,
                    nbMetropolises_int = nbMetropolises)

  context_df <- gen_context(nbCases, quantitatives_mat = quantitatives_mat,
                            qualitatives_list = qualitatives_list,
                            nbQuantitatives_int = nbQuantitatives)

  return(list(contiguity = matCont, repartition = repartition_df,
              context = context_df))
}
# nocov end

# nocov start
#' @describeIn gen_grid simpler version with "Queen" contiguity,
# no mandatory empty zones, and no metropolis, and for context 3
# quantitative variables with different means and variances
#' @keywords internal
#' @importFrom checkmate assertCount assertNumber
#' @examples simple_grid(2L, 3L, 100.0, 123L)
#' @export
simple_grid <- function(x_int, y_int = x_int,
                        avgPersonsPerCell = 100.0,
                        seed = NULL)
{
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)
  assertNumber(avgPersonsPerCell,
               lower = 0.0, finite = TRUE,
               na.ok = FALSE, null.ok = FALSE)


  nbCells <- x_int * y_int
  nbIndividuals <- ceiling(nbCells * avgPersonsPerCell)

  quantitatives_mat <- matrix(c(  0.0, 1.0, # nolint
                                 10.0, 5.0, # nolint
                                -10.0, 9.0), byrow = TRUE, ncol = 2L)

  colnames(quantitatives_mat) <- c("mean", "sd")

  gen_grid(x_int, y_int,
           contiguityType = "Queen",
           nbIndividuals = nbIndividuals,
           nbMinEmptyZones = 0L,
           nbMetropolises = 0L,
           quantitatives_mat = quantitatives_mat,
           qualitatives_list = NULL,
           seed = seed)
}
# nocov end

gen_pb_from_grid <- function(grid, m = 0.0, M = Inf,
                             d = "euclidean",
                             standardQuant = TRUE,
                             binarQual = TRUE,
                             storageMode = "matrix",
                             calculateAllDistances = TRUE,
                             p = 2.0)
{
  data <- normalize_df(grid$context, standardQuant, binarQual)
  d <- assertElementDistance(d, p)
  d <- corresponding_d(d)

  pb <- constructor_pbCon(data = data,
                          sizes = grid$repartition$nbIndividuals,
                          m = m, M = M,
                          d = d,
                          standardQuant = standardQuant,
                          binarQual = binarQual,
                          storageMode = storageMode)

  if (calculateAllDistances)
    pb[]

  pb

}

gen_pb <- function(x_int, y_int, m = 0.0, M = Inf,
                   d = "euclidean",
                   standardQuant = TRUE,
                   binarQual = TRUE,
                   storageMode = "matrix",
                   calculateAllDistances = TRUE,
                   p = NULL,
                   seed = 123L)
{
  grid <- simple_grid(x_int, y_int, 100.0, seed)

  gen_pb_from_grid(grid,
                   m, M,
                   d,
                   standardQuant, binarQual,
                   storageMode,
                   calculateAllDistances, p)

}

# nocov start
survey_simulation <- function(x_int, y_int, nbIndividuals_int)
{
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)
  assertCount(nbIndividuals_int, positive = TRUE)

  nbCases <- x_int * y_int
  mat_cont <- grid_contiguity_matrix(x_int, y_int, contiguityType = "Queen")

  repartition_df <- gen_repartition(nbCases, nbIndividuals_int,
                                    nbMinEmptyZones_int =
                                      ceiling(nbCases * 29.0 / 35.0),
                                    nbMetropolises_int =
                                      ceiling(0.01 * nbCases))

  context_df <- gen_context(nbCases,
                            quantitatives_mat =
                              matrix(c(0.0, 1.0,
                                       10.0, 5.0), byrow = TRUE, ncol = 2L))

  context_df[repartition_df$metropolis, ] <-
    context_df[repartition_df$metropolis, ] + 3.0

  context_df$type <- "urbain"
  context_df[repartition_df$nbIndividuals == 0L, "type"] <- "rural"

  return(list(contiguity = mat_cont, repartition = repartition_df,
              context = context_df))

}
# nocov end


#' Simulation of a survey
#'
#' This grid simulates a survey. It has a dimension of 50 * 70 (3500 cells).
#' Most of the cells are empty (around 80%). The total number of people is
#' 4000. 1% of the cells are metropolises, aggregating 40% of the population.
#' The context consists of 2 quantitative and 1 qualitative variable.
#' @docType data
#' @keywords data
#' @name survey
#' @rdname survey
#' @keywords internal
#' @format A list with the following components:
#' 1. `contiguity`: Contiguity matrix of the grid
#' 2. `repartition`: A dataframe specifying the distribution of sizes
#' across different cells. Contains:
#'    * `nbIndividuals`: number of individuals in the zone. The sum is
#'    `nbElems_int`.
#'    * `emptyFixedZone`: `TRUE` if the zone has been fixed as empty,
#'    `FALSE` otherwise.
#'    * `metropolis`: `TRUE` if the zone is a metropolis, `FALSE` otherwise.
#' 3. `context`: A dataframe where each row corresponds to the context /
#' data related to an individual and each column is a variable.
"survey_grid_50x70"

# nolint start
# set.seed(123L)
# survey_grid_50x70 = survey_simulation(50, 70, 4000)
# use_data(survey_grid_50x70, overwrite = TRUE)
# nolint end
