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

# First simulation with a simple grid where we have "Queen" contiguity,
# no mandatory empty zones, and no metropolis, and for context 3
# quantitative variables with different means and variances
#' @importFrom checkmate assertCount assertDouble
grid_simulation_1 <- function(x_int,
                              y_int,
                              avgPersonsPerCell_int = 100.0)
{
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)
  assertDouble(avgPersonsPerCell_int, len = 1L, any.missing = FALSE)

  if (avgPersonsPerCell_int <= 0.0)
  {
    stop("avgPersonsPerCell_int must be strictly positive")
  }

  nbCells <- x_int * y_int
  contiguity_matrix <- grid_contiguity_matrix(x_int,
                                              y_int,
                                              contiguityType = "Queen")
  repartition_df <-
    gen_repartition(nbCells,
                    nbCells * avgPersonsPerCell_int,
                    nbMinEmptyZones_int = 0L,
                    nbMetropolises_int = 0L)

  context_df <- gen_context(nbCells,
                            quantitatives_mat =
                              matrix(c(0.0, 1.0,
                                       10.0, 5.0,
                                       -10.0, 9.0), byrow = TRUE, ncol = 2L))

  return(list("contiguity" = contiguity_matrix,
              "repartition" = repartition_df,
              "context" = context_df))
}


#' @importFrom checkmate assertCount testString
c3t_grid_simulation <- function(x_int, y_int, m = 0.0, M = Inf,
                                distance = "euclidean",
                                standardQuant = TRUE,
                                binarQual = TRUE,
                                storageMode = "matrix",
                                calculToutesValeurs = TRUE,
                                p = 2.0)
{
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)

  nomData <- paste0("grid_queen_vide0_metropole0_x", x_int,
                    "_y", y_int,
                    "_indivMoy100_quant3_qual0")
  tryCatch(data(list = nomData, package = "c3t", envir = environment()),
           warning = function(cond) stop("grid not found"))

  grid <- get(nomData, envir = environment())
  rm(nomData)

  data <- normalize_df(grid$context, standardQuant, binarQual)

  if (is.function(distance))
    d <- distance
  else if (testString(distance) && distance %in%
           c("euclidean", "manhattan", "minkowski"))
  {
    d <- switch(distance,
                "manhattan" = distance_manhattan,
                "euclidean" = distance_euclidienne,
                "minkowski" = function_distance_minkowski(p))
  }
  else
  {
    stop("`distance` is incorrect")
  }

  pb <- constructor_pbCon(data = data,
                          sizes = grid$repartition$nbIndividuals,
                          m = m, M = M,
                          d = d,
                          standardQuant = standardQuant,
                          binarQual = binarQual,
                          storageMode = storageMode)

  if (calculToutesValeurs)
   pb[]

  pb
}

#' Grid Simulation
#'
#' Simulates a territory in the form of a grid with its contiguity,
#' a distribution of individuals across the cells, and a context for
#' different sizes.
#' @inheritParams grid_contiguity_matrix
#' @inheritParams gen_context
#' @inheritParams gen_repartition
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
#' @keywords internal
#' @export
grid_simulation <- function(x_int, y_int, contiguityType = "Queen",
                            nbIndividuals = x_int * y_int,
                            nbMinEmptyZones = 0L,
                            nbMetropolises = 0L,
                            propInMetropolises = 0.4,
                            quantitatives_mat = c(0.0, 1.0),
                            qualitatives_list = NULL,
                            nbQuantitatives = 0L)
{
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)

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

  return(list("contiguity" = matCont, "repartition" = repartition_df,
              "context" = context_df))
}

#' A few grids to apply connectivity constraints.
#'
#' Each grid has a number of rows and a number of columns. An average of
#' 100 people per case are dispatched around the grid. No case are fixed as
#' empty or a a metropolis. The context is made of 2 quantitative variables and
#' one qualitative.
#' @docType data
#' @keywords data
#' @name grids
#' @rdname grids
#' @format for each grid a list with the following components:
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
NULL

# nolint start



# #' @rdname grids
# #' @keywords internal
# "grid_queen_vide0_metropole0_x2_y3_indivMoy100_quant3_qual0"

# set.seed(123L)
# grid_queen_vide0_metropole0_x2_y3_indivMoy100_quant3_qual0 = grid_simulation_1(2, 3)
# # use_data(grid_queen_vide0_metropole0_x2_y3_indivMoy100_quant3_qual0, overwrite = TRUE)
# save(grid_queen_vide0_metropole0_x2_y3_indivMoy100_quant3_qual0,
#      file = "simulations/grid_queen_vide0_metropole0_x2_y3_indivMoy100_quant3_qual0.Rdata")


# #' @rdname grids
# #' @keywords internal
# "grid_queen_vide0_metropole0_x3_y3_indivMoy100_quant3_qual0"

# set.seed(123L)
# # grid_queen_vide0_metropole0_x3_y3_indivMoy100_quant3_qual0 = grid_simulation_1(3, 3)
# # use_data(grid_queen_vide0_metropole0_x3_y3_indivMoy100_quant3_qual0, overwrite = TRUE)
# save(grid_queen_vide0_metropole0_x3_y3_indivMoy100_quant3_qual0,
#      file = "simulations/grid_queen_vide0_metropole0_x3_y3_indivMoy100_quant3_qual0.Rdata")


# #' @rdname grids
# #' @keywords internal
# "grid_queen_vide0_metropole0_x4_y4_indivMoy100_quant3_qual0"

# set.seed(123L)
# grid_queen_vide0_metropole0_x4_y4_indivMoy100_quant3_qual0 = grid_simulation_1(4, 4)
# # use_data(grid_queen_vide0_metropole0_x4_y4_indivMoy100_quant3_qual0, overwrite = TRUE)
# save(grid_queen_vide0_metropole0_x4_y4_indivMoy100_quant3_qual0,
#     file = "simulations/grid_queen_vide0_metropole0_x4_y4_indivMoy100_quant3_qual0.Rdata")


#' @rdname grids
#' @keywords internal
"grid_queen_vide0_metropole0_x7_y7_indivMoy100_quant3_qual0"

# set.seed(123L)
# grid_queen_vide0_metropole0_x7_y7_indivMoy100_quant3_qual0 = grid_simulation_1(7, 7)
# use_data(grid_queen_vide0_metropole0_x7_y7_indivMoy100_quant3_qual0, overwrite = TRUE)
# # save(grid_queen_vide0_metropole0_x7_y7_indivMoy100_quant3_qual0,
# #      file = "simulations/grid_queen_vide0_metropole0_x7_y7_indivMoy100_quant3_qual0.Rdata")

#' @rdname grids
#' @keywords internal
"grid_queen_vide0_metropole0_x20_y20_indivMoy100_quant3_qual0"

# set.seed(123L)
# grid_queen_vide0_metropole0_x20_y20_indivMoy100_quant3_qual0 = grid_simulation_1(20, 20)
# use_data(grid_queen_vide0_metropole0_x20_y20_indivMoy100_quant3_qual0, overwrite = TRUE)
# # save(grid_queen_vide0_metropole0_x20_y20_indivMoy100_quant3_qual0,
# #      file = "simulations/grid_queen_vide0_metropole0_x20_y20_indivMoy100_quant3_qual0.Rdata")

#' @rdname grids
#' @keywords internal
"grid_queen_vide0_metropole0_x30_y50_indivMoy100_quant3_qual0"

# set.seed(123L)
# grid_queen_vide0_metropole0_x30_y50_indivMoy100_quant3_qual0 = grid_simulation_1(30, 50)
# use_data(grid_queen_vide0_metropole0_x30_y50_indivMoy100_quant3_qual0, overwrite = TRUE)
# # save(grid_queen_vide0_metropole0_x30_y50_indivMoy100_quant3_qual0,
# #      file = "simulations/grid_queen_vide0_metropole0_x30_y50_indivMoy100_quant3_qual0.Rdata")

#' @rdname grids
#' @keywords internal
"grid_queen_vide0_metropole0_x50_y70_indivMoy100_quant3_qual0"

# set.seed(123L)
# grid_queen_vide0_metropole0_x50_y70_indivMoy100_quant3_qual0 = grid_simulation_1(50, 70)
# use_data(grid_queen_vide0_metropole0_x50_y70_indivMoy100_quant3_qual0, overwrite = TRUE)
# # save(grid_queen_vide0_metropole0_x50_y70_indivMoy100_quant3_qual0,
# #      file = "simulations/grid_queen_vide0_metropole0_x50_y70_indivMoy100_quant3_qual0.Rdata")


# set.seed(123L)
# grid_queen_vide0_metropole0_x175_y200_indivMoy100_quant3_qual0 = grid_simulation_1(175, 200)
# use_data(grid_queen_vide0_metropole0_x175_y200_indivMoy100_quant3_qual0, overwrite = TRUE)
# save(grid_queen_vide0_metropole0_x175_y200_indivMoy100_quant3_qual0,
#      file = "simulations/grid_queen_vide0_metropole0_x175_y200_indivMoy100_quant3_qual0.Rdata")

# nolint end

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

  return(list("contiguity" = mat_cont, "repartition" = repartition_df,
              "context" = context_df))

}


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
