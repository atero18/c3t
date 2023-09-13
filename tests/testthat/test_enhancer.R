rm(list = ls())

test_that("Test enhancment with classic criteria",
{
  set.seed(123L)
  x <- 4L
  y <- 5L
  nbIndividuals <- 100L
  nbCasesVides <- 3L
  nbMetropolises <- 2L
  nbVariablesQuant <- 2L
  grille <- gen_grid(x, y, nbIndividuals = nbIndividuals,
                     nbMinEmptyZones = nbCasesVides,
                     nbMetropolises = nbMetropolises,
                     nbQuantitatives = nbVariablesQuant)

  data <- grille$context
  individus <- grille$repartition$nbIndividuals
  contiguite <- grille$contiguity

  criteresAmelioration <- c("AHC", "CHI")

  m <- 5.0
  M <- 40.0

  d <- "euclidean"

  regionalisation <- c(1L, 2L, 2L, 1L, 1L, 1L, 3L, 1L, 3L, 3L,
                       1L, 3L, 1L, 1L, 1L, 3L, 1L, 1L, 1L, 3L)

  set.seed(123L)
  resEnhance <- enhance_feasible(regionalisation,
                                 contiguity = contiguite,
                                 d = d, data = data,
                                 sizes = individus,
                                 m = m, M = M,
                                 enhanceCriteria = criteresAmelioration,
                                 linkages = c("single", "complete"),
                                 parallele = FALSE,
                                 verbose = FALSE)

  expect_list(resEnhance,
              any.missing = FALSE, all.missing = FALSE,
              len = 2L)

  expect_named(resEnhance, c("initialValues", "results"))

  expect_tibble(resEnhance$results, nrows = 3L, ncols = 6L)
  expect_set_equal(colnames(resEnhance$results),
                   c("criterion", "linkage", "statut",
                     "iterations", "regionalisationOpti",
                     "dt_CHI"))

  expect_number(resEnhance$initialValues, lower = 0.0)
  expect_named(resEnhance$initialValues, "CHI")
  expect_integer(resEnhance$results$iterations, lower = 0L,
                 any.missing = FALSE, all.missing = FALSE)

  expect_double(resEnhance$results$dt_CHI,
                lower = -resEnhance$initialValues["CHI"],
                any.missing = FALSE, all.missing = FALSE)

})
