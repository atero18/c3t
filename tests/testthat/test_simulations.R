rm(list = ls())

test_that("Repartition of individuals",
{
  n <- 100L
  N <- 1000L
  nbEmptyZones <- 10L
  nbMetropolises <- 7L
  pptInMet <- 0.5

  set.seed(123L)
  repartition <- gen_repartition(n, N, nbEmptyZones, nbMetropolises, pptInMet)

  expect_data_frame(repartition)
  expect_identical(ncol(repartition), 3L)
  expect_setequal(colnames(repartition),
                  c("nbIndividuals", "emptyFixedZone", "metropolis"))
  expect_identical(nrow(repartition), n)
  expect_integerish(repartition$nbIndividuals, lower = 0L)
  expect_identical(sum(repartition$nbIndividuals), N)
  expect_identical(sum(repartition$emptyFixedZone), nbEmptyZones)
  expect_identical(unique(repartition[repartition$emptyFixedZone,
                                      "nbIndividuals"]), 0L)
  expect_gte(sum(repartition$nbIndividuals == 0L), nbEmptyZones)
  expect_identical(sum(repartition$metropolis), nbMetropolises)
  expect_none(repartition$emptyFixedZone & repartition$metropolis)
})

test_that("Context generation",
{
  n <- 100L
  nbQuantitatives <- 2L
  qualitatives <- list(v1 = c("one", "two", "three"))

  set.seed(123L)
  context <- gen_context(n,
                         qualitatives_list = qualitatives,
                         nbQuantitatives_int = 2L)

  expect_data_frame(context, nrows = n,
                               ncols = nbQuantitatives + length(qualitatives))
  expect_factor(context$v1, levels = qualitatives$v1)
})
