rm(list = ls())

test_that("French (fr)",
{
  skip_if_not_installed("withr")

  withr::local_language("fr")

  distances <- matrix(c(0.0, 1.0,
                        1.0, 0.0), nrow = 2L)

  expect_error(constructor_DistMat(distances = distances, storageMode = "dsf"),
               regexp = "Mode de stockage de donn\u00e9es non reconnu")
})
