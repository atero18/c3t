rm(list = ls())

test_that("Dunn index exists",
{
  expect_in("Dunn", available_criteria())
})

# Test calculating the Dunn index

test_that("2 clusters, 1 element per cluster",
{
  M <- matrix(c(0.0, 100.0,
                100.0, 0.0), nrow = 2L, byrow = TRUE)

  partition <- c(1L, 2L)

  dunn <- dunn_index(M, partition, linkage = "single")
  expect_type(dunn, "list")
  expect_equal(dunn$diameters, c(0.0, 0.0), ignore_attr = TRUE)
  expect_identical(dunn$compactness, 0.0)
  expect_true(is.matrix(dunn$linkageDistances))
  expect_length(dunn$linkageDistances, 4L)

  interDistMin <- matrix(c(Inf, 100.0,
                           100.0, Inf), nrow = 2L)
  expect_true(all(dunn$linkageDistances == interDistMin))

  expect_identical(dunn$separation, 100.0)
  expect_identical(dunn$D, Inf)

})

test_that("3 clusters, 2 elements per cluster",
{
  # nolint start: implicit_integer_linter
  M <- matrix(c(0, 1, 2, 3, 4, 5,
                1, 0, 1, 2, 3, 4,
                2, 1, 0, 1, 2, 3,
                3, 2, 1, 0, 1, 2,
                4, 3, 2, 1, 0, 1,
                5, 4, 3, 2, 1, 0), nrow = 6L)
  # nolint end


  partition <- c(1L, 1L, 2L, 2L, 3L, 3L)

  dunn <- dunn_index(M, partition, linkage = "single")

  # Diameter of each cluster: 1 by construction
  expect_equal(dunn$diameters, c(1L, 1L, 1L), ignore_attr = TRUE)
  expect_identical(dunn$compactness, 1.0)

  # Tests with min linkage distance
  dunnMin <- dunn
  interDistMin <- matrix(c(Inf, 1.0, 3.0,
                           1.0, Inf, 1.0,
                           3.0, 1.0, Inf), nrow = 3L)

  expect_true(all(dunnMin$linkageDistances == interDistMin))
  expect_identical(dunnMin$separation, 1.0)
  expect_identical(dunnMin$D, 1.0)

  # Tests with max linkage distance
  dunnMax <- dunn_index(M, partition, linkage = "complete")
  interDistMax <- matrix(c(Inf, 3.0, 5.0,
                           3.0, Inf, 3.0,
                           5.0, 3.0, Inf), nrow = 3L)

  expect_true(all(dunnMax$linkageDistances == interDistMax))
  expect_identical(dunnMax$separation, 3.0)
  expect_identical(dunnMax$D, 3.0)

})

test_that("3 clusters, 2 elements per cluster - indirect version",
{
  # nolint start: implicit_integer_linter
  M <- matrix(c(0, 1, 2, 3, 4, 5,
                1, 0, 1, 2, 3, 4,
                2, 1, 0, 1, 2, 3,
                3, 2, 1, 0, 1, 2,
                4, 3, 2, 1, 0, 1,
                5, 4, 3, 2, 1, 0), nrow = 6L)
  # nolint end


  partition <- c(1L, 1L, 2L, 2L, 3L, 3L)

  # Tests with min linkage distance
  dunnMin <- clustering_criterion(partition, "Dunn",
                                  M, linkage = "single")

  expect_identical(dunnMin, 1.0)

  # Tests with max linkage distance
  dunnMax <- clustering_criterion(partition, "Dunn",
                                  M, linkage = "complete")

  expect_identical(dunnMax, 3.0)

})

# Test updating the Dunn index
test_that("Updating Dunn index - min linkage distance",
{
  # nolint start: implicit_integer_linter
  M <- matrix(c(0, 1, 2, 3, 4, 5,
                1, 0, 1, 2, 3, 4,
                2, 1, 0, 1, 2, 3,
                3, 2, 1, 0, 1, 2,
                4, 3, 2, 1, 0, 1,
                5, 4, 3, 2, 1, 0), nrow = 6L)
  # nolint end


  partition <- c(1L, 1L, 2L, 2L, 3L, 3L)

  dunn <- dunn_index(M, partition, linkage = "single")

  # Transfer the 4th element from cluster 2 to cluster 3
  dunnUpdated <- update_dunn_index(dunn, M, partition, 2L, 3L, 4L)
  expect_type(dunnUpdated, "list")
  expect_equal(dunnUpdated$diameters, c(1.0, 0.0, 2.0), ignore_attr = TRUE)
  expect_identical(dunnUpdated$compactness, 2.0)

  interDistMinUpdated <- matrix(c(Inf, 1.0, 2.0,
                                  1.0, Inf, 1.0,
                                  2.0, 1.0, Inf), nrow = 3L)

  expect_true(all(dunnUpdated$linkageDistances == interDistMinUpdated))
  expect_identical(dunnUpdated$separation, 1.0)
  expect_identical(dunnUpdated$D, 1.0 / 2.0)
})

test_that("Updating Dunn index - max linkage distance",
{
  # nolint start: implicit_integer_linter
  M <- matrix(c(0, 1, 2, 3, 4, 5,
                1, 0, 1, 2, 3, 4,
                2, 1, 0, 1, 2, 3,
                3, 2, 1, 0, 1, 2,
                4, 3, 2, 1, 0, 1,
                5, 4, 3, 2, 1, 0), nrow = 6L)
  # nolint end

  partition <- c(1L, 1L, 2L, 2L, 3L, 3L)

  dunn <- dunn_index(M, partition, linkage = "complete")

  # Transfer the 4th element from cluster 2 to cluster 3
  dunnUpdated <- update_dunn_index(dunn, M, partition, 2L, 3L, 4L)
  expect_type(dunnUpdated, "list")
  expect_equal(dunnUpdated$diameters, c(1.0, 0.0, 2.0), ignore_attr = TRUE)
  expect_identical(dunnUpdated$compactness, 2.0)

  interDistMaxUpdated <- matrix(c(Inf, 2.0, 5.0,
                                  2.0, Inf, 3.0,
                                  5.0, 3.0, Inf), nrow = 3L)

  expect_true(all(dunnUpdated$linkageDistances == interDistMaxUpdated))
  expect_identical(dunnUpdated$separation, 2.0)
  expect_identical(dunnUpdated$D, 1.0)
})

# Test contiguous Dunn calculation
test_that("Calculating contiguous Dunn index - min linkage distance",
{
  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F, F, F, F,
                         T, F, T, F, F, F,
                         F, T, F, F, F, F,
                         F, F, F, F, T, F,
                         F, F, F, T, F, T,
                         F, F, F, F, T, F),
                       nrow = 6L, byrow = TRUE)
  # nolint end

  # nolint start: implicit_integer_linter
  distMatrix <- matrix(c(0, 1, 2, 3, 3, 2,
                         1, 0, 3, 2, 4, 7,
                         2, 3, 0, 4, 5, 1,
                         3, 2, 4, 0, 6, 3,
                         3, 4, 5, 6, 0, 5,
                         2, 7, 1, 3, 5, 0), nrow = 6L)
  # nolint end

  sizes <- 1L:6L
  m <- 3.0
  M <- 6.0
  pb <- constructor_pbCon(distMatrix, contiguity, sizes, m = m, M = M)

  partition <- c(1L, 1L, 3L, 2L, 2L, 3L)

  dunn <- pb$quality_partition(partition, "Dunn",
                               linkage = "single", connected = TRUE)

  dunnCC1 <- dunn_index(pb[1L:3L, 1L:3L],
                        partition[1L:3L], linkage = "single",
                        valueOnly = TRUE)
  dunnCC2 <- dunn_index(pb[4L:6L, 4L:6L],
                        partition[4L:6L], linkage = "single",
                        valueOnly = TRUE)

  expect_identical(dunn, mean(c(dunnCC1, dunnCC2)))
})
