rm(list = ls())

# Medoid
test_that("With one or two elements the medoid is the first",
{
  # One element
  expect_identical(medoid(as.matrix(12.0)), 1L)

  # Two elements
  distances <- matrix(c(0.0, 5.0,
                        5.0, 0.0), nrow = 2L)
  expect_identical(medoid(distances), 1L)
})

test_that("Correct value is returned",
{

  distances <- matrix(c( 0.0, 5.0, 10.0, # nolint; spaces_inside_linter
                         5.0, 0.0,  3.0,
                        10.0, 3.0,  0.0), nrow = 3L)
  colnames(distances) <- rownames(distances) <- c("2", "x", "1")

  # Without name
  expect_identical(medoid(distances, byName = FALSE), 2L)

  # With name
  expect_identical(medoid(distances, byName = TRUE), "x")
})

# Test medoids_partition function
test_that("medoids_partition function returns correct result", {
  # Test when the partition is empty
  partition <- NULL
  distances <- matrix(0.0, nrow = 0L, ncol = 0L)
  expect_error(medoids_partition(distances, partition))

  # Test when the partition has only 1 cluster
  partition <- c(1L, 1L, 1L)
  distances <- matrix(c(0.0, 1.0, 2.0,
                        1.0, 0.0, 3.0,
                        2.0, 3.0, 0.0), nrow = 3L, ncol = 3L)
  realMedoid <- c("1" = 1L)
  expect_identical(medoids_partition(distances, partition), realMedoid)

  # Test when the partition has multiple clusters
  n <- 6L
  partition <- c(1L, 1L, 1L, 2L, 2L, 2L)
  distances <- gen_distances(n)
  medoid1 <- which.min(rowSums(distances[1L:3L, 1L:3L]))
  medoid2 <- which.min(rowSums(distances[4L:6L, 4L:6L])) + 3L
  realMedoids <- c("1" = medoid1, "2" = medoid2)
  expect_identical(medoids_partition(distances, partition), realMedoids)
})
