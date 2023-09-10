rm(list = ls())

COLNAMES <- c("x1", "x2", "cluster1", "cluster2")
expect_matrix_boundary <- function(x, nrows = 0L,
                                   partition = NULL, cluster_int = NULL)
{
  assertCount(nrows)

  any.missing <- nrows == 0L
  all.missing <- nrows == 0L
  ncols <- length(COLNAMES)
  expect_matrix(x,
                mode = "integer",
                nrows = nrows, ncols = ncols,
                any.missing = any.missing, all.missing = all.missing,
                null.ok = FALSE)

  expect_identical(colnames(x), COLNAMES)

  if (nrow(x) > 0L && !is.null(cluster_int))
  {
    expect_all(x[, 3L] == cluster_int | x[, 4L]) # nolint: object_usage_linter
    expect_xor(partition[x[, 1L]] == cluster_int, # nolint: object_usage_linter
               partition[x[, 2L]] == cluster_int)
  }

}

test_that("Test interior_boundary_cluster function",
{

  # Test case 1: empty cluster
  partition <- c(1L, 1L, 1L, 1L)
  contiguity <- complete_contiguity_matrix(4L)
  cluster_int <- 2L
  expect_warning({actual_result <-
    interior_boundary_cluster(contiguity, cluster_int, partition)})
  expect_matrix_boundary(actual_result, nrows = 0L, partition, cluster_int)

  # Test case 2: cluster with one element
  partition <- c(1L, 2L, 2L, 3L)
  contiguity <- complete_contiguity_matrix(4L)
  cluster_int <- 3L
  actual_result <- interior_boundary_cluster(contiguity, cluster_int, partition)
  expect_matrix_boundary(actual_result, nrows = 3L, partition, cluster_int)

  # Test case 3: cluster with multiple elements
  partition <- c(1L, 2L, 2L, 3L, 1L, 3L)
  contiguity <- matrix(FALSE, nrow = 6L, ncol = 6L)
  contiguity[2L, 3L] <- contiguity[3L, 2L] <- TRUE
  contiguity[3L, 4L] <- contiguity[4L, 3L] <- TRUE
  contiguity[5L, 6L] <- contiguity[6L, 5L] <- TRUE
  cluster_int <- 2L
  actual_result <- interior_boundary_cluster(contiguity, cluster_int, partition)
  expect_matrix_boundary(actual_result, nrows = 1L, partition, cluster_int)

})

test_that("Test interior_boundary function",
{

  # Test case 1: contiguity matrix
  partition <- c(1L, 1L, 1L, 1L)
  contiguity <- matrix(FALSE, nrow = 4L, ncol = 4L)
  contiguity[2L, 3L] <- contiguity[3L, 2L] <- TRUE
  actual_result <- interior_boundary(contiguity, partition)
  expect_matrix_boundary(actual_result, nrows = 0L)
  # Test case 2: contiguity graph
  g <- contiguity_matrix_to_graph(contiguity)
  partition <- c(1L, 1L, 2L, 2L, 3L)
  actual_result <- interior_boundary(g, partition)
  expect_matrix_boundary(actual_result, nrows = 1L)
})
