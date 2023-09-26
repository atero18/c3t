rm(list = ls())

# Test that update_partition returns a numeric vector
test_that("update_partition returns a numeric vector",
{
  partition <- c(1L, 2L, 3L, 1L, 3L)
  new_partition <- update_partition(partition, 1L, 3L)
  expect_vector(partition, size = 5L)
  expect_numeric(new_partition)

})

# Test that update_partition updates the partition correctly
test_that("update_partition updates the partition correctly",
{
  partition <- c(1L, 2L, 3L, 1L, 3L)
  new_partition <- update_partition(partition, 1L, 3L)
  expect_equal(new_partition, c(1, 2, 1, 1, 1))
})

# Test that update_partition throws an error if cluster1 or cluster2 is not in partition
test_that("update_partition throws an error if cluster1 or cluster2 is not in partition",
{
  partition <- c(1L, 2L, 3L, 1L, 3L)
  expect_error(update_partition(partition, 1L, 4L))
})

# Test that update_sizes returns a numeric vector
test_that("update_sizes returns a numeric vector",
{
  sizes <- c(10.0, 20.0, 30.0, 40.0)
  cluster1 <- 1L
  cluster2 <- 3L
  newSizes = update_sizes(sizes, cluster1, cluster2)
  expect_vector(newSizes, size = 4L)
  expect_numeric(newSizes)
})

# Test that update_sizes updates the sizes correctly
test_that("update_sizes updates the sizes correctly",
{
  sizes <- c(10.0, 20.0, 30.0, 40.0)
  cluster1 <- 1L
  cluster2 <- 3L
  newSizes <- update_sizes(sizes, cluster1, cluster2)
  expect_equal(newSizes, c(40.0, 20.0, NaN, 40.0))
})
