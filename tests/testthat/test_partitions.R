rm(list = ls())

test_that("Are two partitions equivalent?",
{
  partition1 <- c(1L, 1L, 2L, 1L)

  expect_true(equivalent_partitions(partition1, partition1))

  partition2 <- c(3L, 3L, 2L, 3L)
  expect_true(equivalent_partitions(partition1, partition2))

  partition3 <- c("x", "x", "y", "x")
  expect_true(equivalent_partitions(partition1, partition3))

  partition4 <- c(1L, 1L, 2L, 2L)
  expect_false(equivalent_partitions(partition1, partition4))

  partition5 <- c(1L, 1L, 2L)
  expect_false(equivalent_partitions(partition1, partition5))

})
