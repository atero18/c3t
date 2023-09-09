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

test_that("merge_cc_paritions returns merged partition vector",
{
  # Test case 1: Only one partition, should return the same partition
  components <- c(1L, 1L, 1L, 1L)
  partitions <- list(1L:4L)
  result <- merge_cc_partitions(components, partitions)
  expect_identical(result, 1L:4L)

  # Test case 2: Multiple partitions, should merge them correctly
  components <- c(1L, 2L, 1L, 2L, 3L, 3L)
  partitions <- list(c(1L, 3L), c(2L, 2L), c(2L, 4L))
  expected <- c(1L, 2L, 3L, 2L, 4L, 5L)
  result <- merge_cc_partitions(components, partitions)
  expect_identical(result, expected)

  # Test case 3: Incorrect number of elements in partitions
   components <- c(1L, 1L, 1L, 1L)
   partitions <- list(1L:3L)  # Incorrect number of elements
   expect_error(merge_cc_partitions(components, partitions))

   # Test case 4: Empty components vector
   components <- integer(0L)
   partitions <- list()
   expect_error(merge_cc_partitions(components, partitions))

   # Test case 5: Empty partitions list
   components <- 1L:5L
   partitions <- list()
   expect_error(merge_cc_partitions(components, partitions))
})
