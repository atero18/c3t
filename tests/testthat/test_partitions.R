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


test_that("is_feasible_solution",
{

  # # Test case 1: Valid partition with all constraints
  # partition_valid <- c(1, 1, 2, 3, 2)
  # expect_true(is_feasible_solution(partition_valid))
  #
  # # Test case 2: Valid partition with contiguity constraint
  # partition_valid_contiguity <- c(1, 2, 2, 3, 1)
  # contiguity_valid <- matrix(c(F, 1, 1, 1, 0), nrow = 2)
  # expect_true(is_feasible_solution(partition_valid_contiguity, contiguity_valid))
  #
  # # Test case 3: Invalid partition with contiguity constraint
  # partition_invalid_contiguity <- c(1, 2, 2, 3, 1)
  # contiguity_invalid <- matrix(c(0, 1, 1, 0, 0), nrow = 2)
  # expect_false(is_feasible_solution(partition_invalid_contiguity, contiguity_invalid))Sure! Here are some more test cases:
  #
  # # Test case 4: Valid partition with size constraints
  # partition_valid_sizes <- c(1, 2, 2, 3, 1)
  # sizes_valid <- c(2, 1, 3, 1, 2)
  # expect_true(is_feasible_solution(partition_valid_sizes, sizes = sizes_valid, m = 1, M = 4))
  #
  # # Test case 5: Invalid partition with size constraints
  # partition_invalid_sizes <- c(1, 2, 2, 3, 1)
  # sizes_invalid <- c(2, 1, 3, 1, 2)
  # expect_false(is_feasible_solution(partition_invalid_sizes, sizes = sizes_invalid, m = 4, M = 10))
  #
  #
  # # Test case 6: Valid partition with minimum size constraint
  # partition_valid_min_size <- c(1, 2, 2, 3, 1)
  # sizes_valid_min_size <- c(2, 1, 3, 1, 2)
  # expect_true(is_feasible_solution(partition_valid_min_size, sizes = sizes_valid_min_size, m = 2))
  #
  # # Test case 7: Invalid partition with minimum size constraint
  # partition_invalid_min_size <- c(1, 2, 2, 3, 1)
  # sizes_invalid_min_size <- c(2, 1, 3, 1, 2)
  # expect_false(is_feasible_solution(partition_invalid_min_size, sizes = sizes_invalid_min_size, m = 4))
  #
  #
  # # Test case 8: Valid partition with maximum size constraint
  # partition_valid_max_size <- c(1, 2, 2, 3, 1)
  # sizes_valid_max_size <- c(2, 1, 3, 1, 2)
  # expect_true(is_feasible_solution(partition_valid_max_size, sizes = sizes_valid_max_size, M = 4))
  #
  # # Test case 9: Invalid partition with maximum size constraint
  # partition_invalid_max_size <- c(1, 2, 2, 3, 1)
  # sizes_invalid_max_size <- c(2, 1, 3, 1, 2)
  # expect_false(is_feasible_solution(partition_invalid_max_size, sizes = sizes_invalid_max_size, M = 2))
})
