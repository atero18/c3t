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


  n <- 5L
  # Test case 1: Valid partition with all constraints
  partition_valid <- c(1L, 1L, 2L, 3L, 2L)
  expect_true(is_feasible_solution(partition_valid))

  partition <- c(1L, 2L, 2L, 3L, 1L)
  # Test case 2: Valid partition with contiguity constraint
  contiguity_valid <- closest_neighbor_contiguity(n, loop = TRUE)
  expect_true(is_feasible_solution(partition, contiguity_valid))

  # Test case 3: Invalid partition with contiguity constraint
  contiguity_invalid <- closest_neighbor_contiguity(n, loop = FALSE)
  expect_false(is_feasible_solution(partition, contiguity_invalid))

  sizes <- c(2.0, 1.0, 3.0, 1.0, 2.0)

  # Test case 4: Valid partition with size constraints
  contiguity_valid <- closest_neighbor_contiguity(n, loop = TRUE)
  expect_true(is_feasible_solution(partition, sizes = sizes,
                                   m = 1.0, M = 4.0))

  # Test case 5: Invalid partition with size constraints
  expect_false(is_feasible_solution(partition, sizes = sizes,
                                    m = 4.0, M = 10.0))


  # Test case 6: Valid partition with minimum size constraint
  sizes_valid_min_size <- c(2.0, 1.0, 3.0, 2.0, 2.0)
  expect_true(is_feasible_solution(partition, sizes = sizes_valid_min_size,
                                   m = 2.0))


  # Test case 7: Invalid partition with minimum size constraint
  expect_false(is_feasible_solution(partition, sizes = sizes, m = 4.0))


  # Test case 8: Valid partition with maximum size constraint
  expect_true(is_feasible_solution(partition, sizes = sizes, M = 4.0))

  # Test case 9: Invalid partition with maximum size constraint
  sizes_invalid_max_size <- c(2.0, 1.0, 3.0, 1.0, 2.0)
  expect_false(is_feasible_solution(partition, sizes = sizes, M = 2.0))
})
