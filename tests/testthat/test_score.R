rm(list = ls())

test_that("score_constraints_min calculates score correctly",
{
  clustersSizes <- c(5L, 6L, 3L, 2L, 4L)

  # Test case 1: m = 0, all clusters respect the constraint
  expect_equal(score_constraints_min(clustersSizes, m = 0L), 0.0)

  # Test case 2: m = 3, some clusters do not respect the constraint
  expect_equal(score_constraints_min(clustersSizes, m = 3L), 1.0)

  # Test case 3: m = 6, all clusters do not respect the constraint
  expect_equal(score_constraints_min(clustersSizes, m = 6L), 10.0)
})

test_that("score_constraints_max calculates score correctly",
{
  clustersSizes <- c(5L, 6L, 3L, 2L, 4L)

  # Test case 1: M = Inf, all clusters respect the constraint
  expect_equal(score_constraints_max(clustersSizes, M = Inf), 0.0)

  # Test case 2: M = 6, all clusters respect the constraint
  expect_equal(score_constraints_max(clustersSizes, M = 6L), 0.0)

  # Test case 3: M = 3, some clusters do not respect the constraint
  expect_equal(score_constraints_max(clustersSizes, M = 3L), 6.0)

  # Test case 4: M = 1, all clusters do not respect the constraint
  expect_equal(score_constraints_max(clustersSizes, M = 1L),
               sum(clustersSizes - 1.0))


})

test_that("score_constraints calculates score correctly",
{
  clustersSizes <- c(5L, 6L, 3L, 2L, 4L)

  # Test case 1: m = 0, M = Inf, all clusters respect the constraint
  expect_equal(score_constraints(clustersSizes, m = 0L, M = Inf), 0.0)

  # Test case 2: m = 3, M = 5, some clusters do not respect the constraint
  expect_equal(score_constraints(clustersSizes, m = 3L, M = 5L), 2.0)

  # Test case 3: m = 7, M = 10, all clusters do not respect the constraint
  expect_equal(score_constraints(clustersSizes, m = 7L, M = 10L), 15.0)
})
