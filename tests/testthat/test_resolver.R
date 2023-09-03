rm(list = ls())

# Strict Regularization (without size constraint relaxation)

test_that("No exchanges to be made",
{
  m <- 2.0
  sizes <- rep(1L, 4L)
  M <- 2.0

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F, F,
                         T, F, T, F,
                         F, T, F, T,
                         F, F, T, F), byrow = TRUE, nrow = 4L)
  # nolint end

  distances <- matrix(0.0, nrow = 4L, ncol = 4L)
  regionalisation <- c(1L, 1L, 2L, 2L)

  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation)
  expect_identical(res$status, "already_feasible")
  expect_identical(res$itFusions, 0L)
  expect_identical(res$itTransfers, 0L)
  expect_identical(res$regionalisation, regionalisation)
  expect_length(res$resolvedRegions, 0L)
  expect_null(res$finalSmallRegions)
  expect_null(res$finalBigRegions)

})

test_that("Cannot make any exchanges",
{
  # The first cluster lacks individuals (below min value)
  # and the second has too many (above max value). But if
  # the second gives its transferrable element to the first,
  # then the first cluster will have too many elements.
  m <- 2.0
  M <- 3.0
  sizes <- c(1.0, 0.0, 3.0, 3.0)

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F, F,
                         T, F, T, F,
                         F, T, F, T,
                         F, F, T, F), nrow = 4L, byrow = TRUE)
  # nolint end

  distances <- matrix(0.0, nrow = 4L, ncol = 4L)
  regionalisation <- c(1L, 1L, 2L, 2L)

  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation)

  expect_identical(res$status, "unresolvable")
  expect_identical(res$itFusions, 0L)
  expect_identical(res$itTransfers, 0L)
  expect_identical(res$regionalisation, regionalisation)
  expect_length(res$resolvedRegions, 0L)
  expect_identical(res$finalSmallRegions, 1L)
  expect_identical(res$finalBigRegions, 2L)

})

test_that("Impossible regularization",
{
  m <- 3.0
  M <- 3.0
  sizes <- c(0.0, 2.0, 2.0, 0.0)

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F, F,
                         T, F, T, F,
                         F, T, F, T,
                         F, F, T, F), nrow = 4L, byrow = TRUE)
  # nolint end

  distances <- matrix(0.0, nrow = 4L, ncol = 4L)
  regionalisation <- c(1L, 1L, 2L, 2L)

  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation)

  expect_identical(res$status, "unresolvable")
  expect_identical(res$itFusions, 0L)
  expect_identical(res$itTransfers, 0L)
  expect_identical(res$regionalisation, regionalisation)
  expect_null(res$resolvedRegions)
  expect_identical(res$finalSmallRegions, c(1L, 2L))
  expect_null(res$finalBigRegions)

})

test_that("Makes a single exchange",
{

  m <- 0.0
  M <- 2.0
  sizes <- rep(1.0, 4L)

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F, F,
                         T, F, T, F,
                         F, T, F, T,
                         F, F, T, F), nrow = 4L, byrow = TRUE)
  # nolint end

  distances <- matrix(0.0, nrow = 4L, ncol = 4L)

  regionalisation <- c(1L, 1L, 1L, 2L)
  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation)

  expect_identical(res$status, "fully_resolved")
  expect_identical(res$itTransfers, 1L)
  expect_identical(res$itFusions, 0L)
  expect_identical(res$regionalisation, c(1L, 1L, 2L, 2L))
  expect_length(res$resolvedRegions, 1L)
  expect_null(res$finalSmallRegions)
  expect_null(res$finalBigRegions)

})

test_that("Partial regularization",
{
  # The first cluster is too small while the second is of
  # appropriate size. The 3rd element will be given to the
  # first cluster, but this won't be enough: the only element
  # that can be given next is the 4th, which if removed from
  # the second cluster will result in a size too small for it.
  m <- 4.0
  M <- 10.0
  sizes <- c(1.0, 1.0, 1.0, 3.0, 1.0)

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F, F, F,
                         T, F, T, F, F,
                         F, T, F, T, F,
                         F, F, T, F, T,
                         F, F, F, T, F), byrow = TRUE, nrow = 5L)
  # nolint end

  distances <- matrix(0.0, nrow = 5L, ncol = 5L)
  regionalisation <- c(1L, 1L, 2L, 2L, 2L)
  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation)
  expect_identical(res$status, "partially_resolved")
  expect_identical(res$itTransfers, 1L)
  expect_identical(res$regionalisation, c(1L, 1L, 1L, 2L, 2L))
  expect_length(res$resolvedRegions, 0L)
  expect_identical(res$finalSmallRegions, 1L)
  expect_null(res$finalBigRegions)

})

test_that("Singleton transfer allowed",
{
  m <- 2.0
  M <- Inf
  sizes <- c(1.0, 2.0, 1.0)

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F,
                         T, F, T,
                         F, T, F), nrow = 3L, byrow = TRUE)
  # nolint end

  distances <- matrix(0.0, nrow = 3L, ncol = 3L)
  regionalisation <- c(1L, 2L, 2L)

  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation)

  expect_identical(res$status, "fully_resolved")
  expect_identical(res$itFusions, 0L)
  expect_identical(res$itTransfers, 1L)
  expect_identical(res$regionalisation, c(1L, 1L, 1L))
  expect_identical(res$resolvedRegions, 1L)
  expect_null(res$finalSmallRegions)
  expect_null(res$finalBigRegions)

})

test_that("Singleton transfer not allowed",
{
  m <- 2.0
  M <- Inf
  sizes <- c(1.0, 2.0, 1.0)

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F,
                         T, F, T,
                         F, T, F), nrow = 3L, byrow = TRUE)
  # nolint end

  distances <- matrix(0.0, nrow = 3L, ncol = 3L)
  regionalisation <- c(1L, 2L, 2L)

  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation,
                            allowSingletonTransfer = FALSE)
  expect_identical(res$status, "unresolvable")
  expect_identical(res$itTransfers, 0L)
  expect_identical(res$itFusions, 0L)
  expect_identical(res$regionalisation, c(1L, 2L, 2L))
  expect_null(res$resolvedRegions)
  expect_identical(res$finalSmallRegions, 1L)
  expect_null(res$finalBigRegions)

})

# Regularization with Max Relaxation

test_that("Impossible without max relaxation, works with",
{
  m <- 3.0
  M <- 3.0
  sizes <- c(2.0, 2.0, 0.0)

  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F,
                         T, F, T,
                         F, T, F), nrow = 3L, byrow = TRUE)
  # nolint end

  distances <- matrix(0.0, nrow = 3L, ncol = 3L)
  regionalisation <- c(1L, 2L, 2L)

  # Without relaxation
  res <- resolve_unfeasible(distances, contiguity,
                            sizes, m = m, M = M,
                            regionalisation = regionalisation)
  expect_identical(res$status, "unresolvable")
  expect_identical(res$itTransfers, 0L)
  expect_identical(res$itFusions, 0L)
  expect_identical(res$regionalisation, regionalisation)
  expect_null(res$resolvedRegions)
  expect_identical(res$finalSmallRegions, c(1L, 2L))
  expect_null(res$finalBigRegions)


  # With relaxation
  res <- resolve_unfeasible(distances, contiguity, sizes,
                            m = m, M = 4.0 / 3.0 * M,
                            regionalisation = regionalisation)
  expect_identical(res$status, "fully_resolved")
  expect_identical(res$itTransfers, 1L)
  expect_identical(res$itFusions, 0L)
  expect_identical(res$regionalisation, c(1L, 1L, 1L))
  expect_identical(res$resolvedRegions, c(1L, 2L))
  expect_null(res$finalSmallRegions)
  expect_null(res$finalBigRegions)

})
