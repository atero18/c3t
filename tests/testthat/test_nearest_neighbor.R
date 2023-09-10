rm(list = ls())

# Nearest Neighbor Calculation

# -- The results returned in a simple case are correct
test_that("3 elements",
{
  n <- 3L
  M <- matrix(c(0.0, 1.0, 2.0,
                1.0, 0.0, 3.0,
                2.0, 3.0, 0.0), nrow = n)

  nn <- nearest_neighbor(M, inner = TRUE)
  expect_vector(nn, size = n)
  expect_named(nn, as.character(seq_len(n)))
  expect_equal(nn, c(2L, 1L, 1L), ignore_attr = TRUE)
})

test_that("3 elements with names",
{
  n <- 3L
  M <- matrix(c(0.0, 1.0, 2.0,
                1.0, 0.0, 3.0,
                2.0, 3.0, 0.0), nrow = n)

  rownames(M) <- c("a", "b", "c")
  colnames(M) <- c("x", "y", "z")

  nn <- nearest_neighbor(M, inner = TRUE)
  expect_vector(nn, size = n)
  expect_named(nn, rownames(M))
  expect_equal(nn, colnames(M)[c(2L, 1L, 1L)], ignore_attr = TRUE)
})

test_that("5 elements with subset",
{
  n <- 5L
  # nolint start: implicit_integer_linter
  M <- matrix(c(0,1,2,3,4,
                1,0,4,5,6,
                2,4,0,6,7,
                3,5,6,0,8,
                4,6,7,8,0), nrow = n)
  # nolint end

  nn <- nearest_neighbor(M, inner = TRUE, subsetPoints = 2L:4L)
  expect_vector(nn, size = 3L)
  expect_named(nn, as.character(2L:4L))
  expect_equal(nn, c(1L, 1L, 1L), ignore_attr = TRUE)

  nn <- nearest_neighbor(M, inner = TRUE, subsetNeighbors = 2L:4L)
  expect_vector(nn, size = n)
  expect_named(nn, as.character(seq_len(n)))
  expect_equal(nn, c(2L, 3L, 2L, 2L, 2L), ignore_attr = TRUE)
})

test_that("Vector",
{
  n <- 3L
  M <- matrix(c(0.0, 1.0, 2.0,
                1.0, 0.0, 3.0,
                2.0, 3.0, 0.0), nrow = n)

  nn <- nearest_neighbor(M[1L, ], inner = FALSE, subsetNeighbors = 2L:3L)
  expect_vector(nn, size = 1L)
  expect_named(nn, "1")
  expect_equal(nn, 2L, ignore_attr = TRUE)
})

test_that("Contiguity Constraint",
{
  n <- 3L
  M <- matrix(c(0.0, 1.0, 2.0,
                1.0, 0.0, 3.0,
                2.0, 3.0, 0.0), nrow = n)

  # nolint start: T_and_F_symbol_linter
  C = matrix(c(TRUE, FALSE, TRUE,
               FALSE, TRUE, TRUE,
               TRUE, TRUE, TRUE), nrow = n)

  # nolint end

  nn <- nearest_neighbor(M, inner = TRUE, contiguity = C)
  expect_vector(nn, size = n)
  expect_named(nn, as.character(seq_len(n)))
  expect_equal(nn, c(3L, 3L, 1L), ignore_attr = TRUE)
})

test_that("Inner Importance",
{
  n <- 3L
  M <- matrix(c(0.0, 1.0, 2.0,
                1.0, 0.0, 3.0,
                2.0, 3.0, 0.0), nrow = n)

  nn <- nearest_neighbor(M, inner = TRUE)
  expect_length(nn, n)
  expect_equal(nn, c(2L, 1L, 1L), ignore_attr = TRUE)

  nn <- nearest_neighbor(M, inner = FALSE)
  expect_equal(nn, seq_len(n), ignore_attr = TRUE)

})

# Nearest Neighbor Update

test_that("Update - no linkage",
{
  n <- 5L

  # nolint start: implicit_integer_linter
  M <- matrix(c(0,1,2,3,4,
                1,0,4,5,6,
                2,4,0,6,7,
                3,5,6,0,8,
                4,6,7,8,0), nrow = n)
  # nolint end

  nn <- nearest_neighbor(M, inner = TRUE)

  c1 <- 2L
  c2 <- nn[c1] # 1L

  M[c2, -c2] <-  M[-c2, c2] <- NaN

  updateNN <- update_nearest_neighbor(M, c1, c2, nn)
  expect_vector(updateNN, size = n)
  expect_named(updateNN, as.character(seq_len(n)))
  expect_true(is.nan(updateNN[c2]))
  expect_all(updateNN[-c2] != c2)
  expect_true(updateNN[c1] != c2)

  updateNN2 <- nearest_neighbor(M, inner = TRUE)
  expect_identical(updateNN, updateNN2)
})

test_that("Update - single linkage",
{
  n <- 5L

  # nolint start: implicit_integer_linter
  M <- matrix(c(0,1,2,3,4,
                1,0,4,5,6,
                2,4,0,6,7,
                3,5,6,0,8,
                4,6,7,8,0), nrow = n)
  # nolint end

  nn <- nearest_neighbor(M, inner = TRUE)

  minLink <- min_distance(M)
  c1 <- minLink[1L]
  c2 <- minLink[2L]

  newM <- update_single_linkage(M, c1, c2)

  updateNN <- update_nearest_neighbor(newM, c1, c2, nn, linkage = "single")

  expect_vector(updateNN, size = n)
  expect_named(updateNN, as.character(seq_len(n)))
  expect_true(is.nan(updateNN[c2]))
  expect_all(updateNN[-c2] != c2)
  expect_true(updateNN[c1] != c2)

  updateNN2 <- nearest_neighbor(newM, inner = TRUE)
  expect_identical(updateNN, updateNN2)
})

test_that("Update - complete linkage",
{
  n <- 5L

  # nolint start: implicit_integer_linter
  M <- matrix(c(0,1,2,3,4,
                1,0,4,5,6,
                2,4,0,6,7,
                3,5,6,0,8,
                4,6,7,8,0), nrow = n)
  # nolint end

  nn <- nearest_neighbor(M, inner = TRUE)

  minLink <- min_distance(M)
  c1 <- minLink[1L]
  c2 <- minLink[2L]

  newM <- update_complete_linkage(M, c1, c2)

  updateNN <- update_nearest_neighbor(newM, c1, c2, nn, linkage = "complete")

  expect_vector(updateNN, size = n)
  expect_named(updateNN, as.character(seq_len(n)))
  expect_true(is.nan(updateNN[c2]))
  expect_all(updateNN[-c2] != c2)
  expect_true(updateNN[c1] != c2)

  updateNN2 <- nearest_neighbor(newM, inner = TRUE)
  expect_identical(updateNN, updateNN2)
})
