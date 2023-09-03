rm(list = ls())

# nolint start: T_and_F_symbol_linter
true_mat_rook <- matrix(c(F, T, T, F,
                          T, F, F, T,
                          T, F, F, T,
                          F, T, T, F),
                        ncol = 4L, nrow = 4L, byrow = TRUE)

# nolint end

# Creating a Contiguity Matrix
test_that("Creating a Contiguity Matrix",
{
  expect_no_error(contiguityMat(true_mat_rook))
})

matrix_contiguity <- contiguityMat(true_mat_rook)

# Accessing Data
test_that("Accessing Data",
{
  expect_true(matrix_contiguity[1L, 2L])
  expect_false(matrix_contiguity[2L, 3L])
  expect_false(matrix_contiguity[4L, 1L])
})

# Preservation of Symmetry
test_that("Preservation of Symmetry",
{
  expect_identical(matrix_contiguity[1L, 4L], matrix_contiguity[4L, 1L])
  expect_identical(matrix_contiguity[1L, 3L], matrix_contiguity[3L, 1L])
})

# The returned matrix in its entirety is correct
test_that("Complete matrix is correct",
{
  expect_equal(matrix_contiguity[], true_mat_rook, ignore_attr = TRUE)
})

# Adjacency Test
test_that("Adjacency",
{
  expect_true(are_adjacent(matrix_contiguity, 1L, 1L))
  expect_true(are_adjacent(matrix_contiguity, 1L, 2L))
  expect_false(are_adjacent(matrix_contiguity, 1L, 4L))
})
