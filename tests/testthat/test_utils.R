rm(list = ls())

test_that("which_na returns the correct indices of NA values",
{
  # Test case 1: mat with no NA values
  n <- 3L
  mat1 <- matrix(seq_len(n^2L), nrow = n)
  result1 <- which_na(mat1)
  expect_matrix(result1, nrows = 0L, ncols = 2L)

  # Test case 2: mat with NA values, removeSymmetry = TRUE,
  # and removeDiagonal = TRUE
  mat2 <- matrix(c(1.0,  NA, 3.0,
                   4.0, 5.0, 6.0,
                    NA, 8.0, 9.0),
                 nrow = 3L, ncol = 3L, byrow = TRUE)
  result2 <- unname(which_na(mat2,
                             removeSymmetry = TRUE, removeDiagonal = TRUE))
  expect_matrix(result2, nrows = 2L, ncols = 2L)
  expect_identical(result2[2L, ], c(1L, 2L))
  expect_identical(result2[1L, ], c(3L, 1L))

  # Test case 3: mat with NA values, removeSymmetry = FALSE,
  # and removeDiagonal = TRUE

  result3 <- unname(which_na(mat2,
                             removeSymmetry = FALSE, removeDiagonal = TRUE))
  expect_matrix(result3, nrows = 2L, ncols = 2L)
  expect_identical(result3[2L, ], c(1L, 2L))
  expect_identical(result3[1L, ], c(3L, 1L))

  # Test case 4: mat with NA values, removeSymmetry = TRUE,
  # removeDiagonal = FALSE
  result4 <- unname(which_na(mat2,
                             removeSymmetry = TRUE, removeDiagonal = FALSE))
  expect_matrix(result4, nrows = 2L, ncols = 2L)
  expect_identical(result4[2L, ], c(1L, 2L))
  expect_identical(result4[1L, ], c(3L, 1L))

  # Test case 5: mat with NA values, removeSymmetry = FALSE,
  # removeDiagonal = FALSE
  result5 <- unname(which_na(mat2,
                             removeSymmetry = FALSE, removeDiagonal = FALSE))
  expect_matrix(result5, nrows = 2L, ncols = 2L)
  expect_identical(result5[2L, ], c(1L, 2L))
  expect_identical(result5[1L, ], c(3L, 1L))

  # Test case 6: mat with only diagonal NA values, removeDiagonal = TRUE
  # nolint start: spaces_inside_linter
  mat6 <- matrix(c( NA, 2.0, 3.0,
                   4.0,  NA, 6.0,
                   7.0, 8.0,  NA), nrow = 3L)
  # nolint end
  result6 <- unname(which_na(mat6, removeDiagonal = TRUE))
  expect_matrix(result6, nrows = 0L, ncols = 2L)

  # Test case 7: mat with only diagonal NA values, removeDiagonal = FALSE
  result7 <- unname(which_na(mat6, removeDiagonal = FALSE))
  expect_matrix(result7, nrows = 3L, ncols = 2L)
  expect_identical(result7[1L, ], c(1L, 1L))
  expect_identical(result7[2L, ], c(2L, 2L))
  expect_identical(result7[3L, ], c(3L, 3L))

  # Test case 7: empty matrix
  mat8 <- matrix()
  result8 <- which_na(mat8)
  expect_matrix(result8, nrows = 0L, ncols = 2L)
})

test_that(paste0("which_na throws an error if removeSymmetry argument is not",
                 "a logical value"),
{
  # Test case 1: removeSymmetry argument is a character
  expect_error(which_na(matrix(), removeSymmetry = "TRUE"))

  # Test case 2: removeSymmetry argument is a numeric value
  expect_error(which_na(matrix(), removeSymmetry = 1L))
})

test_that(paste0("which_na throws an error if removeDiagonal argument is not",
                 " a logical value"),
{
  # Test case 1: removeDiagonal argument is a character
  expect_error(which_na(matrix(), removeDiagonal = "TRUE"))

  # Test case 2: removeDiagonal argument is a numeric value
  expect_error(which_na(matrix(), removeDiagonal = 1L))
})

test_that(paste0("which_na throws an error if removeSymmetry argument is not",
                 " a logical value"),
{
  # Test case 1: removeSymmetry argument is a character
  expect_error(which_na(matrix(), removeSymmetry = "TRUE"))

  # Test case 2: removeSymmetry argument is a numeric value
  expect_error(which_na(matrix(), removeSymmetry = 1L))
})

test_that(paste0("which_na throws an error if removeDiagonal argument is not",
                 " a logical value"),
{
  # Test case 1: removeDiagonal argument is a character
  expect_error(which_na(matrix(), removeDiagonal = "TRUE"))

  # Test case 2: removeDiagonal argument is a numeric value
  expect_error(which_na(matrix(), removeDiagonal = 1L))
})

test_that("which_na returns the correct results for large matrices",
{

  n <- 1000L
  # Test case 1: Large matrix with random NA values
  set.seed(123L)
  mat <- matrix(sample(c(1L:5L, NA_real_), size = n^2L, replace = TRUE),
                nrow = n)
  result <- which_na(mat, removeSymmetry = FALSE, removeDiagonal = FALSE)
  expect_matrix(result, nrows = sum(is.na(mat)))

  # Test case 2: Large matrix with no NA values
  mat2 <- matrix(seq_len(n^2L), nrow = n)
  result2 <- which_na(mat2)
  expect_matrix(result2, nrows = 0L, ncols = 2L)
})
