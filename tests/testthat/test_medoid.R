rm(list = ls())

# Medoid
test_that("With one or two elements the medoid is the first",
{
  # One element
  expect_identical(medoid(as.matrix(12.0)), 1L)

  # Two elements
  distances <- matrix(c(0.0, 5.0,
                        5.0, 0.0), nrow = 2L)
  expect_identical(medoid(distances), 1L)
})

test_that("Correct value is returned",
{

  distances <- matrix(c( 0.0, 5.0, 10.0, # nolint; spaces_inside_linter
                         5.0, 0.0,  3.0,
                        10.0, 3.0,  0.0), nrow = 3L)
  colnames(distances) <- rownames(distances) <- c("2", "x", "1")

  # Without name
  expect_identical(medoid(distances, byName = FALSE), 2L)

  # With name
  expect_identical(medoid(distances, byName = TRUE), "x")
})
