rm(list = ls())

# Problem Creation

X <- matrix(c(0.0, 2.0, 3.0,
              2.0, 0.0, 4.0,
              3.0, 4.0, 0.0), nrow = 3L, ncol = 3L)

# nolint start: T_and_F_symbol_linter
contiguity <- matrix(c(T, F, T,
                       F, T, T,
                       T, T, T), nrow = 3L)
# nolint end

sizes <- 1L:3L

test_that("With matrix (distances + contiguity), no constraints",
{
  expect_no_error(constructor_pbCon(X,
                                    contiguity = contiguity,
                                    sizes = sizes))
})

test_that("With contiguity graph",
{
  cont_graph <- contiguity_matrix_to_graph(contiguity)
  expect_no_error(constructor_pbCon(X,
                                    contiguity = cont_graph,
                                    sizes = sizes))
})

test_that("Without sizes",
{
  expect_no_error({pb <- constructor_pbCon(X, contiguity)})
  expect_identical(pb$N, as.double(nrow(X)))
})

# Constraint Elements Verification
test_that("Constraints must be numbers",
{
  expect_error(constructor_pbCon(X, contiguity, m = "a"))
  expect_error(constructor_pbCon(X, contiguity, M = 1L:3L))
})

test_that("Constraints must have correct values",
{
  expect_error(constructor_pbCon(X, contiguity,
                                 m = sum(sizes)) + 1.0)
  expect_error(constructor_pbCon(X, contiguity,
                                 M = 0.0))
  expect_error(constructor_pbCon(X, contiguity,
                                 sizes = sizes,
                                 M = max(sizes) - 1.0))
  expect_error(constructor_pbCon(X, contiguity,
                                 sizes = sizes,
                                 m = sum(sizes) + 2.0,
                                 M = max(sizes) + 1.0))
})

# Contiguity Verification
test_that("Contiguity information is correct",
{
  pb <- constructor_pbCon(X, contiguity)
  expect_true(are_adjacent(pb, 1L, 1L))
  expect_false(are_adjacent(pb, 1L, 2L))
  expect_true(are_adjacent(pb, 1L, 3L))
})

# Presence of Constraints Verification
test_that("Check for the presence of constraints",
{
  pb <- constructor_pbCon(X, contiguity)
  expect_false(pb$hasMinConstraint())
  expect_false(pb$hasMaxConstraint())
  expect_false(pb$hasSizeConstraint())

  pb <- constructor_pbCon(X, contiguity, m = 2.0)
  expect_true(pb$hasMinConstraint())
  expect_false(pb$hasMaxConstraint())
  expect_true(pb$hasSizeConstraint())

  pb <- constructor_pbCon(X, contiguity, M = 2.0)
  expect_false(pb$hasMinConstraint())
  expect_true(pb$hasMaxConstraint())
  expect_true(pb$hasSizeConstraint())

  pb <- constructor_pbCon(X, contiguity, m = 1.0, M = 2.0)
  expect_false(pb$hasMinConstraint())
  expect_true(pb$hasMaxConstraint())
  expect_true(pb$hasSizeConstraint())
})

# Testing problem splitting into connected subproblems
test_that("Splitting a problem into connected subproblems",
{
  # nolint start: T_and_F_symbol_linter
  contiguity <- matrix(c(F, T, F, F,
                         T, F, F, F,
                         F, F, F, T,
                         F, F, T, F), nrow = 4L, byrow = TRUE)
  # nolint end

  # nolint start: implicit_integer_linter
  distanceMatrix <- matrix(c(0, 1, 2, 3,
                             1, 0, 3, 2,
                             2, 3, 0, 4,
                             3, 2, 4, 0), nrow = 4L)
  # nolint end
  #
  sizes <- 1L:4L
  m <- 3L
  M <- 6.0

  pb <- constructor_pbCon(distanceMatrix, contiguity, sizes, m = m, M = M)
  connected_subproblems <- pb$split_connected_spb()
  expect_type(connected_subproblems, "list")
  expect_length(connected_subproblems, 2L)
  # nolint start: undesirable_function_linter
  . <- sapply(connected_subproblems, expect_s4_class, class = "pbCon")
  . <- sapply(connected_subproblems,
              function(sPb) expect_true(sPb$m %in% c(0.0, m)))
  . <- sapply(connected_subproblems,
              function(sPb) expect_true(sPb$M %in% c(M, Inf)))
  . <- sapply(connected_subproblems,
              function(sPb) expect_identical(nrow(sPb), 2L))
  . <- sapply(connected_subproblems,
              function(sPb) expect_length(V(sPb$contiguity), 2L))
  . <- sapply(connected_subproblems,
              function(sPb) expect_length(E(sPb$contiguity), 1L))
  # nolint end
  expect_identical(connected_subproblems[[1L]][], distanceMatrix[1L:2L, 1L:2L])
  expect_identical(connected_subproblems[[2L]][], distanceMatrix[3L:4L, 3L:4L])
  expect_identical(connected_subproblems[[1L]]$sizes, sizes[1L:2L])
  expect_identical(connected_subproblems[[2L]]$sizes, sizes[3L:4L])
  expect_named(connected_subproblems[[1L]], 1L:2L)
  expect_named(connected_subproblems[[2L]], 3L:4L)
})
