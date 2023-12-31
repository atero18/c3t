rm(list = ls())

X <- matrix(c(0.0, 2.0, 3.0,
              2.0, 0.0, 4.0,
              3.0, 4.0, 0.0), nrow = 3L, ncol = 3L)
rownames(X) <- colnames(X) <- 1L:3L

# Creation using an existing distance matrix and through a function
test_that("Creation using a distance matrix",
{
  expect_no_error(constructor_DistMat(distances = X))
})

M <- constructor_DistMat(X)

test_that("Creation using a matrix via as",
{
  expect_no_error(as(X, "DistMat"))
  expect_true(all.equal(M, as(X, "DistMat")))
})

d <- euclidean_distance
data <- data.frame(a = 1L:3L, b = 5L:7L)

test_that("Creation using a distance function",
{
  expect_no_error(constructor_DistMat(d = d, data = data))
})

test_that("Using a function requires data",
{
  expect_error(constructor_DistMat(d = d))
})

test_that("The provided data frame must not be empty",
{
  expect_error(constructor_DistMat(d = d,
                                  data = data.frame(a = numeric(0L),
                                                    b = numeric(0L))))
})

# Property checks

test_that("DistMat is a matrix",
{
  expect_true(is.matrix(M))
})

test_that("DistMat is symmetric",
{
  expect_symmetric(M)
})

test_that("Returned dimensions are correct",
{
  expect_length(dim(M), 2L)
  expect_identical(dim(M), c(3L, 3L))
  expect_identical(nrow(M), 3L)
  expect_identical(ncol(M), 3L)
})

test_that("DistMat is numeric",
{
  expect_true(is.numeric(M)) # nolint: expect_type_linter
})

# Value checks

test_that("Accessing data one by one",
{
  expect_identical(M[1L, 1L], 0.0)
  expect_identical(M[1L, 2L], X[1L, 2L])
  expect_identical(M[2L, 1L], X[1L, 2L])
})

test_that("Accessing a single row or column",
{
  expect_all(M[1L, ] == X[1L, ])
  expect_all(M[, 1L] == X[, 1L])
})

test_that("Returned sub-matrix is correct",
{
  expect_all(M[1L, 1L:2L] == X[1L, 1L:2L])
  expect_all(M[1L:2L, 1L] == X[1L:2L, 1L])
})

test_that("Returned matrix is indeed the distance matrix",
{
  expect_identical(M[], X)
})

# Equality operations
test_that("Equality (distance matrices)",
{
  expect_equal(M, M, ignore_attr = TRUE)
  U <- M$copy()
  U[1L, 2L] <- 5.0
  expect_false(isTRUE(all.equal(M, U)))
  U$d <- euclidean_distance
  U$data <- data
  expect_false(isTRUE(all.equal(M, U)))
})

test_that("Equality (distance matrix vs. matrix)",
{
  expect_true(all.equal(M, M[]))
  expect_true(all.equal(M[], M))
  U <- X
  U[1L, 2L] <- M[1L, 2L] + 1.0
  expect_false(isTRUE(all.equal(M, U)))
})


# Assignment checks

test_that("Assignments are verified",
{
  expect_no_error({M[1L, 1L] <- 0.0})
  expect_error({M[1L, 1L] <- 5.0})
  expect_no_error({M[1L, 2L] <- 2.0})
  expect_error({M[1L, 2L] <- -2.0})
})

test_that("Symmetry is preserved after assignment",
{
  newValue <- 7.0
  M[1L, 2L] <- newValue
  expect_identical(M[2L, 1L], newValue)
})

test_that("Correct calculation using a distance function",
{
  mat <- constructor_DistMat(d = d, data = data)
  expect_identical(mat[1L, 2L], euclidean_distance(data[1L, ], data[2L, ]))
  expect_identical(mat[2L, 1L], mat[1L, 2L])
})

# Check calculation using a distance function
test_that("The used function should return positive values",
{
  dNeg <- function(x, y) {-1.0}
  mat <- constructor_DistMat(d = dNeg, data = data)
  expect_error(mat[1L, 2L])
})

test_that("Comparison for function calculation with existing",
{
  n <- 10L
  dataVals <- data.frame(a = seq(1L, n),
                         b = seq(-n + 1L, 0L),
                         c = seq(2L * n + 1L, 3L * n))

  d <- euclidean_distance
  M <- constructor_DistMat(d = d, data = dataVals)
  expect_identical(M[], as.matrix(dist(dataVals)))
})

test_that("Linkage distances do not depend of storage mode",
{
  cppModes <- LINKAGES$cppMode
  cppModes <- cppModes[!is.na(cppModes)]
  set.seed(123L)
  for (n in seq(10L, 100L, 10L))
  {
    distances <- gen_distances(n)
    distances_vec <- as(distances, "SymVMat")
    distances_mat <- as(distances, "SymMMat")

    partition <- sample(seq_len(n), replace = TRUE, size = n)
    partition <- standardize_partition(partition)
    k <- nbClusters(partition)

    nbTries <- sample(seq_len(n), size = 1L)

    distancesToCalc <-
      matrix(sample(seq_len(k), replace = TRUE, size = 2L * nbTries),
             nrow = nbTries, ncol = 2L)

    equality <- distancesToCalc[, 1L] == distancesToCalc[, 2L]

    if (any(equality))
      distancesToCalc <- distancesToCalc[!equality, , drop = FALSE]

    for (mode in cppModes)
    {
      resVec <- calculDistancesInter_Cpp(distances_vec, partition,
                                         distancesToCalc,
                                         mode)

      resMat <- calculDistancesInter_Cpp(distances_mat, partition,
                                         distancesToCalc,
                                         mode)

      expect_equal(resVec, resMat, tolerance = 10L^-6L)

    }

  }

})
