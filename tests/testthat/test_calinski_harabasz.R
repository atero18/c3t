rm(list = ls())

# Initial computation
test_that("Correct format of returned data - valueOnly = TRUE",

{
  n <- 3L
  data <- data.frame(a = seq_len(n), b = seq_len(n))
  partition <- c(1L, 2L, 2L)

  res <- calinski_harabasz(data, partition, valueOnly = TRUE)
  expect_vector(res)
  expect_type(res, "double")
  expect_length(res, 1L)
  expect_gte(res, 0.0)
})

test_that("Correct format of returned data - valueOnly = FALSE",
{
  n <- 3L
  data <- data.frame(a = seq_len(n), b = seq_len(n))
  partition <- c(1L, 2L, 2L)
  k <- length(unique(partition))


  res <- calinski_harabasz(data, partition, valueOnly = FALSE)
  expect_list(res)

  expect_vector(res$k, ptype = integer(), size = 1L)
  expect_identical(res$k, k)

  nbElementsClusters <- as.vector(table(partition))
  names(nbElementsClusters) <- 1L:2L
  expect_identical(res$nbElementsClusters, nbElementsClusters)

  expect_vector(res$betweenClusters, ptype = double(), size = k)
  expect_all(res$betweenClusters >= 0.0)
  expect_named(res$betweenClusters, as.character(1L:2L))

  expect_vector(res$B, ptype = double(), size = 1L)
  expect_gte(res$B, 0.0)

  expect_vector(res$withinCluster, ptype = double(), size = k)
  expect_all(res$withinCluster >= 0.0)
  expect_named(res$withinCluster, as.character(1L:2L))

  expect_vector(res$W, ptype = double(), size = 1L)
  expect_gte(res$W, 0.0)

  expect_vector(res$CHI, ptype = double(), size = 1L)
  expect_gte(res$CHI, 0.0)
})

test_that("Index value is independent of return type",
{
  n <- 3L
  data <- data.frame(a = seq_len(n), b = seq_len(n))
  partition <- c(1L, 2L, 2L)

  ICHvalueOnly <- calinski_harabasz(data, partition, valueOnly = TRUE)
  ICHComplete  <- calinski_harabasz(data, partition, valueOnly = FALSE)$CHI

  expect_identical(ICHComplete, ICHvalueOnly)
})

test_that("Equality with fpc package",
{
  skip_if_not_installed("fpc")

  set.seed(123L)
  data <- gen_pb(20L, 20L)$data
  n <- nrow(data)

  nbEssais <- 30L
  resfpc <- double(nbEssais)
  resc3t <- double(nbEssais)

  while (nbEssais > 0L)
  {

    nbElements <- sample(2L:n, size = 1L)
    partition <- sample(seq_len(nbElements), replace = TRUE, size = nbElements)
    partition <- standardize_partition(partition)
    resfpc[nbEssais] <- fpc::calinhara(data[seq_len(nbElements), ], partition)
    resc3t[nbEssais] <- calinski_harabasz(data[seq_len(nbElements), ],
                                          partition,
                                          valueOnly = TRUE)
    nbEssais <- nbEssais - 1L
  }

  expect_equal(resc3t, resfpc, tolerance = testthat_tolerance())
})

test_that("Result for specific values",
{
  n <- 3L
  data <- data.frame(a = seq_len(n), b = seq_len(n))

  # When n = 1
  partition <- 1L
  expect_identical(calinski_harabasz(data[1L, ], partition, valueOnly = TRUE),
                   Inf)

  # When k = 1
  partition <- rep(1L, n)

  expect_identical(calinski_harabasz(data, partition, valueOnly = TRUE), 0.0)

  # When k = n
  partition <- seq_len(n)

  expect_identical(calinski_harabasz(data, partition, valueOnly = TRUE), NaN)
})

test_that("Partition normalization is not needed",
{
  n <- 10L
  data <- data.frame(a = seq_len(n), b = seq_len(n))

  set.seed(123L)
  partition <- sample(seq_len(n), size = n, replace = TRUE)
  partition[1L] <- n + 1L

  idsClusters <- sort(unique(partition))

  resWithoutNormalization <- calinski_harabasz(data, partition,
                                               valueOnly = FALSE)

  ICHWithoutNormalization <- resWithoutNormalization$CHI

  expect_named(resWithoutNormalization$withinCluster,
               as.character(idsClusters))
  expect_named(resWithoutNormalization$betweenClusters,
               as.character(idsClusters))
  expect_named(resWithoutNormalization$nbElementsClusters,
               as.character(idsClusters))

  partition <- standardize_partition(partition)
  ICHWithNormalization <- calinski_harabasz(data, partition,
                                            valueOnly = TRUE)

  expect_identical(ICHWithoutNormalization, ICHWithNormalization)
})

# Update
test_that("Update ICH - Case when donor still has elements",
{
  n <- 5L
  data <- data.frame(a = seq_len(n), b = seq_len(n))

  partition <- c(1L, 1L, 3L, 3L, 2L)

  initialICH <- calinski_harabasz(data, partition, valueOnly = FALSE)

  donor <- 1L
  givenElements <- 2L
  receiver <- 2L


  updateICH <- update_calinski_harabasz(initialICH, donor, receiver,
                                        givenElements,
                                        partition, data)

  # Update with pbCon
  pb <- constructor_pbCon(d = "euclidean", data = data)
  updateICH2 <- pb$update_quality_partition(partition, "CHI", initialICH,
                                            donor, receiver, givenElements)

  partition[givenElements] <- receiver
  withoutUpdateICH <- calinski_harabasz(data, partition, valueOnly = FALSE)


  expect_identical(updateICH, withoutUpdateICH)
  expect_identical(updateICH, updateICH2)
})

test_that("Update ICH - Case when donor no longer has elements",
{
  n <- 5L
  data <- data.frame(a = seq_len(n), b = seq_len(n))

  partition <- c(2L, 1L, 3L, 3L, 2L)
  k <- length(unique(partition))

  initialICH <- calinski_harabasz(data, partition, valueOnly = FALSE)

  donor <- 1L
  givenElements <- 2L
  receiver <- 2L


  updateICH <- update_calinski_harabasz(initialICH, donor, receiver,
                                        givenElements,
                                        partition, data)


  partition[givenElements] <- receiver
  withoutUpdateICH <- calinski_harabasz(data, partition, valueOnly = FALSE)

  expect_identical(updateICH, withoutUpdateICH)
})
