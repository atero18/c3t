rm(list = ls())

test_that("Unitary partition",
{
  pb <- gen_pb(7L, 7L)
  n <- pb$n()
  partition <- gen_initial_partition(pb, mode = "unitary")
  expect_integerish(partition, len = n)
  expect_identical(partition, seq_len(n))
})

test_that("Generate the righ number of fusions",
{
  pb <- gen_pb(7L, 7L)
  n <- pb$n()

  set.seed(123L)
  for (nbFusions in seq_len(n - 2L))
  {
    partition <- gen_initial_partition(pb, mode = "random", nbFusions)
    expect_integerish(partition, len = n)
    expect_identical(nbClusters(partition), n - nbFusions)
    expect_true(pb$isRegionalisation(partition))
  }
})

test_that("Generate the right number of partitions",
{
  pb <- gen_pb(7L, 7L)
  n <- pb$n()

  modes <- c("unitary", "random", "random", "unitary")
  nbFusions <- 10L
  set.seed(123L)

  partitions <- gen_initial_partitions(pb, modes, nbFusions)

  expect_list(partitions, len = 3L)

  expect_named(partitions, c("unitary", "random_1", "random_2"))

  nbClusters <- c(n, n - nbFusions, n - nbFusions)
  for (i in seq_along(partitions))
  {
    expect_integer(partitions[[i]], len = n,
                   any.missing = FALSE, all.missing = FALSE)
    expect_identical(nbClusters(partitions[[i]]), nbClusters[i])
  }
})
