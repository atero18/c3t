rm(list = ls())

set.seed(1234L)

# Comparison with existing (implies a set where all elements are
# pairwise contiguous, meaning the contiguity constraint disappears)

n <- 5L
cont <- complete_contiguity_matrix(n)
indivs <- seq_len(n)
M <- sum(indivs)

donnees <- gen_context(n, quantitatives_mat = matrix(c( 10.0,  2.0, # nolint: spaces_inside_linter
                                                        -2.0,  4.0,
                                                      1000.0, 40.0),
                                                    ncol = 2L))
d <- euclidean_distance
dist_res <- dist(donnees, method = "euclidean")
pb <- constructor_pbCon(d = d, data = donnees,
                        sizes = indivs, contiguity = cont)

test_that("Linkage distance - single linkage",
{
  resPerso <- AHR_single(pb, linkage = "single")
  calculR <- hclust(dist_res, method = "single")

  for (k in 2L:4L)
    expect_identical(cutree(resPerso, k), cutree(calculR, k))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))

})

test_that("Linkage distance - complete linkage - diameter",
{
  resPerso <- AHR_single(pb, linkage = "complete")
  calculR <- hclust(dist_res, method = "complete")

  for (k in 2L:4L)
    expect_identical(cutree(resPerso, k), cutree(calculR, k))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))
})

test_that("Linkage distance - centroid linkage",
{
  resPerso <- AHR_single(pb, linkage = "centroid")
  calculR <- hclust(dist_res, method = "centroid")

  for (k in 2L:4L)
    expect_identical(cutree(resPerso, k), cutree(calculR, k))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))

})

test_that("Linkage distance - average linkage",
{
  resPerso <- AHR_single(pb, linkage = "average")
  calculR <- hclust(dist_res, method = "average")

  for (k in 2L:4L)
    expect_identical(cutree(resPerso, k), cutree(calculR, k))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))
})

test_that("Elements distance (`d`) - Minkowski p = 3",
{
  d <- function_distance_minkowski(3.0)
  pb <- constructor_pbCon(d = d, data = donnees,
                          sizes = indivs, contiguity = cont)
  dist_res <- dist(donnees, method = "minkowski", p = 3.0)
  resPerso <- AHR_single(pb, linkage = "single")
  calculR <- hclust(dist_res, method = "single")

  for (k in 2L:4L)
    expect_identical(cutree(resPerso, k), cutree(calculR, k))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))

})

test_that("Elements distance (`d`) - Manhattan",
{
  d <- distance_manhattan
  pb <- constructor_pbCon(d = d, data = donnees,
                          sizes = indivs, contiguity  = cont)
  dist_res <- dist(donnees, method = "manhattan")
  resPerso <- AHR_single(pb, linkage = "single")
  calculR <- hclust(dist_res, method = "single")

  for (k in 2L:4L)
    expect_identical(cutree(resPerso, k), cutree(calculR, k))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))
})

test_that("Elements distance (`d`) - supremum distance",
{
  d <- distance_sup
  pb <- constructor_pbCon(d = d, data = donnees,
                          sizes = indivs, contiguity = cont)

  dist_res <- dist(donnees, method = "maximum")
  resPerso <- AHR_single(pb, linkage = "single")
  calculR <- hclust(dist_res, method = "single")

  for (k in 2L:4L)
    expect_identical(cutree(resPerso, k), cutree(calculR, k))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))
})

# Verification of the contiguity constraint

test_that("Checking contiguity on a small set",
{
  n <- 3L
  # nolint start: spaces_inside_linter
  matDistances <- matrix(c( 0.0, 50.0, 1.0,
                            50.0,  0.0, 2.0,
                            1.0,  2.0, 0.0), nrow = 3L)
  # nolint end


  distances <- constructor_DistMat(matDistances)
  # nolint start: spaces_inside_linter
  data <- gen_context(n, quantitatives_mat = matrix(c(    0.0,  1.0,
                                                         -2.0,  4.0,
                                                       1000.0, 40.0),
                                                    ncol = 2L))
  # nolint end

  # nolint start: T_and_F_symbol
  mat_cont = matrix(c(T,T,F,
                      T,T,F,
                      F,F,T), nrow = 3L, ncol = 3L)
  # nolint end


  # 1 and 2 are contiguous, and 3 is not contiguous with anyone. Thus, we have
  # two connected components, {1,2} and {3}
  sizes <- seq_len(n)

  # Without the contiguity constraint, it is preferable to link 1 to 3
  # and leave 2 alone. However, with the constraint, there is no choice: the
  # only possible linkage is {1,2} and {3}.
  pb <- constructor_pbCon(distances = matDistances, data = data,
                          sizes = sizes, contiguity = mat_cont)

  resPerso <- AHR_single(pb)
  expect_identical(cutree(resPerso, 2L), c(1L, 1L, 2L))

  for (Partition in resPerso$partitions)
    expect_true(pb$isFeasibleSolution(Partition$partition))
})

test_that("Checking contiguity on a larger set",
{
  pb7x7 <- gen_pb(7L, 7L, m = 0.0, M = Inf, d = "euclidean")

  resPerso <- AHR_single(pb7x7)

  for (k in names(resPerso))
    expect_true(pb7x7$isRegionalisation(resPerso[[k]]$partition))
})


# Verification of the maximum constraint

test_that("With an excessively large max constraint, no cluster can be formed",
{
  n <- nrow(donnees)
  sizes <- rep(1L, n)
  M <- 1.0
  pb <- constructor_pbCon(d = d, data = donnees,
                          sizes = sizes, contiguity = cont, M = M)

  expect_length(AHR_single(pb), 1L)
})

test_that("Functioning of the grid",
{
  pb7x7 <- gen_pb(7L, 7L, d = "euclidean")

  d <- pb7x7$d
  context <- pb7x7$data
  contiguity <- pb7x7$contiguity
  m <- 3.0
  M <- Inf

  res <- AHR(d = d, data = context, contiguity = contiguity, m = m, M = M,
             parallele = FALSE, verbose = FALSE)

  expect_list(res, types = "tbl", len = 3L)
  expect_named(res, c("grid", "initialPartitions", "results"))

})

# Same result independently of the complete linkage calculation
test_that("Same result independently of the entire linkage calculation",
{
  pb7x7 <- gen_pb(7L, 7L, d = "euclidean", calculateAllDistances = FALSE)
  pb7x7BIS <- pb7x7$copy(shallow = FALSE)

  partition <- seq_len(49L)
  ##partition <- rep(seq_len(7L), each = 7L)
  res1 <- AHR_single(pb7x7, linkage = "single",
                     partitionInit = partition, calculDistComplet = FALSE)
  res2 <- AHR_single(pb7x7BIS, linkage = "single",
                     partitionInit = partition, calculDistComplet = TRUE)

  expect_true(all.equal(res1, res2))
  for (i in names(res1))
    expect_identical(res1[[i]]$partition,
                     res2[[i]]$partition)
})
