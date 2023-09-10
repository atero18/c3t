rm(list = ls())

# Strict Regularization (without size constraint relaxation)

NOMINCONSTRAINT <- 0.0
NOMAXCONSTRAINT <- Inf

expect_resolver_return <-
  function(m, M, sizes,
           regionalisation,
           status, itFusions = 0L, itTransfers = 0L,
           newRegionalisation = regionalisation,
           nbResolvedRegions = 0L,
           finalSmallRegions = NULL, finalBigRegions = NULL,
           contiguity = closest_neighbor_contiguity(n),
           distances = empty_distance_matrix(n),
           allowSingletonTransfert = TRUE)
{
  n <- length(sizes)
  set.seed(123L)
  res <- resolve_unfeasible(distances, contiguity, sizes, m = m, M = M,
                            regionalisation = regionalisation,
                            allowSingletonTransfert = allowSingletonTransfert)


  initialProblematicRegions <- c(res$initialSmallRegions, res$initialBigRegions)
  nbInitialProblematicRegions <- length(initialProblematicRegions)

  expect_integer(initialProblematicRegions,
                 any.missing = FALSE, all.missing = FALSE,
                 unique = TRUE, null.ok = TRUE)

  AVAILABLESTATUS <- c("already_feasible", "unresolvable",
                       "partially_resolved", "fully_resolved")
  expect_choice(res$status, AVAILABLESTATUS, null.ok = FALSE)
  expect_identical(res$status, status)


  if (status %in% c("already_feasible", "unresolvable"))
  {
    expect_identical(res$itFusions, 0L)
    expect_identical(res$itTransfers, 0L)
    expect_length(res$resolvedRegions, 0L)
    expect_identical(res$initialSmallRegions, res$finalSmallRegions)
    expect_identical(res$initialBigRegions, res$finalBigRegions)
    expect_identical(res$regionalisation, regionalisation)

    if (status == "already_feasible")
    {
      expect_length(res$initialSmallRegions, 0L)
      expect_length(res$initialBigRegions, 0L)
      expect_null(res$finalSmallRegions)
      expect_null(res$finalBigRegions)
    }
    else if (status == "unresolvable")
    {
      expect_gt(nbInitialProblematicRegions, 0L)
    }
  }
  else
  {
    expect_count(res$itFusions)
    expect_identical(res$itFusions, itFusions)
    expect_count(res$itTransfers)
    expect_identical(res$itTransfers, itTransfers)

    expect_any(res$regionalisation != regionalisation) # nolint: object_usage_linter
    expect_identical(res$regionalisation, newRegionalisation)
    expect_length(res$resolvedRegions, nbResolvedRegions)

    expect_length(res$resolvedRegions, nbResolvedRegions)

    expect_gt(res$itTransfers + res$itFusions, 0L)

    expect_gt(nbInitialProblematicRegions, 0L)

    expect_subset(res$resolvedRegions, initialProblematicRegions)

    if (status == "fully_resolved")
    {
      expect_length(res$resolvedRegions, nbInitialProblematicRegions)
      expect_null(res$finalSmallRegions)
      expect_null(res$finalBigRegions)
    }
    else if (status == "partially_resolved")
    {
      expect_gt(length(res$finalSmallRegions) +
                  length(res$finalBigRegions), 0L)
      expect_subset(res$finalSmallRegions, res$initialSmallRegions,
                    empty.ok = TRUE)
      expect_identical(res$finalSmallRegions, finalSmallRegions)
      expect_subset(res$finalBigRegions, res$initialBigRegions,
                    empty.ok = TRUE)
      expect_identical(res$finalBigRegions, finalBigRegions)

    }

  }


  expect_subset(res$finalSmallRegions, res$initialSmallRegions, empty.ok = TRUE)
  expect_identical(res$finalSmallRegions, finalSmallRegions)
  expect_subset(res$finalBigRegions, res$initialBigRegions, empty.ok = TRUE)
  expect_identical(res$finalBigRegions, finalBigRegions)

}
test_that("No exchanges to be made",
{
  n <- 4L
  m <- 2.0
  sizes <- one_person_per_case(n)
  M <- 2.0

  regionalisation <- c(1L, 1L, 2L, 2L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "already_feasible",
                         itFusions = 0L, itTransfers = 0L,
                         newRegionalisation = regionalisation,
                         nbResolvedRegions = 0L,
                         finalSmallRegions = NULL,
                         finalBigRegions = NULL)

})

test_that("Cannot make any exchanges",
{
  # The first cluster lacks individuals (below min value)
  # and the second has too many (above max value). But if
  # the second gives its transferrable element to the first,
  # then the first cluster will have too many elements.
  n <- 4L
  m <- 2.0
  M <- 3.0
  sizes <- c(1.0, 0.0, 3.0, 3.0)

  regionalisation <- c(1L, 1L, 2L, 2L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "unresolvable",
                         itFusions = 0L, itTransfers = 0L,
                         newRegionalisation = regionalisation,
                         nbResolvedRegions = 0L,
                         finalSmallRegions = 1L,
                         finalBigRegions = 2L)

})

test_that("Impossible regularization",
{
  n <- 4L
  m <- 3.0
  M <- 3.0
  sizes <- c(0.0, 2.0, 2.0, 0.0)


  regionalisation <- c(1L, 1L, 2L, 2L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "unresolvable",
                         itFusions = 0L, itTransfers = 0L,
                         newRegionalisation = regionalisation,
                         nbResolvedRegions = 0L,
                         finalSmallRegions = c(1L, 2L),
                         finalBigRegions = NULL)

})

test_that("Makes a single transfert",
{

  n <- 4L
  m <- 0.0
  M <- 2.0
  sizes <- one_person_per_case(n)

  regionalisation <- c(1L, 1L, 1L, 2L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "fully_resolved",
                         itFusions = 0L, itTransfers = 1L,
                         newRegionalisation = c(1L, 1L, 2L, 2L),
                         nbResolvedRegions = 1L,
                         finalSmallRegions = NULL,
                         finalBigRegions = NULL)

})

test_that("Makes more than one transfert",
{

  n <- 6L
  m <- 2.0
  M <- NOMAXCONSTRAINT
  sizes <- one_person_per_case(n)


  regionalisation <- c(1L, 2L, 2L, 2L, 2L, 3L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "fully_resolved",
                         itFusions = 0L, itTransfers = 2L,
                         newRegionalisation = c(1L, 1L, 2L, 2L, 3L, 3L),
                         nbResolvedRegions = 2L,
                         finalSmallRegions = NULL,
                         finalBigRegions = NULL)

})

test_that("Partial regularization",
{
  # The first cluster is too small while the second is of
  # appropriate size. The 3rd element will be given to the
  # first cluster, but this won't be enough: the only element
  # that can be given next is the 4th, which if removed from
  # the second cluster will result in a size too small for it.
  n <- 5L
  m <- 4.0
  M <- 10.0
  sizes <- c(1.0, 1.0, 1.0, 3.0, 1.0)

  regionalisation <- c(1L, 1L, 2L, 2L, 2L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "partially_resolved",
                         itFusions = 0L, itTransfers = 1L,
                         newRegionalisation = c(1L, 1L, 1L, 2L, 2L),
                         nbResolvedRegions = 0L,
                         finalSmallRegions = 1L,
                         finalBigRegions = NULL)
})

test_that("Singleton transfer allowed",
{
  n <- 3L
  m <- 2.0
  M <- NOMAXCONSTRAINT
  sizes <- c(1.0, 2.0, 1.0)

  regionalisation <- c(1L, 2L, 2L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "fully_resolved",
                         itFusions = 0L, itTransfers = 1L,
                         newRegionalisation = c(1L, 1L, 1L),
                         nbResolvedRegions = 1L,
                         finalSmallRegions = NULL,
                         finalBigRegions = NULL,
                         allowSingletonTransfert = TRUE)

})

test_that("Singleton transfer not allowed",
{
  n <- 3L
  m <- 2.0
  M <- NOMAXCONSTRAINT
  sizes <- c(1.0, 2.0, 1.0)

  regionalisation <- c(1L, 2L, 2L)

  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "unresolvable",
                         itFusions = 0L, itTransfers = 0L,
                         newRegionalisation = c(1L, 2L, 2L),
                         nbResolvedRegions = 0L,
                         finalSmallRegions = 1L,
                         finalBigRegions = NULL,
                         allowSingletonTransfert = FALSE)

})

# Regularization with Max Relaxation

test_that("Impossible without max relaxation, works with",
{
  n <- 3L
  m <- 3.0
  M <- 3.0
  sizes <- c(2.0, 2.0, 0.0)

  regionalisation <- c(1L, 2L, 2L)

  # Without relaxation
  expect_resolver_return(m, M, sizes, regionalisation,
                         status = "unresolvable",
                         itFusions = 0L, itTransfers = 0L,
                         newRegionalisation = regionalisation,
                         nbResolvedRegions = 0L,
                         finalSmallRegions = c(1L, 2L),
                         finalBigRegions = NULL)


  # With relaxation
  expect_resolver_return(m, 4.0 / 3.0 * M, sizes, regionalisation,
                         status = "fully_resolved",
                         itFusions = 0L, itTransfers = 1L,
                         newRegionalisation = c(1L, 1L, 1L),
                         nbResolvedRegions = 2L,
                         finalSmallRegions = NULL,
                         finalBigRegions = NULL)

})
