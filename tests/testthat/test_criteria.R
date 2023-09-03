rm(list = ls())

# Optimal partitions
test_that("optimal partitions gives the right order",
{
  values <- c(0.3, 0.8, 0.2)
  names(values) <- c("x", "y", "z")
  rightOrder <- c(2L, 1L, 3L)
  rightNames <- names(values)[rightOrder]

  # When byName = FALSE
  dunn <- optimal_partitions(values, "Dunn", byName = FALSE)
  expect_identical(dunn, rightOrder)
  ch <- optimal_partitions(values, "CHI", byName = FALSE)
  expect_identical(ch, rightOrder)

  if ("Silhouette" %in% available_criteria())
  {
    silhouette <- optimal_partitions(values, "Silhouette", byName = FALSE)
    expect_identical(silhouette, rightOrder)
  }

  # When byName = TRUE
  dunn <- optimal_partitions(values, "Dunn", byName = TRUE)
  expect_identical(dunn, rightNames)
  ch <- optimal_partitions(values, "CHI", byName = TRUE)
  expect_identical(ch, rightNames)

  if ("Silhouette" %in% available_criteria())
  {
    silhouette <- optimal_partitions(values, "Silhouette", byName = TRUE)
    expect_identical(silhouette, rightNames)
  }
})
