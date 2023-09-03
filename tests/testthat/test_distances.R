rm(list = ls())

# Jump distance
# 2 clusters of 3 elements with a context of 2 variables
# One element is present in both clusters
dataTest <- data.frame(q1 = c(3.5, 2.6, 1.7, -1.2, 1.7, 3.9),
                       q2 = c(1.0, 0.0, 6.0,  5.0, 6.0, 4.0)) # nolint: spaces_inside_linter


partition <- c(1L, 1L, 1L, 2L, 2L, 2L)
distance <- function(i, j) {norm(dataTest[i, ] - dataTest[j, ], type = "2")}
distances_mat <- matrix(0.0, nrow = 6L, ncol = 6L)
for (i in 1L:5L)
{
  for (j in (i + 1L):6L)
    distances_mat[i, j] <- distances_mat[j, i] <- distance(i, j)
}

# Since an element is present in both clusters, the expected value is 0
test_that("Single linkage",
{
  expect_identical(single_linkage(distances_mat[1L:3L, ][, 4L:6L]), 0.0)
})

# For max jump distance, the value should be > 0 as the sets are not equal
test_that("Complete linkage",
{
  expect_gt(complete_linkage(distances_mat[1L:3L, ][, 4L:6L]), 0.0)
  expect_identical(complete_linkage(distances_mat[1L:3L, ][, 4L:6L]),
                   max(distances_mat[1L:3L, ][, 4L:6L])) # nolint: spaces_inside_linter
})
