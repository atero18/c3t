rm(list = ls())


# Test of update with no NA
#
test_that("Update single linkage",
{
  # nolint start: implicit_integer_linter
  linkDist <- matrix(c(0, 1, 2, 3, 4,
                       1, 0, 4, 3, 5,
                       2, 4, 0, 4, 4,
                       3, 3, 4, 0, 2,
                       4, 5, 4, 2, 0), nrow = 5L)
  # nolint end

  c1 <- 3L
  c2 <- 5L
  updateLink <- update_single_linkage(linkDist, c1, c2)
  expect_true(is.matrix(updateLink))
  expect_identical(nrow(updateLink), nrow(linkDist))
  expect_true(isSymmetric(updateLink))
  expect_true(all(updateLink >= 0.0, na.rm = TRUE))
  idsDiffc1c2 <- setdiff(seq_len(nrow(linkDist)), c(c1, c2))
  expect_identical(linkDist[idsDiffc1c2, idsDiffc1c2],
                   updateLink[idsDiffc1c2, idsDiffc1c2])
  idsDiffc1 <- setdiff(seq_len(nrow(linkDist)), c1)
  idsDiffc2 <- setdiff(seq_len(nrow(linkDist)), c2)
  expect_true(all(is.nan(updateLink[c2, idsDiffc2])))
  expect_identical(updateLink[c1, idsDiffc1c2],
                   pmin(linkDist[c1, idsDiffc1c2], linkDist[c2, idsDiffc1c2]))

})

test_that("Update complete linkage",
{
  # nolint start: implicit_integer_linter
  linkDist <- matrix(c(0, 1, 2, 3, 4,
                       1, 0, 4, 3, 5,
                       2, 4, 0, 4, 4,
                       3, 3, 4, 0, 2,
                       4, 5, 4, 2, 0), nrow = 5L)
  # nolint end

  c1 <- 3L
  c2 <- 5L
  updateLink <- update_complete_linkage(linkDist, c1, c2)
  expect_true(is.matrix(updateLink))
  expect_identical(nrow(updateLink), nrow(linkDist))
  expect_true(isSymmetric(updateLink))
  expect_true(all(updateLink >= 0.0, na.rm = TRUE))
  idsDiffc1c2 <- setdiff(seq_len(nrow(linkDist)), c(c1, c2))
  expect_identical(linkDist[idsDiffc1c2, idsDiffc1c2],
                   updateLink[idsDiffc1c2, idsDiffc1c2])
  idsDiffc1 <- setdiff(seq_len(nrow(linkDist)), c1)
  idsDiffc2 <- setdiff(seq_len(nrow(linkDist)), c2)
  expect_true(all(is.nan(updateLink[c2, idsDiffc2])))
  expect_identical(updateLink[c1, idsDiffc1c2],
                   pmax(linkDist[c1, idsDiffc1c2], linkDist[c2, idsDiffc1c2]))

})


test_that("Update average linkage",
{
  # nolint start: implicit_integer_linter
  linkDist <- matrix(c(0, 1, 2, 3, 4,
                       1, 0, 4, 3, 5,
                       2, 4, 0, 4, 4,
                       3, 3, 4, 0, 2,
                       4, 5, 4, 2, 0), nrow = 5L)
  # nolint end


 linkDist <- linkDist[1L:3L, 1L:3L]
 partition <- c(1L, 2L, 1L, 2L, 3L)
 nbElems <- unname(table(partition))

 for (c1 in 1L:3L)
 {
   c2 <- c1 + 1L
   while (c2 <= 3L)
   {
     updateLink <- update_average_linkage(linkDist,
                                          c1, c2,
                                          partitionBefore = partition)
     expect_true(is.matrix(updateLink))
     expect_identical(nrow(updateLink), nrow(linkDist))
     expect_true(isSymmetric(updateLink))
     expect_true(all(updateLink >= 0.0, na.rm = TRUE))
     idsDiffc1c2 <- setdiff(seq_len(nrow(linkDist)), c(c1, c2))
     expect_identical(linkDist[idsDiffc1c2, idsDiffc1c2],
                      updateLink[idsDiffc1c2, idsDiffc1c2])
     idsDiffc1 <- setdiff(seq_len(nrow(linkDist)), c1)
     idsDiffc2 <- setdiff(seq_len(nrow(linkDist)), c2)
     expect_true(all(is.nan(updateLink[c2, idsDiffc2])))
     expect_identical(updateLink[c1, idsDiffc1c2],
                      nbElems[c1] / (nbElems[c1] + nbElems[c2]) *
                        linkDist[c1, idsDiffc1c2] +
                        nbElems[c2] / (nbElems[c1] + nbElems[c2]) *
                        linkDist[c2, idsDiffc1c2])
     c2 <- c2 + 1L
   }

 }

})

test_that("Update some Lance-Williams linkage",
{

  # nolint start: implicit_integer_linter
  linkDist <- matrix(c(0, 1, 2, 3, 4,
                       1, 0, 4, 3, 5,
                       2, 4, 0, 4, 4,
                       3, 3, 4, 0, 2,
                       4, 5, 4, 2, 0), nrow = 5L)
  # nolint end

  c1 <- 2L
  c2 <- 4L

  a1 <- 3.0
  a2 <- 5.3
  b <- 1.3
  g <- 2L

  updateLink <- lance_williams_update(linkDist, c1, c2, a1, a2, b, g)
  expect_true(is.matrix(updateLink))
  expect_identical(nrow(updateLink), nrow(linkDist))
  expect_true(isSymmetric(updateLink))
  expect_true(all(updateLink >= 0.0, na.rm = TRUE))
  idsDiffc1c2 <- setdiff(seq_len(nrow(linkDist)), c(c1, c2))
  expect_identical(linkDist[idsDiffc1c2, idsDiffc1c2],
                   updateLink[idsDiffc1c2, idsDiffc1c2])
  idsDiffc1 <- setdiff(seq_len(nrow(linkDist)), c1)
  idsDiffc2 <- setdiff(seq_len(nrow(linkDist)), c2)
  expect_true(all(is.nan(updateLink[c2, idsDiffc2])))
  expect_identical(updateLink[c1, idsDiffc1c2],
                   a1 * linkDist[c1, idsDiffc1c2] +
                   a2 * linkDist[c2, idsDiffc1c2] +
                    b * linkDist[c1, c2] +
                    g * abs(linkDist[c1, idsDiffc1c2] -
                            linkDist[c2, idsDiffc1c2]))

})

# Specify useless clusters
##
# Calculate missing distances when needed
test_that("Calculate distance when needed - single linkage",
{
  # nolint start: implicit_integer_linter
  elemDist <- matrix(c(0, 1, 2, 3, 4,
                       1, 0, 4, 3, 5,
                       2, 4, 0, 4, 4,
                       3, 3, 4, 0, 2,
                       4, 5, 4, 2, 0), nrow = 5L)
  # nolint end

  linkDist <- elemDist[1L:3L, 1L:3L]

  c1 <- 2L
  c2 <- 3L

  linkDist[1L, 2L] <- linkDist[2L, 1L] <- NA_real_

  partition <- c(1L, 2L, 3L, 2L, 1L)
  updateLink <- update_single_linkage(linkDist, c1, c2,
                                      elemDistances = elemDist,
                                      partitionBefore = partition)

  expect_true(is.matrix(updateLink))
  expect_identical(nrow(updateLink), nrow(linkDist))
  expect_true(isSymmetric(updateLink))
  expect_true(all(updateLink >= 0.0, na.rm = TRUE))
  expect_false(is.na(updateLink[1L, 2L]))
  expect_identical(updateLink[1L, 2L],
                   min(elemDist[which(partition == 1L),
                                which(partition %in% c(c1, c2))]))
})

test_that("Don't calcule distances if not necessary - single linkage",
{
  # nolint start: implicit_integer_linter
  elemDist <- matrix(c(0, 1, 2, 3, 4,
                       1, 0, 4, 3, 5,
                       2, 4, 0, 4, 4,
                       3, 3, 4, 0, 2,
                       4, 5, 4, 2, 0), nrow = 5L)
  # nolint end

  c1 <- 1L
  c2 <- 2L
  linkDist <- elemDist
  linkDist[3L, 4L] <- linkDist[4L, 3L] <- NA_real_
  partition <- seq_len(5L)

  updateLink <- update_single_linkage(linkDist, c1, c2,
                                      elemDistances = elemDist,
                                      partitionBefore = partition)

  expect_true(isSymmetric(updateLink))
  expect_true(is.na(updateLink[3L, 4L]))
})
