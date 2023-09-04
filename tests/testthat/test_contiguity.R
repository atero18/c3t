rm(list = ls())


# Checking the contiguity matrix for a simple grid
# 2x2 Grid
x <- y <-  2L
# Rook Contiguity
mat_rook <- grid_contiguity_matrix(x, y, "Rook")
testContiguityMatrix(mat_rook)
# nolint start: T_and_F_symbol
true_mat_rook <- matrix(c(T, T, T, F,
                          T, T, F, T,
                          T, F, T, T,
                          F, T, T, T), ncol = 4L, nrow = 4L, byrow = TRUE)
# nolint end
test_that("Rook 2*2",
{
  expect_true(all(mat_rook == true_mat_rook))
})

# Queen Contiguity
# For "Queen" case, every cell is contiguous to every other cell

mat_queen <- grid_contiguity_matrix(x, y, "Queen")
testContiguityMatrix(mat_queen)
test_that("Queen 2*2",
{
  expect_true(all(mat_queen))
})

# 3x2 Grid
x <- 3L
y <- 2L
# Rook Contiguity
mat_rook <- grid_contiguity_matrix(x, y, "Rook")
testContiguityMatrix(mat_rook)
# nolint start: T_and_F_symbol
true_mat_rook <- matrix(c(T, T, F, T, F, F,
                          T, T, T, F, T, F,
                          F, T, T, F, F, T,
                          T, F, F, T, T, F,
                          F, T, F, T, T, T,
                          F, F, T, F, T, T),
                        ncol = 6L, nrow = 6L, byrow = TRUE)
# nolint end
test_that("Rook 3*2",
{
  expect_true(all(mat_rook == true_mat_rook))
})

# Queen Contiguity
mat_queen <- grid_contiguity_matrix(x, y, "Queen")
testContiguityMatrix(mat_queen)
# nolint start: T_and_F_symbol
true_mat_queen <- matrix(c(T, T, F, T, T, F,
                           T, T, T, T, T, T,
                           F, T, T, F, T, T,
                           T, T, F, T, T, F,
                           T, T, T, T, T, T,
                           F, T, T, F, T, T),
                         ncol = 6L, nrow = 6L, byrow = TRUE)
# nolint end
test_that("Queen 3*2",
{
  expect_true(all(mat_queen == true_mat_queen))
})

# Coordinate Conversion Tool (x,y) -> Matrix Position
test_that("Coordinate Conversion -> Rank",
{
  expect_identical(xy_to_rank_grid(1L, 1L, 2L, 3L), 1L)
  expect_identical(xy_to_rank_grid(2L, 1L, 2L, 3L), 2L)
  expect_identical(xy_to_rank_grid(1L, 2L, 2L, 3L), 3L)
})

# Cells (1,1) and (2,1) are indeed contiguous in "Rook" and "Queen" sense
test_that("Contiguity (1,1) (2,1)",
{
  expect_true(true_mat_rook[xy_to_rank_grid(1L, 1L, 3L, 2L),
                            xy_to_rank_grid(2L, 1L, 3L, 2L)])
  expect_true(true_mat_queen[xy_to_rank_grid(1L, 1L, 3L, 2L),
                             xy_to_rank_grid(2L, 1L, 3L, 2L)])
})

# Cells (1,1) and (2,2) are contiguous in "Queen" sense but not "Rook"
test_that("Partial Contiguity (1,1) (2,2)",
{
  expect_false(true_mat_rook[xy_to_rank_grid(1L, 1L, 3L, 2L),
                             xy_to_rank_grid(2L, 2L, 3L, 2L)])
  expect_true(true_mat_queen[xy_to_rank_grid(1L, 1L, 3L, 2L),
                             xy_to_rank_grid(2L, 2L, 3L, 2L)])
})
# Cells (1,1) and (3,2) are not contiguous in "Rook" and "Queen" sense
test_that("Non-Contiguity (1,1) (3,2)",
{
  expect_false(true_mat_rook[xy_to_rank_grid(1L, 1L, 3L, 2L),
                             xy_to_rank_grid(3L, 2L, 3L, 2L)])
  expect_false(true_mat_queen[xy_to_rank_grid(1L, 1L, 3L, 2L),
                              xy_to_rank_grid(3L, 2L, 3L, 2L)])
})


# Matrix Position -> Coordinate Conversion Tool
test_that("Rank Coordinate Conversion, 1st Row",
{
  res <- rank_to_xy_grid(3L, 5L, 2L)
  expect_identical(res$i, 3L)
  expect_identical(res$j, 1L)
})

test_that("Rank Coordinate Conversion, End of Row",
{
  res <- rank_to_xy_grid(5L, 5L, 2L)
  expect_identical(res$i, 5L)
  expect_identical(res$j, 1L)
})

test_that("Rank Coordinate Conversion, New Row",
{
  res <- rank_to_xy_grid(9L, 5L, 2L)
  expect_identical(res$i, 4L)
  expect_identical(res$j, 2L)
})

# Checking Inter-Cluster Contiguity Detection

test_that("Cluster Contiguity",
{
  x <- 2L
  y <- 3L
  contiguity_mat <- grid_contiguity_matrix(x, y, "Queen")
  partition_vec <- c(1L, 1L, 2L, 2L, 3L, 3L)
  M <- clusters_contiguity_matrix(partition_vec, contiguity_mat)

  expect_true(testContiguityMatrix(M, isComplete = TRUE))
  expect_false(M[1L, 3L])
  expect_identical(sum(M), 4L)
})

# Checking Internal Boundary Calculation
test_that("Internal Boundary",
{
  # nolint start: T_and_F_symbol_linter
  M <- matrix(c(F, T, T, F, F,
                T, F, T, F, F,
                T, T, F, T, F,
                F, F, T, F, T,
                F, F, F, T, F), nrow = 5L, byrow = TRUE)
  # nolint end
  region_vec <- c(1L, 1L, 1L, 2L, 2L)
  boundary <- interior_boundary(M, region_vec, removeSymmetry = TRUE)
  expect_identical(nrow(boundary), 1L)
  boundarySym <- interior_boundary(M, region_vec, removeSymmetry = FALSE)
  expect_identical(nrow(boundarySym), 2L)

  expect_true(all(boundarySym[1L, ] == c(3L, 4L, 1L, 2L)) ||
                all(boundarySym[1L, ] == c(4L, 3L, 2L, 1L)))
  expect_true(all(boundarySym[2L, ] == c(3L, 4L, 1L, 2L)) ||
                all(boundarySym[2L, ] == c(4L, 3L, 2L, 1L)))

})

# Checking Articulation Point Detection
test_that("Articulation Points (Single Present)",
{
  # nolint start: T_and_F_symbol_linter
  M <- matrix(c(F, T, T, F, F,
                T, F, T, F, F,
                T, T, F, T, T,
                F, F, T, F, T,
                F, F, T, T, F), nrow = 5L, byrow = 5L)
  # nolint end

  set.seed(123L)
  for (i in setdiff(1L:5L, 3L))
    expect_false(is_articulation_pt(i, M))

  expect_true(is_articulation_pt(3L, M))
})


test_that("Articulation Points (None Present)",
{
  # nolint start: T_and_F_symbol_linter
  M <- matrix(c(F, T, F, F, T,
                T, F, T, F, F,
                F, T, F, T, F,
                F, F, T, F, T,
                T, F, F, T, F), nrow = 5L, byrow = 5L)
  # nolint end

  set.seed(123L)
  for (i in 1L:5L)
    expect_false(is_articulation_pt(i, M))
})

test_that("Articulation Points - Test with Clustering",
{
  # nolint start: T_and_F_symbol_linter
  M <- matrix(c(F, T, F, F, F,
                T, F, T, T, F,
                F, T, F, F, F,
                F, T, F, F, T,
                F, F, F, T, F), nrow = 5L, byrow = TRUE)
  # nolint end

  region_vec <- c(1L, 1L, 1L, 2L, 2L)
  articulation_point <- is_articulation_pt_cluster(2L, M, region_vec)
  expect_true(articulation_point)

})

# Checking Transferable Point Detection
# Checking Internal Boundary Calculation
test_that(paste0("Transferable Points - Transferable Points are the Points on ",
                 "the Internal Boundary"),
{
  # nolint start: T_and_F_symbol_linter
  M <- matrix(c(F, T, T, F, F,
                T, F, T, F, F,
                T, T, F, T, F,
                F, F, T, F, T,
                F, F, F, T, F), nrow = 5L, byrow = TRUE)
  # nolint end
  contiguityGraph <- contiguity_matrix_to_graph(M)

  region_vec <- c(1L, 1L, 1L, 2L, 2L)
  transferable_elements <-
    transferable_elements_pb(region_vec, contiguityGraph = contiguityGraph)

  expect_length(transferable_elements$articulation, 0L)
  expect_setequal(transferable_elements$transferable, c(3L, 4L))
  expect_true(all(transferable_elements$frontier[1L, ] == c(3L, 4L, 1L, 2L)) ||
                all(transferable_elements$frontier[1L, ] == c(4L, 3L, 2L, 1L)))

  expect_true(all(transferable_elements$frontier[2L, ] == c(3L, 4L, 1L, 2L)) ||
                all(transferable_elements$frontier[2L, ] == c(4L, 3L, 2L, 1L)))
})

test_that("Transferable Points - Loss of Symmetry",
{
  # nolint start: T_and_F_symbol_linter
  M <- matrix(c(F, T, F, F, F,
                T, F, T, T, F,
                F, T, F, F, F,
                F, T, F, F, T,
                F, F, F, T, F), nrow = 5L, byrow = TRUE)
  # nolint end
  contiguityGraph <- contiguity_matrix_to_graph(M)

  region_vec <- c(1L, 1L, 1L, 2L, 2L)
  transferable_elements <-
    transferable_elements_pb(region_vec, contiguityGraph = contiguityGraph)

  expect_identical(nrow(transferable_elements$frontier), 2L)
  expect_identical(transferable_elements$articulation, 2L)
  expect_identical(transferable_elements$transferable, 4L)
  expect_true(all(transferable_elements$frontier[1L, ] == c(2L, 4L, 1L, 2L)) ||
                all(transferable_elements$frontier[1L, ] == c(4L, 2L, 2L, 1L)))

  expect_true(all(transferable_elements$frontier[2L, ] == c(2L, 4L, 1L, 2L)) ||
                all(transferable_elements$frontier[2L, ] == c(4L, 2L, 2L, 1L)))
})
