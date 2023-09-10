#' @importFrom checkmate assertCount assertNumber
gen_distances <- function(n, mean = 1.0, sd = 1.0)
{
  assertCount(n, positive = TRUE)
  assertNumber(mean)
  assertNumber(sd, lower = 0.0, finite = TRUE)

  nbDistances <- (n - 1L) * n / 2L

  distances_vec <- rnorm(nbDistances, mean = mean, sd = sd)
  distances_vec <- abs(distances_vec)

  distances_mat <- matrix(0.0, nrow = n, ncol = n)

  distances_mat[lower.tri(distances_mat, diag = FALSE)] <- distances_vec

  distances_mat + t(distances_mat)
}

#' @importFrom checkmate assertCount assertFlag
#' @keywords internal
empty_distance_matrix <- function(n, madeofNA = FALSE)
{
  assertCount(n)
  assertFlag(madeofNA)

  value <- ifelse(madeofNA, NA_real_, 0.0)

  matrix(value, nrow = n, ncol = n)
}
