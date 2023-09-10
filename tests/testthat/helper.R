
# Load libraries useful for tests and helpers -----------------------------

suppressWarnings(
  {
    library("checkmate", quietly = TRUE)
    loadNamespace("purrr")
  }
)


# Creation of new expectations --------------------------------------------


#' @importFrom checkmate assertFlag assertFunction
#' @importFrom testthat expect_true expect_false
#' @importFrom purrr compose
create_expect_with_logical <- function(withinFunction, expectTrue = TRUE)
{
  assertFlag(expectTrue)
  assertFunction(withinFunction)
  expectFun <- ifelse(expectTrue, testthat::expect_true, testthat::expect_false)
  purrr::compose(expectFun, withinFunction, .dir = "backward")
}

expect_any <- create_expect_with_logical(any, expectTrue = TRUE)

expect_none <- create_expect_with_logical(any, expectTrue = FALSE)

expect_all <- create_expect_with_logical(all, expectTrue = TRUE)

#' @importFrom purrr compose
allXOR <- purrr::compose(all, xor)
expect_xor <- create_expect_with_logical(allXOR, expectTrue = TRUE)

expect_anyNA <- create_expect_with_logical(anyNA, expectTrue = TRUE)

expect_noNA <- create_expect_with_logical(anyNA, expectTrue = FALSE)

expect_all_equal <- create_expect_with_logical(all.equal, expectTrue = TRUE)

expect_symmetric <- create_expect_with_logical(isSymmetric, expectTrue = TRUE)

rm(create_expect_with_logical)

