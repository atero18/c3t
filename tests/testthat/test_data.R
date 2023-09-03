rm(list = ls())


# Normalize a df
test_that("Quantitative variables are normalized",
{
  n <- 10L
  df <- data.frame(quant1 = seq(1L, n),
                   quant2 = seq(n + 1L, 2L * n),
                   qual1 = rep_len(c("a", "b", "c"),
                                    length.out = n))

  dfQuant <- normalize_df(df, standardQuant = TRUE, binarQual = FALSE)

  expect_identical(mean(dfQuant$quant1), 0.0)
  expect_identical( var(dfQuant$quant1), 1.0) # nolint: spaces_inside_linter

  expect_identical(mean(dfQuant$quant2), 0.0)
  expect_identical( var(dfQuant$quant2), 1.0) # nolint: spaces_inside_linter

  expect_identical(dfQuant$qual1, df$qual1)

})

test_that("Categorial variables are splitted",
{
  n <- 10L
  df <- data.frame(quant1 = seq(1L, n),
                   quant2 = seq(n + 1L, 2L * n),
                   qual1 = rep_len(c("a", "b", "c"),
                                   length.out = n),
                   qual2 = rep_len(c("x", "y"),
                                   length.out = n))

  dfQual <- normalize_df(df, standardQuant = FALSE, binarQual = TRUE)
  nbQuant <- 2L
  nbCat1 <- length(unique(df$qual1))
  nbCat2 <- length(unique(df$qual2))

  # Expected number of columns
  expect_identical(ncol(dfQual), nbQuant + nbCat1 - 1L + nbCat2 - 1L)

  # Initial categorial variables must be removed
  expect_false("qual1" %in% colnames(dfQual))
  expect_false("qual2" %in% colnames(dfQual))

  # New variables must be made of 0 and 1
  values <- lapply(dfQual[, -seq_len(nbQuant)], unique)
  values <- unique(do.call("c", values))
  expect_identical(values, c(0L, 1L))

  expect_identical(dfQual[, c("quant1", "quant2")], df[, c("quant1", "quant2")])

})

test_that("Categorical variables are also standardized",
{
  n <- 10L
  df <- data.frame(quant1 = seq(1L, n),
                   quant2 = seq(n + 1L, 2L * n),
                   qual1 = rep_len(c("a", "b", "c"),
                                   length.out = n),
                   qual2 = rep_len(c("x", "y"),
                                   length.out = n))

  dfComp <- normalize_df(df, standardQuant = TRUE, binarQual = TRUE)

  nbQuant <- 2L
  nbCat1 <- length(unique(df$qual1))
  nbCat2 <- length(unique(df$qual2))
  nbCols <- nbQuant + nbCat1 - 1L + nbCat2 - 1L

  # Expected number of columns
  expect_identical(ncol(dfComp), nbCols)


  # Check for centering and scaling
  expect_equal(colMeans(dfComp), double(nbCols), ignore_attr = TRUE)
  expect_equal(vapply(dfComp, var, double(1L)), rep_len(1.0, nbCols),
               ignore_attr = TRUE,
               tolerance = testthat_tolerance())

  # dummy columns must have only two distincts values
  uniqueValuesCat <- lengths(lapply(dfComp[, -seq_len(nbQuant)], unique))
  expect_equal(uniqueValuesCat, rep_len(2L, nbCols - nbQuant),
               ignore_attr = TRUE)

})

test_that("Standardization work when there is only one row",
{
  df <- data.frame(quant1 = 1.0, quant2 = 11.0)

  dfQuant <- normalize_df(df, standardQuant = TRUE, binarQual = TRUE)


  expect_identical(dfQuant, df)
})

# Detection of quantitative variables
test_that("Some are quantitatives, some aren't",
{
  n <- 10L
  df <- data.frame(quant1 = seq(1L, n),
                   quant2 = seq(n + 1L, 2L * n),
                   qual1 = rep_len(c("a", "b", "c"),
                                   length.out = n),
                   qual2 = rep_len(c("x", "y"),
                                   length.out = n))

  mask <- numeric_variabes(df)
  typeVariables <- c(TRUE, TRUE, FALSE, FALSE)
  names(typeVariables) <- colnames(df)

  expect_identical(mask, typeVariables)

  expect_false(all_numeric_variables(df))

})

test_that("All are quantitative",
{
  n <- 10L
  df <- data.frame(quant1 = seq(1L, n),
                   quant2 = seq(n + 1L, 2L * n))

  expect_true(all_numeric_variables(df))
})

# Detect normalized data frame
test_that("With qualitative variables",
{
  n <- 10L
  df <- data.frame(quant1 = seq(1L, n),
                   quant2 = seq(n + 1L, 2L * n),
                   qual1 = rep_len(c("a", "b", "c"),
                                   length.out = n))

  dfQuant <- normalize_df(df, standardQuant = TRUE, binarQual = FALSE)

  expect_false(is_df_normalized(df))
  expect_true(is_df_normalized(dfQuant))
})

test_that("With qualitative variables only",
{
  df <- data.frame(qual1 = rep_len(c("a", "b", "c"),
                                   length.out = 10L),
                   qual2 = rep_len(c("x", "y"),
                                   length.out = 10L))

  expect_true(is_df_normalized(df))
})

test_that("With no data",
{
  df <- data.frame()

  expect_true(is_df_normalized(df))
})
