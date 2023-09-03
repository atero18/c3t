#' Normalization of data
#' @param df A data frame of equivalent ([tibble][tibble::tibble-package],
#' data table...)
#' @name normalize_df
NULL

#' @describeIn normalize_df Normalize quantitative and / or qualitative `df`
#' variables.
#' @param standardQuant `TRUE` (default) if quantitative variables of `df` must
#' be normalised / standardised, i.e centered (having a 0 empirical mean)
#' and reduced / scaled (having a 1 unbiased empirical standard deviation).
#' `FALSE` otherwise. If the number of rows is equal to 1 standardisation isn't
#' realized. Standardisation is realized after transformation
#' from qualitative to quantitative, so if `binarQual` is equal to `TRUE`
#' transformed qualitative variables will be normalized.
#' @param binarQual `TRUE` if non-numeric values must be converted in numeric
#' form by creating for each qualitative variable one dummy variable
#' per modality minus one. Original columns will be removed.
#' `FALSE` (default) otherwise. Ignore if all variables
#' are quantitative. If `standardQuant` equals to `TRUE` created dummies
#' will be standardised as well.
#' @returns For `normalize_df` the normalized version of `df`. If
#' `binarQual = TRUE` and some categorical variables were if `df` the number
#' of variables (columns) of `df` might be higher.
#' @importFrom fastDummies dummy_cols
#' @seealso [base::scale()]
#' @seealso [fastDummies::dummy_cols()]
#' @references Kaplan, J. & Schlegel, B. (2023).
#' [fastDummies](https://github.com/jacobkap/fastDummies): Fast Creation
#' of Dummy (Binary) Columns and Rows from Categorical Variables.
#' @export
#' @importFrom checkmate assertFlag assertDataFrame
normalize_df <- function(df, standardQuant = TRUE, binarQual = FALSE)
{
  # Check arguments
  assertFlag(standardQuant)
  assertFlag(binarQual)
  assertDataFrame(df)

  if (length(df) == 0L)
    return(df)

  numericCols <- numeric_variabes(df)

  binarQual <- binarQual && !all(numericCols)

  if (binarQual)
  {
    df <- dummy_cols(df, remove_first_dummy = TRUE, ignore_na = TRUE,
                     remove_selected_columns = TRUE)


    numericCols <- numeric_variabes(df)
  }

  standardQuant <- standardQuant && nrow(df) > 1L && any(numericCols)
  if (standardQuant)
    df[, numericCols] <- scale(df[, numericCols], center = TRUE, scale = TRUE)

  return(df)
}

#' @describeIn normalize_df Give a logical mask indicating which columns
#' of `df` are numeric.
#' @returns For `numeric_variables` a logical vector of length `ncol(df)`
#' indicating which variable is numeric and which is not.
#' @export
#' @importFrom checkmate assertDataFrame
numeric_variabes <- function(df)
{
  # Check arguments
  assertDataFrame(df)

  vapply(df, is.numeric, logical(1L))
}

#' @describeIn normalize_df Give a flag indicating if all variables of `df`
#' are quantitative. Shortcut for `all(numeric_variables(df))`.
#' @returns For `all_numeric_variables`, `TRUE` if all variables of `df`
#' are quantitative, `FALSE` otherwise. Returns `TRUE` by convention if there
#' is no variable at all.
#' @export
all_numeric_variables <- function(df) all(numeric_variabes(df))

#' @describeIn normalize_df Check if the numeric variables of `df`
#' are normalized.
#' @returns For `is_df_normalized`, a flag equals to `TRUE` if all numeric
#' variables are normalized, `FALSE` otherwise. Returns `TRUE` if there is
#' no numeric variables.
#' @export
#' @importFrom stats var
#' @importFrom checkmate assertDataFrame
is_df_normalized <- function(df)
{
  # Check arguments
  assertDataFrame(df)

  if (length(df) == 0L)
    return(TRUE)

  numCols <- vapply(df, is.numeric, logical(1L))

  if (!any(numCols))
    return(TRUE)

  if (nrow(df) == 1L)
    return(all(colMeans(df[, numCols]) == 0.0))


  all(colMeans(df[, numCols]) == 0.0) &&
    all(vapply(df[, numCols], var, numeric(1L)) == 1L)
}
