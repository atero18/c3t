make_exclusions_lintr <- function()
{
  # Files and folders completly ignores
  ignored <- c("inst/tutorials", "tests/testthat.R", "tests/testthat/helper.R")

  exclusions_list <- as.list(ignored)

  testFilenames <- list.files("tests/testthat", pattern = "^test_")

  # Indentation linter is ignore within each test file
  for (testFile in testFilenames)
  {
    pathFile <- paste0("tests/testthat/", testFile)
    exclusions_list[[pathFile]] <- list(indentation_linter = Inf)
  }

  exclusions_list
}


