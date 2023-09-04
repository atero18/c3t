make_exclusions_lintr <- function()
{
  # Files and folders completly ignores
  ignored <- file.path(c("inst/tutorials",
                         "tests/testthat.R",
                         "tests/testthat/helper.R"))

  exclusions_list <- as.list(ignored)

  testFilenames <- list.files(file.path("tests/testthat"), pattern = "^test_")

  # Indentation linter is ignore within each test file
  for (testFile in testFilenames)
  {
    pathFile <- file.path("tests/testthat/", testFile)
    exclusions_list[[pathFile]] <- list(indentation_linter = Inf)
  }

  exclusions_list
}
