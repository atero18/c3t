
if (interactive())
{
  pacmanInstalled <-
    suppressMessages(requireNamespace("pacman", quietly = TRUE))

  if (!pacmanInstalled)
  {
    message("Installing `pacman` package")
    utils::install.packages("pacman", repos = "https://cran.rstudio.com/")
  }

  pkgsToInstall <- list()
  pkgsToLoad <- list()

  # Install package dependencies
  if (file.exists("./DESCRIPTION"))
  {
    if (!pacman::p_isinstalled("remotes"))
      utils::install.packages("remotes", repos = "https://cran.rstudio.com/")

    pkgsToInstall[["dependencies"]] <-
      remotes::local_package_deps(".", dependencies = TRUE)
  }

  pkgsToLoad[["dev"]] <-
    c("devtools", "usethis", "roxygen2", "testthat", "covr")

  pkgsToLoad[["prog"]] <- c("checkmate", "lintr")


  pkgsToInstall[["tests"]] <- c("fpc", "withr")

  pkgsToLoad <- do.call("c", pkgsToLoad)
  pkgsToInstall <- do.call("c", pkgsToInstall)

  pkgsToLoad <- unique(pkgsToLoad)
  pkgsToInstall <- setdiff(pkgsToInstall, pkgsToLoad)

  if (length(pkgsToInstall) > 0L)
  {
    notInstalledPkgs <- !vapply(pkgsToInstall, pacman::p_isinstalled, logical(1L))
    pkgsToInstall <- pkgsToInstall[notInstalledPkgs]
  }


  if (length(pkgsToInstall) > 0L)
    utils::install.packages(pkgsToInstall, repos = "https://cran.rstudio.com/")

  if (length(pkgsToLoad) > 0L)
    pacman::p_load(char = pkgsToLoad, character.only = TRUE)

  rm(list = ls())
}

usethis::use_tidy_description()

options(c3t_verbose = TRUE)
