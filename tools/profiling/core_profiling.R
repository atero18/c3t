library(htmlwidgets)
library(profvis)
library(checkmate)
library(devtools)
library(bench)
library(purrr)

devtools::load_all(".")

mark <- partial(bench::mark,
                filter_gc = FALSE,
                relative = TRUE,
                check = TRUE)

sauvegarde_profiling <- function(profiling, prefixe, nomProfiling)
{
  assertClass(profiling, "profvis")
  assertString(prefixe)
  assertString(nomProfiling)

  htmlwidgets::saveWidget(profiling, paste("profiling/saves/", prefixe,
                                           "_", Sys.Date(), "_", nomProfiling,
                                           ".html", sep = ""))

  invisible(TRUE)
}
