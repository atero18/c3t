rm(list = ls())
loadNamespace("devtools")
loadNamespace("purrr")
loadNamespace("htmlwidgets")
loadNamespace("bench")
loadNamespace("profvis")
loadNamespace("checkmate")

exportAllFunction <- TRUE
devtools::load_all(".", reset = TRUE, export_all = exportAllFunction)

rm(list = ls())

profvis <- profvis::profvis

mark <- purrr::partial(bench::mark,
                       filter_gc = FALSE,
                       relative = TRUE,
                       check = TRUE)

press <- purrr::partial(bench::press,
                        filter_gc = FALSE,
                        check = FALSE)

save_profiling <- function(profiling, prefixe, nomProfiling)
{
  assertClass(profiling, "profvis")
  assertString(prefixe)
  assertString(nomProfiling)

  htmlwidgets::saveWidget(
    profiling,
    file.path("tools/profiling/saves/", prefixe, "_", Sys.Date(), "_",
              nomProfiling, ".html"))

  invisible(TRUE)
}

set.seed(123L)
