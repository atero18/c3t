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

mark <- purrr::partial(bench::mark,
                       filter_gc = FALSE,
                       relative = TRUE,
                       check = TRUE)

press <- purrr::partial(bench::press,
                        filter_gc = FALSE,
                        check = FALSE)

SAVEPROFVIS <- FALSE
save_profiling <- function(profiling, prefixe, nomProfiling)
{
  checkmate::assertClass(profiling, "profvis")
  checkmate::assertString(prefixe)
  checkmate::assertString(nomProfiling)

  SAVEFOLDER <- file.path("./tools/profiling/saves/")

  if (!dir.exists(SAVEFOLDER))
    dir.create(SAVEFOLDER)

  htmlwidgets::saveWidget(
    profiling,
    paste0("tools/profiling/saves/", prefixe, "_", Sys.Date(), "_",
              nomProfiling, ".html"),
    selfcontained = TRUE)

  invisible(profiling)
}

profiling <- function(..., save = SAVEPROFVIS, prefix, name, seed = 123L)
{
  checkmate::assertFlag(save)
  checkmate::assertInt(seed, null.ok = TRUE)

  if (!is.null(seed))
    set.seed(seed)


  if (save)
  {
    checkmate::assertString(prefix, min.chars = 1L)
    checkmate::assertString(name, min.chars = 1L)
  }

  p <- profvis::profvis(...)

  if (save)
    save_profiling(p, prefix, name)

  invisible(p)
}

set.seed(123L)
