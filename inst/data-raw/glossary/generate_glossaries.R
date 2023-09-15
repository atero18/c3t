RELATIVEFOLDERPATH <- file.path("./inst/data-raw/glossary")

#' @importFrom checkmate assertString assertChoice assertFileExists
#' @importFrom dplyr arrange
#' @importFrom stringr str_to_title str_to_sentence
#' @importFrom usethis use_vignette
generate_glossary <- function(lang, definitionsFile = file.path(RELATIVEFOLDERPATH, "glossary.csv"))
{

  assertString(lang)
  lang <- tolower(lang)

  assertFileExists(definitionsFile)
  metadataGlossaries <-
    read.csv2(file.path(RELATIVEFOLDERPATH, "glossaries_metadata.csv"),
             header = TRUE)
  AVAILABLELANG <- metadataGlossaries$language
  assertChoice(lang, AVAILABLELANG)

  metadataGlossary <- metadataGlossaries[metadataGlossaries == lang, ]
  rm(metadataGlossaries)

  definitions <- read.csv2(definitionsFile, header = TRUE)
  definitions <- definitions[definitions$language %in% c(lang, "*", NA), ]

  if (nrow(definitions) == 0L)
    stop("No definitions available for this language")

  definitions$word <- str_to_title(definitions$word)
  definitions$definition <- str_to_sentence(definitions$definition)
  definitions <- arrange(definitions, word)
  definitions$combine = paste0("**", definitions$word, "**", " : ", definitions$definition)

  definitions$firstLetter <- vapply(definitions$word, substr, character(1L), 1L, 1L)

  title <- metadataGlossary$title

  nameVignette <- paste("glossary", lang, sep = "-")


  usethis::use_vignette(nameVignette, title)

  pathVignette <- paste0("vignettes/", nameVignette, ".Rmd")

  glossaryText <- readLines(pathVignette)

  add_line <- function(text, pos = 0L, listItem = TRUE)
  {
    listItem == listItem && pos == 0L
    if (pos > 0L)
      text <- paste(rep("#", pos), text, "\n")
    else if (listItem)
      text <- paste("-", text)

    glossaryText <<- c(glossaryText, text)
  }

  intro <- metadataGlossary$intro
  add_line(intro)

  for (firstLetter in unique(definitions$firstLetter))
  {
    add_line(firstLetter, pos = 1L)

    whichWords <- which(definitions$firstLetter == firstLetter)
    for (line in definitions[whichWords, "combine"])
      add_line(line)

  }

  writeLines(glossaryText, pathVignette, sep = "\n")

}
