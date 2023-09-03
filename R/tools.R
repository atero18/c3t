#' Transforms a vector of partitions into a list where each element
#' is a vector indicating which element of the set is present in each cluster.
#'
#' @param partition A numeric vector representing the partition.
#'
#' @returns A list of vectors, each vector representing a cluster with
#' element indices.
#' @keywords internal
partition_to_list <- function(partition)
{
  # Argument verification
  assertPartition(partition)

  split(seq_along(partition), partition)
}

#' Ensures the uniqueness of labels in a vector by appending a numeric index
#' to duplicate labels.
#'
#' @param vec A vector with labels.
#'
#' @returns A vector with unique labels.
#' @noRd
unicite_labels <- function(vec)
{
  labels <- table(vec)
  for (label in names(labels)[labels > 1.0])
  {
    posElementsLabel <- which(vec == label)
    vec[posElementsLabel] <- paste(label,
                                   seq_along(posElementsLabel),
                                   sep = "_")
  }

  vec
}

#' @importFrom stringi stri_trans_general
#' @importFrom checkmate assertCharacter
#' @keywords internal
character_transformation <- function(x,
                                     keepCase = FALSE,
                                     keepAccents = FALSE)
{
  assertCharacter(x)
  if (length(x) == 0L)
    return(character(0L))

  # Set all characters in lower case
  if (!keepCase)
    x <- tolower(x)

  # Remove accents from the characters
  if (!keepAccents)
  {
    transformer <- "Latin-ASCII"
    x <- stri_trans_general(x, transformer)
  }

  x
}

#' @importFrom checkmate assertFlag assertList
match_str <- function(x, choices,
                      checkNames = TRUE,
                      keepCase = FALSE,
                      keepAccents = FALSE)
{

  if (length(x) == 0L)
    return(integer(0L))

  x <- character_transformation(x, keepCase, keepAccents)

  if (is.list(choices))
  {
    assertList(choices, types = "character",
               min.len = 1L,
               any.missing = FALSE,
               all.missing = FALSE)

    assertFlag(checkNames)
    checkNames <- checkNames && !is.null(names(choices))



    if (checkNames)
    {
      matching <-
        match(x, character_transformation(names(choices),
                                          keepCase,
                                          keepAccents))

      naMask <- is.na(matching)
      if (!any(naMask))
        return(matching)
    }
    else
    {
      matching <- rep(NA_integer_, length(x))
      naMask <- rep(TRUE, length(x))
    }

    lengths <- lengths(choices)

    if (any(lengths == 0L))
    {
      stop("All character vectors must contain value")
    }

    k <- 1L

    while (any(naMask) && k <= length(choices))
    {
      mask <- x[naMask] %in% character_transformation(choices[[k]],
                                                      keepCase,
                                                      keepAccents)

      if (any(mask))
      {
        matching[naMask][mask] <- k
        naMask[mask] <- FALSE
      }

      k <- k + 1L
    }
  }
  else
  {
    matching <- match(x, character_transformation(choices,
                                                  keepCase,
                                                  keepAccents))
  }

  matching
}

#' @importFrom checkmate assertFlag
which_na <- function(mat, removeSymmetry = TRUE, removeDiagonal = TRUE)
{
  assertFlag(removeSymmetry)
  assertFlag(removeDiagonal)

  maskNA <- is.na(mat)

  whichNA <- which(maskNA, arr.ind = TRUE)

  if (length(whichNA) == 0L)
    return(whichNA)

  if (removeSymmetry)
  {
    maskSuperior <- whichNA[, 1L] > whichNA[, 2L]
    if (any(maskSuperior))
    {
      whichNA[maskSuperior, 2L:1L] <-
        whichNA[maskSuperior, 1L:2L, drop = FALSE]

      whichNA <- unique(whichNA, MARGIN = 1L)
    }
  }

  if (removeDiagonal)
  {
    maskEquality <- whichNA[, 1L] == whichNA[, 2L]

    if (any(maskEquality))
      whichNA <- whichNA[!maskEquality, drop = FALSE]
  }

  whichNA
}
