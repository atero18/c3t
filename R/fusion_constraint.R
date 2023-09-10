
# Find information about fusion constraints -------------------------------


#' Available fusion constraints
#'
#' What parameters are available for fusion constraints?
#' @name available_fusion_constraints
#' @returns a character vector with available choices.
NULL

#' @describeIn available_fusion_constraints return implemented
#' fusion constraints.
#' @family available parameters
#' @keywords internal
#' @export
available_fusion_constraints <- function()
{
  c(NA_character_, FUSIONCONSTRAINTS$fusionConstraint)
}

#' @importFrom checkmate assertFlag
corresponding_fusion_const <- function(fusionConstraints, simplify = TRUE)
{
  if (is.null(fusionConstraints))
    return(NA_character_)

  else if (length(fusionConstraints) == 0L)
    return(character(0L))

  assertFlag(simplify)

  fusionConstraints <-
    ifelse(fusionConstraints == FALSE, NA_character_, fusionConstraints) # nolint: redundant_equals_linter

  if (all(is.na(fusionConstraints)))
  {
    if (simplify)
      return(NA_character_)
    else
      return(fusionConstraints)

  }

  naMask <- is.na(fusionConstraints)

  if (any(naMask))
  {
    correspondingNA <- NA_character_
    fusionConstraints <- fusionConstraints[!is.na(fusionConstraints)]
  }

  else
    correspondingNA <- NULL

  correspondingString <-
    FUSIONCONSTRAINTS$fusionConstraint[match_str(fusionConstraints,
                                                 FUSIONCONSTRAINTS$names)]

  if (anyNA(correspondingString))
    correspondingString[is.na(correspondingString)] <- NaN

  if (simplify)
    correspondingString <- unique(correspondingString)

  c(correspondingNA, correspondingString)
}

simplify_fusion_constraints <- function(fusionConstraints, pb = NULL,
                                        partition = NULL)
{
  fusionConstraints <-
    corresponding_fusion_const(fusionConstraints, simplify = TRUE)

  uselessMode <-
    FUSIONCONSTRAINTS$fusionConstraint[FUSIONCONSTRAINTS$needsMinConstraint]
  if (inherits(pb, "pbCon") &&
      any(fusionConstraints %in% uselessMode) &&
      (!pb$hasMinConstraint() ||
       (!is.null(partition) && all(.clusters_sizes(partition, pb$sizes) >=
                                   pb$m))))
  {
    fusionConstraints <- setdiff(fusionConstraints, uselessMode)

  }

  fusionConstraints
}

#' @importFrom checkmate assertFlag assertVector
assertFusionConstraints <- function(fusionConstraints, unique = FALSE)
{
  assertFlag(unique)
  if (unique)
    max.len <- 1L
  else
    max.len <- NULL
  assertVector(fusionConstraints, min.len = 1L,
               max.len = max.len,
               any.missing = TRUE, all.missing = TRUE)

  fusionConstraints <-
    corresponding_fusion_const(fusionConstraints, simplify = FALSE)

  if (any(is.nan(fusionConstraints)))
  {
    stop("Some fusion constraint have not been recognized")
  }

  TRUE
}

assertFusionConstraint <-
  function(fusionConstraint) assertFusionConstraints(fusionConstraint,
                                                     unique = TRUE)


#' @describeIn available_fusion_constraints return implemented
#' fusion constraint modes.
#' @family available parameters
#' @keywords internal
#' @export
available_fusion_modes <- function()
{
  FUSIONCONSTRAINTMODES
}

#' @importFrom checkmate assertFlag
corresponding_fusion_mode <- function(modes, simplify = TRUE)
{
  if (is.null(modes))
    return(NA_character_)

  else if (length(modes) == 0L)
    return(character(0L))

  assertFlag(simplify)

  modes <- ifelse(modes == FALSE, NA_character_, modes) # nolint: redundant_equals_linter

  if (all(is.na(modes)))
  {
    if (simplify)
      return(NA_character_)
    else
      return(modes)

  }

  naMask <- is.na(modes)

  if (any(naMask))
  {
    correspondingNA <- NA_character_
    modes <- modes[is.na(modes)]
  }
  else
    correspondingNA <- NULL

  correspondingString <-
    FUSIONCONSTRAINTMODES[match_str(modes, FUSIONCONSTRAINTMODES)]

  if (anyNA(correspondingString))
    correspondingString[is.na(correspondingString)] <- NaN

  if (simplify)
    correspondingString <- unique(correspondingString)

 c(correspondingNA, correspondingString)
}

#' @importFrom checkmate assertVector
assertFusionConstraintModes <- function(fusionConstraintModes)
{
  assertVector(fusionConstraintModes, min.len = 1L,
               any.missing = TRUE, all.missing = TRUE)

  if (any(is.nan(corresponding_fusion_mode(fusionConstraintModes,
                                           simplify = FALSE))))
  {
    stop("Some fusion constraint mode have not been recognized")
  }
}

#' @importFrom tidyr expand_grid
#' @importFrom dplyr distinct
fusion_constraint_mode_grid <- function(constraints, modes, pb = NULL,
                                        partition = NULL)
{
  assertFusionConstraints(constraints)
  constraints <- simplify_fusion_constraints(constraints, pb, partition)

  assertFusionConstraintModes(modes)
  modes <- corresponding_fusion_mode(modes)
  modes <- unique(modes)

  grid <- expand_grid(fusionConstraint = constraints,
                      fusionConstraintMode = modes)

  if (any(constraints %in% c(NA_character_, "sizePenalty")))
    noModeMask <- grid$fusionConstraint %in% c(NA_character_, "sizePenalty")
  else
    noModeMask <- logical(nrow(grid))

  if (any(noModeMask))
  {
    grid[noModeMask, "fusionConstraintMode"] <- NA_character_
  }

  if (!all(noModeMask) && anyNA(modes))
  {
    grid[!noModeMask && is.na(grid$modes), "fusionConstraint"] <-
      NA_character_
  }

  distinct(grid)
}
