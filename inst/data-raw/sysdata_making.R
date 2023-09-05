#' @include distances_lilnkage.R
#' @importFrom tibble tibble
LINKAGES <- # nolint: object_name_linter
  tibble(linkage = c("Single", "Complete", "Centroid",
                     "Average", "Medoid", "Hausdorff"),
         names = list(c("single", "saut min", "saut_min"),
                      c("complete", "saut max", "saut_max",
                        "diametre", "diametres"),
                      c("centroid", "barycentre", "barycentres",
                        "barycentres", "UPGMC"),
                      c("average", "moyenne", "moyennes",
                        "mean", "UPGMA"),
                      c("medoid", "medoide", "medoides"),
                      c("hausdorff", "haus")),
         fun = list(single_linkage, complete_linkage,
                    centroid_linkage,
                    mean, medoid, distance_Hausdorff),
         needsd = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
         needsData = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
         hclust = c("single", "complete", "centroid",
                    "average", NA_character_, NA_character_),
         cppMode = c("min", "max", NA_character_,
                     "mean", NA_character_,
                     "hausdorff"),
         unitary = TRUE)

description <-
  c("single linkage distance (minimum)",
    "complete linkage distance (maximum)",
    paste0("distance between cluster centroids. Note: ",
           "`d` must be specified here, and variables in ",
           "`data` must all be quantitative."),
    "average distance between points within clusters (mean)",
    "distance between cluster medoids",
    "minimum radius for each neighborhood cluster constains the other cluster.")

#' @importFrom tibble add_column
LINKAGES <- add_column(LINKAGES,
                       description = description,
                       .after = "linkage")

rm(description)

#' @include fusion_constraint.R

#' @title List of available constraints fusion types
#' @name fusion_constraints
#' @rdname fusion_constraints
#' @keywords internal
NULL

#' @importFrom tibble tibble
FUSIONCONSTRAINTS <- # nolint: object_name_linter
  tibble(fusionConstraint = c("singletonMin",
                              "singletonAny", "minMin", "sizePenalty"),
         names = list("singletonMin",
                      "singletonAny",
                      "minMin",
                      c("sizePenalty", "penaliteTaille")),
         needsMode = c(TRUE, TRUE, TRUE, FALSE),
         needsMinConstraint = c(TRUE, TRUE, TRUE, FALSE))

#' @describeIn fusion_constraints Authorized modes of fusion constraint
#' when one is activated.
#' @keywords internal
FUSIONCONSTRAINTMODES <- c("deterministic", "random") # nolint: object_name_linter

usethis::use_data(LINKAGES, FUSIONCONSTRAINTS, FUSIONCONSTRAINTMODES,
                  internal = TRUE, overwrite = TRUE)

rm(LINKAGES, FUSIONCONSTRAINTS, FUSIONCONSTRAINTMODES)
