source("tools/profiling/setup_profiling.R") # nolint

prefixeSauvegarde <- "fusionTree"


profvis(AHR_fusion_tree(pb7x7))
