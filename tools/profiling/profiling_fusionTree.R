rm(list = ls())
source("profiling/core_profiling.R")

prefixeSauvegarde = "fusionTree"


profvis(AHR_fusion_tree(pb7x7))
