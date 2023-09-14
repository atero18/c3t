
<!-- README.md is generated from README.Rmd. Please edit that file -->

# c3t : Clustering with Connectivity and Size Constraints

<!-- badges: start -->

[![R-CMD-check](https://github.com/atero18/c3t/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/atero18/c3t/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/atero18/c3t/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/atero18/c3t/actions/workflows/test-coverage.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/atero18/c3t/branch/main/graph/badge.svg)](https://app.codecov.io/gh/atero18/c3t?branch=main)
<!-- badges: end -->

## 📒 Table of Contents

- [📒 Table of Contents](#-table-of-contents)
- [🪛 Installation](#-installation)
- [📍 Overview](#-overview)
- [🚀 Getting Started](#-getting-started)

## 🪛 Installation

[![Licence](https://img.shields.io/github/license/atero18/c3t?style&color=5D6D7E)](GPL%203)

You can install the development version of c3t from
[GitHub](https://github.com/atero18/c3t) with:

``` r
# install.packages("devtools")
devtools::install_github("atero18/c3t")
```

## 📍 Overview

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

The `c3t` package in R is a powerful tool for regionalization and
clustering tasks with various constraints. It offers a range of
functions to facilitate these tasks, and its capabilities include:

### 1. Hierarchical Regionalization (AHR)

- Solve regionalization problems while respecting constraints such as
  minimum and maximum size.
- Apply a modified Hierarchical Agglomerative Clustering (HAC)
  algorithm.
- Choose from various linkage methods to define cluster proximity.
- Apply constraints on the merging process to ensure specific criteria
  are met.

### 2. Improving Feasible Solutions

- Enhance feasible solutions while preserving constraints.
- Optimize solutions according to user-defined criteria.
- Explore different linkage methods for clustering.

### 3. Addressing Unfeasible Solutions

- Convert infeasible solutions into feasible ones.
- Adjust cluster sizes to meet constraints.
- Maintain contiguity and other specified criteria.

### 4. Flexible Parameter Customization

- Choose distance measures, such as Euclidean distance.
- Set minimum and maximum size constraints for regions.
- Specify initial partitions for improved results.
- Apply various criteria for evaluation, including Calinski-Harabasz
  index (CHI).

### 5. Versatility in Constraints

- Customize constraints based on your specific requirements.
- Ensure contiguity and size constraints are met.
- Explore different fusion constraint modes.

### 6. Detailed Documentation

- Comprehensive documentation to guide users through each function.
- Examples and explanations of available options.

The `c3t` package is a versatile and robust solution for addressing
regionalization and clustering problems with constraints, making it a
valuable tool for spatial analysis and geographic research.

For detailed information on function usage and options, please refer to
the package documentation.

## 🚀 Getting Started

Before we begin, ensure you have the `c3t` package and the required
dependencies installed. You can do this by running the following
commands:

``` r
library(c3t)
library(dplyr)
```

### Creating a Grid

Let’s start by creating a fictitious grid, which will serve as an
example for using different functions. We can generate this grid using
the `gen_grid` function. You can specify the number of individuals on
the grid, empty cells, and metropolitan areas.

``` r
set.seed(123L)
x <- 4L
y <- 5L
nbIndividuals <- 100L
nbCasesVides <- 3L
nbMetropolises <- 2L
nbVariablesQuant <- 2L

grille <- gen_grid(x, y, nbIndividuals = nbIndividuals,
                   nbMinEmptyZones = nbCasesVides,
                   nbMetropolises = nbMetropolises,
                   nbQuantitatives = nbVariablesQuant)

data <- grille$context
individus <- grille$repartition$nbIndividuals
contiguite <- grille$contiguity
```

### Hierarchical Regionalization (AHR)

To address regionalization problems with constraints like minimum and
maximum size, you can use the `AHR` function. This function applies a
modified Hierarchical Agglomerative Clustering (HAC) algorithm while
respecting contiguity and size constraints.

``` r
resRAH <- AHR(contiguity = contiguite,
              d = "euclidean", data = data,
              sizes = individus,
              m = 5.0, M = 40.0,
              criteria = "CHI",
              fusionConstraints = available_fusion_constraints(),
              fusionConstraintModes = available_fusion_modes(),
              parallele = FALSE)
#> ℹ Heure de début : 2023-09-14 14:28:27.852455
#> ℹ 45 CAH à évaluer
#> → 362 partitions non triviales obtenues
#> ✔ 62 partitions faisables obtenues
#> → 19 redondances ont été supprimées
#> → Calcul du critère CHI
#> → Temps d exécution : 8.30705404281616
```

The function returns a list of feasible solutions, and you can select
the one that best suits your needs.

### Improving Feasible Solutions

Once you have a feasible solution, you can further enhance it while
preserving constraints using the `enhance_feasible` function.

``` r
resEnhance <- enhance_feasible(regionalisation = resRAH$results$partition[[1L]],
                               contiguity = contiguite,
                               d = "euclidean", data = data,
                               sizes = individus,
                               m = 5.0, M = 40.0,
                               enhanceCriteria = c("AHC", "CHI"),
                               linkages = c("single", "complete"),
                               parallele = FALSE,
                               verbose = TRUE)
#> → Évaluation de 3 améliorations
#> → Calcul de 1 critère d'évaluation sur la partition initiale
#> → Calcul de 1 critère d'évaluation sur les 3 partitions améliorées
```

This function allows you to improve your solution according to specified
criteria.

### Addressing Unfeasible Solutions

In cases where a feasible solution cannot be obtained, the
`resolve_unfeasible` function attempts to provide a feasible solution
from an infeasible one.

``` r
regInfaisable <- c(1L, 2L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 4L,
                   4L, 2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)

resolution <- resolve_unfeasible(contiguity = contiguite,
                                 sizes = individus,
                                 data = data,
                                 d = "euclidean", m = 5.0, M = 40.0,
                                 regionalisation = regInfaisable,
                                 verbose = TRUE)
#> → Transfert des éléments un à un
#> ✔ Partition totalement résolue
```

This function aims to transform an infeasible solution into a feasible
one while respecting constraints.

## Conclusion

The `c3t` package offers a variety of tools for regionalization and
clustering with constraints. Explore the documentation and experiment
with different functions to suit your specific use case.

For more details on each function and available options, refer to the
package documentation :

``` r
help(package = "c3t")
```
