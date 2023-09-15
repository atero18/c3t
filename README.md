
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
- [✅ Conclusion](#-conclusion)
- [📚 References](#-references)

## 🪛 Installation

[![Licence](https://img.shields.io/github/license/atero18/c3t?style&color=5D6D7E)](GPL%203)

You can install the development version of `c3t` from
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
- Apply a modified [Agglomerative Hierarchical
  Clustering](https://www.wikiwand.com/en/Hierarchical_clustering#introduction)
  (AHC) algorithm.
- Choose from various linkage methods to define cluster proximity.
- Apply constraints on the merging process to ensure specific criteria
  are met.

### 2. Improving Feasible Solutions

- Enhance feasible solutions while preserving constraints.
- Optimize solutions according to user-defined criteria.
- Explore different linkage methods for clustering.

### 3. Addressing Unfeasible Solutions

- Convert unfeasible solutions into feasible ones.
- Adjust cluster sizes to meet constraints.
- Maintain contiguity and other specified criteria.

### 4. Flexible Parameter Customization

- Choose distance measures, such as Euclidean distance.
- Set minimum and maximum size constraints for regions.
- Specify initial partitions for improved results.
- Apply various criteria for evaluation, including Caliński-Harabasz
  index (CHI) (Caliński and Harabasz 1974).

### 5. Versatility in Constraints

- Customize constraints based on your specific requirements.
- Ensure contiguity and size constraints are met.
- Explore different fusion constraint modes.

## 🚀 Getting Started

Before we begin, ensure you have the `c3t` package and the required
dependencies installed. You can do this by running the following
commands:

``` r
# load the c3t package
library(c3t)
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
nbEmptyZones <- 3L
nbMetropolises <- 2L
nbVariablesQuant <- 2L

grid <- gen_grid(x, y, nbIndividuals = nbIndividuals,
                 nbMinEmptyZones = nbEmptyZones,
                 nbMetropolises = nbMetropolises,
                 nbQuantitatives = nbVariablesQuant)

data <- grid$context
individuals <- grid$repartition$nbIndividuals
contiguity <- grid$contiguity
```

### Hierarchical Regionalization (AHR)

To address regionalization problems with constraints like minimum and
maximum size, you can use the `AHR` function. This function applies a
modified Hierarchical Agglomerative Clustering (HAC) algorithm while
respecting contiguity and size constraints.

``` r
resAHR <- AHR(contiguity = contiguity,
              d = "euclidean", data = data,
              sizes = individuals,
              m = 5.0, M = 40.0,
              criteria = "CHI",
              fusionConstraints = available_fusion_constraints(),
              fusionConstraintModes = available_fusion_modes(),
              parallel = FALSE)
#> ℹ Starting time: 2023-09-15 08:45:50.187542
#> ℹ 45 AHC to evaluate
#> → 362 non-trivial regionalisations obtained
#> ✔ 62 feasible partitions obtained
#> → 19 redundancies have been removed.
#> → Calculation of the CHI criterion
#> → Execution time: 11.103709936142
```

The function returns a list of feasible solutions, and you can select
the one that best suits your needs.

### Improving Feasible Solutions

Once you have a feasible solution, you can further enhance it while
preserving constraints using the `enhance_feasible` function.

``` r
resEnhance <-
  enhance_feasible(regionalisation = resAHR$results$partition[[1L]],
                   contiguity = contiguity,
                   d = "euclidean", data = data,
                   sizes = individuals,
                   m = 5.0, M = 40.0,
                   enhanceCriteria = c("AHC", "CHI"),
                   linkages = c("single", "complete"),
                   parallel = FALSE,
                   verbose = TRUE)
#> → Evaluation of the 3 enhancements
#> → Calculation of 1 evaluation criteria on the initial partition
#> → Calcul of 1 evaluation criteria on the 3 enhanced partitions
```

This function allows you to improve your solution according to specified
criteria.

### Addressing Unfeasible Solutions

In cases where a feasible solution cannot be obtained, the
`resolve_unfeasible` function attempts to provide a feasible solution
from an unfeasible one.

``` r
unfeasibleReg <- c(1L, 2L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 4L,
                   4L, 2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)

resolution <- resolve_unfeasible(contiguity = contiguity,
                                 sizes = individuals,
                                 data = data,
                                 d = "euclidean", m = 5.0, M = 40.0,
                                 regionalisation = unfeasibleReg,
                                 verbose = TRUE)
#> → Transfert of elements one-by-one
#> ✔ fully resolved partition
```

This function aims to transform an unfeasible solution into a feasible
one while respecting constraints.

## ✅ Conclusion

The `c3t` package offers a variety of tools for regionalization and
clustering with constraints. Explore the documentation and experiment
with different functions to suit your specific use case.

For more details on each function and available options, refer to the
package documentation :

``` r
help(package = "c3t")
```

For any issue or feedback please go to the [Issue
page](https://github.com/atero18/c3t/issues) of the repository. You can
as well discuss about the project in the [Discussion
page](https://github.com/atero18/c3t/discussions).

Please note that this project is released with a [Contributor Code of
Conduct](https://atero18.github.io/c3t/CODE_OF_CONDUCT.html). By
participating in this project you agree to abide by its terms.

## 📚 References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-calinski_dendrite_1974" class="csl-entry">

Caliński, T., and J. Harabasz. 1974. “A Dendrite Method for Cluster
Analysis.” *Communications in Statistics - Theory and Methods* 3 (1):
1–27. <https://doi.org/10.1080/03610927408827101>.

</div>

<div id="ref-christine_algorithme_2000" class="csl-entry">

Christine, Marc, and Michel Isnard. 2000. “Un Algorithme de Regroupement
d’unités Statistiques Selon Certains Critères de Similitude.” *Insee
Méthodes*, 50.
<http://jms-insee.fr/2000/S03_4_ACTE_CHRISTINE-ISNARD_JMS2000.PDF>.

</div>

</div>
