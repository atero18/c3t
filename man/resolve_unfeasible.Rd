% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resolver.R
\name{resolve_unfeasible}
\alias{resolve_unfeasible}
\title{Resolve an unfeasible regionalisation}
\usage{
resolve_unfeasible(
  distances = NULL,
  contiguity = NULL,
  sizes = NULL,
  d = NULL,
  data = NULL,
  m = 0,
  M = Inf,
  standardQuant = FALSE,
  binarQual = FALSE,
  regionalisation,
  maxItTransfers = Inf,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{distances}{The distance matrix of the problem. This can be omitted
if a distance function \code{d} and data context \code{data} are provided. If only
\code{distances} is provided, all distances must be present. (distance matrix)}

\item{contiguity}{A contiguity matrix or an \code{igraph} contiguity graph. If not
provided, the problem is considered completely contiguous (all elements are
neighbors of each other).}

\item{sizes}{Represents the size of each element. By default, it is set
to \code{1} for each element (the size of a cluster becomes its cardinal).
All data must be positive or zero. (positive real numeric vector)}

\item{d}{Distance function between elements. This can be omitted if
\code{distances} is already indicated. If present, \code{data} must also be
specified. Some classical distances are available, it is recommended to use
them rather than a personal function for optimisation reasons :
\itemize{
\item "\code{euclidean}": Euclidean distance.
\item "\code{manhattan}" : Manhattan distance.
\item "\code{minkowski}" : Minkowski's distance. In that case a value for p >= 1
must be specified.
}

(function or string)}

\item{data}{A data.frame where each row represents data related to an
element. This can be omitted if \code{d} is omitted. Present variables can
be quantitative or qualitative. If qualitative variables are present,
some distances may not be used. Possibility of standardising variables and
transforming qualitative variables into binary variables (one-hot encoding)
using \code{standardQuant} and \code{binarQual}. (data.frame)}

\item{m}{Minimum size constraint. Must be positive or zero and small enough
for the problem to be feasible. Default is \code{0} (no constraint).
(positive number)}

\item{M}{Maximum size constraint. Must be positive, superior or equal to \code{m}
and large enough for the problem to be feasible.
Default is \code{Inf} (no constraint). (positive number)}

\item{standardQuant}{\code{TRUE} if the variables in \code{data} should be
standardised (i.e., centered and scaled), \code{FALSE} (default) otherwise.
Standardisation is applied after the possible binarization of qualitative
variables (see \code{binarQual}). (flag)}

\item{binarQual}{\code{TRUE} if qualitative variables should be binarized (one-hot
encoding), for example, to make the data set compatible with common distances
or to standardize these variables. \code{FALSE} (default) otherwise. (flag)}

\item{regionalisation}{A partition checking the contiguity constraint but
not checking at least the min or the max size constraint.}

\item{maxItTransfers}{positive integer indicating the maximum number
of elements can be transferred.}

\item{verbose}{Logical indicating whether to display progress messages.
Default is TRUE.}

\item{...}{used for development.}
}
\description{
This algorithm tries to transform an unfeasible solution that check
connectivity constraint but not one or more size constraints. Even if
the problem is feasible the algorithm doesn't ensure to succeed.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
}
\examples{
set.seed(123L)
grid  <- simple_grid(4L, 5L)
data <- grid$context
sizes <- grid$repartition$nbIndividuals
contiguity <- grid$contiguity

m <- 200.0
M <- 800.0

d <- "euclidean"

unfeasibleSolution <-
  c(1L, 2L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 4L,
    4L, 2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)

resolve_unfeasible(contiguity = contiguity,
                   sizes = sizes,
                   data = data,
                   d = d, m = m, M = M,
                   regionalisation = unfeasibleSolution,
                   verbose = FALSE)
}
