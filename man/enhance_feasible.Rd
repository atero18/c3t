% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhancer.R
\name{enhance_feasible}
\alias{enhance_feasible}
\title{Improvement of an existing feasible solution}
\usage{
enhance_feasible(
  regionalisation,
  distances = NULL,
  contiguity = NULL,
  sizes = NULL,
  d = NULL,
  data = NULL,
  m = 0,
  M = Inf,
  standardQuant = FALSE,
  binarQual = FALSE,
  enhanceCriteria = c("AHC", "Silhouette", "Dunn"),
  linkages = "saut max",
  evaluationCriteria = enhanceCriteria,
  maxIt = Inf,
  parallel = TRUE,
  nbCores = detectCores() - 1L,
  verbose = TRUE
)
}
\arguments{
\item{regionalisation}{feasible regionalisation to optimize.}

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

\item{enhanceCriteria}{A vector of criteria used for the enhancement of the
actual feasible solution.
Currently available choices are those in \code{\link[=available_criteria]{available_criteria()}},
plus "AHC" (depends of the \code{linkages} parameter). Compared to others AHC
doesn't improve a global criterion but do this
locally, hoping to reduce computing time. Regarding to this criterion
a feasible solution, built by move a unique element from a cluster to another
is better if the element is closer to the other cluster than it's actual
(depending of some linkage).}

\item{linkages}{Vector of linkage distances used when a criterion ("\code{Dunn}",
"\code{AHC}") needs it.}

\item{evaluationCriteria}{criteria used for comparison after enhancement.
They are evaluated on each feasible solution given by each criterion used
for enhancement. Must be a vector composed of the available criteria in c3t.
For the Dunn index there will be one criterion per linkage given.
See \code{\link[=available_criteria]{available_criteria()}}.}

\item{maxIt}{maximum number of allowed iterations. Default is \code{Inf}.
(strictly positive integer)}

\item{parallel}{Logical indicating whether to use parallel processing.
Default is TRUE.}

\item{nbCores}{Number of CPU cores to use for parallel processing
(sockets method). Default is one less than the detected number of cores.}

\item{verbose}{Logical indicating whether to display progress messages.
Default is TRUE.}
}
\value{
a \link[tibble:tibble-package]{tibble} with one row per try.
For each row the following variables:
\itemize{
\item \code{criterion}: name of the criterion used for improvement.
\item \code{linkage}: type of linkage distance used (\code{NA} if this argument is
irrelevant for the actual criterion).
\item \code{sampleSize}: size of the sample for the calculation of the criterion
(\code{NA} if irrelevant).
\item \code{statut}: state of improvement. Indicates whether an improvement
could be made or not.
\item \code{iterations}: number of improving iterations performed.
\item \code{regionalisationOpti}: the new regionalisation. Identical to the input
argument if no improvement could be made.
\item one column per criterion indicated in  \code{critereEvaluation}. If some of
those criteria use a linkage distance, there will be one column per linkage
distance given in \code{linkage} and per criterion.
}
}
\description{
Algorithm to improve (according to a certain criterion) a
solution that is feasible for a certain classification problem with
connectivity and size constraints. \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
}
\references{
Marc Christine and Michel Isnard. "Un algorithme de regroupement
d'unités statistiques selon certains critères de similitudes"
Insee Méthodes, 2000, p. 50`
}
