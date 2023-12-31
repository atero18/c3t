% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/articulation_points.R
\name{articulation_points}
\alias{articulation_points}
\alias{articulation_pts_cluster}
\alias{articulation_pts_clusters}
\title{Find articulation points}
\usage{
articulation_pts_cluster(contiguity, partition, cluster)

articulation_pts_clusters(
  contiguity,
  partition,
  clusters = clustersIDs(partition)
)
}
\arguments{
\item{contiguity}{A contiguity matrix or a contiguity graph
(package \code{igraph}). (contiguity matrix or \code{igraph} graph)}

\item{partition}{A partition. For \code{is_articulation_pt_cluster},
must be a regionalisation.
(vector of strictly positive integers)}

\item{cluster}{The cluster identifier, not empty for \code{partition}.
(strictly positive integer)}

\item{clusters}{List of non-empty clusters for \code{partition} for which to
determine their articulation point(s). Default is all clusters.}
}
\value{
A vector composed of indices of elements that are articulation
points. (vector of strictly positive integers)

For \code{articulation_pts_clusters} a vector composed of indices of
elements from clusters
in \code{clusters} that are articulation points.
(vector of strictly positive integers)
}
\description{
find articulation points of a partition.
}
\section{Functions}{
\itemize{
\item \code{articulation_pts_cluster()}: Returns the list of articulation points
in a cluster resulting from a partition.

\item \code{articulation_pts_clusters()}: Returns the list of points that are
articulation points for different clusters.

}}
\seealso{
\code{\link[igraph:articulation_points]{igraph::articulation_points()}}
}
\keyword{internal}
