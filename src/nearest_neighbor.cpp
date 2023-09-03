#include <Rcpp.h>

#include "nearest_neighbor.h"

//' @title Find the nearest neighbor of some points
//' @description For a set of points (`subsetPoints`), find the nearest
//' neighbor in another set (`subsetNeighbors`). Distances are given in a
//' `p' x n'` matrix (`distances`), with `p' >= p` (p the number of points)
//' and `n' >= n` (n the number of neighbors). If points and neighbors come from
//' the same set, put inner to `TRUE` for assuring the nearest neighbor of a
//' point will not be itself.
//' @param distances distances between points and neighbors. A row is for a
//' point and a column a neighbor.
//' @param subsetPoints A vector of strictly positive integers giving the
//' position of the points the nearest neighbor must be found. Positions
//' starts from 1 like in R. Must be at least of length one and it would be
//' better if values are unique.
//' @param subsetNeighbors A vector of strictly positive integers giving the
//' position of the elements of the set that can be consider has a neighbor
//' for each point of `subsetPoints`. Positions starts from 1 like in R.
//' Must be at least of length one and it would be
//' better if values are unique.
//' @param contiguity In the case there are contiguity constraints
//' between points. If an element of `subsetPoints` and one
//' of `subsetNeighbors` are not contiguous then the second element cannot
//' be a neighbor of the first one. If no contiguity constraints must be
//' considered a matrix with 0 rows and 0 columns must be sent. Otherwise
//' must be a logical matrix with same number of rows and columns than
//' `distante`.
//' @returns a vector of length `p` with for each point the position of its
//' nearest neighbor (0 if it doesn't exist). Positions starts from 1 like in R.
//' @keywords internal
//' @name cpp_nearest_neighbor


// Neutral contiguity matrix for nearest neighbor research
const Rcpp::LogicalMatrix noContiguityConstraints(0, 0);

//' @describeIn cpp_nearest_neighbor case when `distances` is a matrix
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector nearest_neighbor_matrix(const Rcpp::NumericMatrix& distances,
                                            const Rcpp::NumericVector& subsetPoints,
                                            const Rcpp::NumericVector& subsetNeighbors,
                                            const Rcpp::LogicalMatrix& contiguity,
                                            const bool& inner)
{
  const unsigned int nbPoints = subsetPoints.size();
  Rcpp::NumericVector nns(nbPoints);

  // Check if there is a contiguity constraint or not
  bool contiguityConstraint = contiguity.size() > 0;

  unsigned int point = 0;
  double minDist = 0.0;
  unsigned int noNeighbor = 0;
  unsigned int nn = noNeighbor;
  for(unsigned int k = 0; k < nbPoints; ++k)
  {
    point = subsetPoints(k) - 1;
    minDist = R_PosInf;
    for(auto j: subsetNeighbors)
    {
      --j;

      // If `distances` is a matrix distance for a set and
      // elements in `subsetPoints` are part of this set,
      // we might don't want to consider the diagonal which
      // will have necessaraly the minimum distance because it will
      // be between a point and itself
      if(inner && j == point)
        continue;

      if(!Rcpp::NumericVector::is_na(distances(point, j)) &&
         (!contiguityConstraint || contiguity(point, j)) &&
         (distances(point, j) < minDist || nn == noNeighbor))
      {
        nn = j + 1;
        minDist = distances(point, j);
      }
    }

    nns(k) = nn;
  }

  return nns;
}
