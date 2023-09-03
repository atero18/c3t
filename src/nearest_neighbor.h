#ifndef NEAREST_NEIGHBOR_H
#define NEAREST_NEIGHBOR_H


Rcpp::NumericVector nearest_neighbor_matrix(const Rcpp::NumericMatrix& distances,
                                            const Rcpp::NumericVector& subsetPoints,
                                            const Rcpp::NumericVector& subsetNeighbors,
                                            const Rcpp::LogicalMatrix& contiguity,
                                            const bool& inner);

#endif
