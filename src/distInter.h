#ifndef DISTINTER_H
#define DISTINTER_H

unsigned int posElementSymVMat(unsigned int i,
                               unsigned int j,
                               const int& dim,
                               const bool& aDefautDiag);

float compSMSymVMat(const arma::vec& values,
                    const unsigned int& dim,
                    const arma::uvec& lignes,
                    const arma::uvec& colonnes,
                    const Rcpp::String& comp,
                    const bool& aDefautDiag,
                    const double& defaultDiag);

float compSMSymMMat(const arma::mat& values,
                    const arma::uvec& lignes,
                    const arma::uvec& colonnes,
                    const Rcpp::String& comp);

Rcpp::NumericVector distanceInterSymVMat(const arma::vec& values,
                                         const arma::uvec& partition,
                                         const Rcpp::NumericMatrix& indexs,
                                         const unsigned int& dim,
                                         const bool& aDefautDiag = true,
                                         const double& defaultDiag = 0.0,
                                         const Rcpp::String& comp = "min");

Rcpp::NumericVector distanceInterSymMMat(const arma::mat& values,
                                         const arma::uvec& partition,
                                         const Rcpp::NumericMatrix& indexs,
                                         const Rcpp::String& comp = "min");

#endif
