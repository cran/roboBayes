#include "Unit.h"
#include "Model.h"
#include "Parameters.h"

Rcpp::List roboBayes(const arma::mat& datapts,
                     const arma::mat& covariates,
                     Unit& prevUnit,
                     Model& model0,
                     arma::mat& y_store,
                     arma::mat& x_store,
                     arma::mat& allcov,
                     arma::uword time,
                     const Parameters& parms);
