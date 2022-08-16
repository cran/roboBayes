#ifndef UNIT_H
#define UNIT_H
#include <RcppArmadillo.h>
#include "Model.h"

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

class Unit {

  public:

  Model pars;
  arma::rowvec jtR;
  arma::rowvec R;
  arma::urowvec RL;
  arma::mat RFull;
  arma::mat Rm;

  arma::uvec truncRind;

  arma::uvec outliers;

  arma::umat cpInds;
  arma::vec lastLs;

  arma::uvec modelCP;
  arma::cube modelB;
  arma::cube modelV;

};

#endif
