#ifndef MODEL_H
#define MODEL_H
#include <RcppArmadillo.h>
#include "Parameters.h"

//[[Rcpp::depends(RcppArmadillo)]]

/*
  Class to hold sufficient statistics and hyperparameters
*/


class Model {

  public : 

  arma::cube B;
  arma::cube V;
  arma::cube Lambda;
  arma::cube XTX;
  arma::cube XTY;
  arma::cube YTY;

  arma::vec nu;

  Model(){};

  Model(Parameters& parms);
  Model(Rcpp::List& pars);

  Rcpp::List toList();

  void join(const Model mod);
  void subset(arma::uvec idx);
  void remove(arma::uvec idx);

};


#endif
