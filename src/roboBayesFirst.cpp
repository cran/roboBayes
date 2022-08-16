#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "roboBayes.h"
#include "Model.h"
#include "Parameters.h"
#include "Unit.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/* 
  This function is used when the hyperparameters and sufficient statistics
  were not provided by the user
  @param datapts A matrix object {n x d}. The observations for n time points.
  @param covariates A matrix object {n x k}. The covariates for n time points.
  @param paramsList A List object containing analysis settings.

  Function returns a List as defined in function roboBayes(...).
*/

//
// [[Rcpp::export]]
List roboBayesFirst(const arma::mat& datapts,
                    const arma::mat& covariates,
                    const Rcpp::List& paramsList) {

  Parameters parms(paramsList);

  Unit prevUnit;

  // model hyperparameters initialization

  // model coefficients
  Model model0;
  model0.B.zeros(parms.k, parms.d, 1);

  // model variance
  model0.V.zeros(parms.d, parms.d, 1);
  model0.V.fill(0.1);
  model0.V.slice(0).diag().ones();

  // scale parameter of the matrix Normal.
  model0.Lambda.zeros(parms.k, parms.k, 1);
  model0.Lambda.slice(0).diag().fill(0.01);

  // degrees of freedom of the Inverse Wishart.
  model0.nu.zeros(1);
  model0.nu = parms.d - 1.0 + 0.1;

  // sufficient statistics initialization

  model0.XTX.zeros(parms.k, parms.k, 1);
  model0.XTY.zeros(parms.k, parms.d, 1);
  model0.YTY.zeros(parms.d, parms.d, 1);

  prevUnit.pars = model0;

  // joint distribution
  prevUnit.jtR = ones<rowvec>(1);

  // posterior distribution
  prevUnit.R = ones<rowvec>(1);

  // retained run lengths
  prevUnit.RL = ones<urowvec>(1);

  if (parms.getR) {
    // if user requested matrix of posterior distribution for all time points
    // initialize matrix to number of time points in data
    prevUnit.RFull = zeros<mat>(datapts.n_rows+1, datapts.n_rows+1);
    prevUnit.RFull(0,0) = 1.0;
  } else {
    // if not keeping matrix of posterior distribution, create place holders
    // for these data units
    prevUnit.RFull = zeros<mat>(0, 0);
  }

  mat y_store, x_store;

  if (parms.getOutliers || parms.getModels) {
    // if user requested matrix of probabilities used for outlier elimination
    // or model parameters for each change point, initialize matrices
    prevUnit.Rm = ones<mat>(1, parms.Lm);
    y_store = zeros(parms.Lm, parms.d);
    x_store = zeros(parms.Lm, parms.k);

  } else {
    // if not keeping matrix of probabilities used for outlier elimination
    // or model parameters create place holders the data units
    prevUnit.Rm = ones<mat>(1, 0);
    y_store = zeros(0, parms.d);
    x_store = zeros(0, parms.k);
  }

  // time points identified as outliers
  prevUnit.outliers = zeros<uvec>(0);

  // indices of runs retained
  prevUnit.truncRind = ones<uvec>(1);

  // change points identified
  prevUnit.cpInds = ones<umat>(parms.cpthresh.n_elem, 0);

  // probability that a change occurred in the previous Lsearch time points
  prevUnit.lastLs = zeros<vec>(0);

  // models for change points -- these act as placeholders if getModels = FALSE
  prevUnit.modelCP = zeros<uvec>(0);
  prevUnit.modelB = zeros<cube>(parms.k, parms.d, 0);
  prevUnit.modelV = zeros<cube>(parms.d, parms.d, 0);

  // run length dependent covariates for retained run lengths
  mat allcov(parms.kt,0);

  // initial time is 1
  uword time = 1;

  // run analysis
  List result = roboBayes(datapts, covariates, prevUnit, model0, 
                          y_store, x_store, allcov, time, parms);

  return result;

}
