#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "Parameters.h"
#include "roboBayes.h"
#include "Model.h"
#include "Unit.h"

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* 
  This function is used when the hyperparameters and sufficient statistics
  were provided by the user in initial step or when doing a continuation
  step
  @param datapts A matrix object {n x d}. The observations for n time points.
  @param covariates A matrix object {n x k}. The covariates for n time points.
  @param RoboBayes A List object. Results of a previous analysis, or the
    initialization steps
  @param paramsList A List object containing analysis settings.

  Function returns a List as defined in function roboBayes(...).
*/

//
// [[Rcpp::export]]
List roboBayesPrior(const arma::mat& datapts,
                    const arma::mat& covariates,
                    const Rcpp::List& RoboBayes,
                    const Rcpp::List& paramsList) {

  Parameters parms(paramsList);

  Unit prevUnit;

  // copy initial model hyperparameters from previous time step
  // and convert to internal class structure
  List p0List = as<List>(RoboBayes["model0"]);
  Model model0(p0List);

  // copy model hyperparameters and sufficient statistics from previous time
  // step and convert into internal class structure
  List pTList = as<List>(RoboBayes["pars"]);
  Model parsT(pTList);
  prevUnit.pars = parsT;

  // joint distribution
  prevUnit.jtR = as<rowvec>(RoboBayes["jtR"]);

  // posterior distribution
  prevUnit.R = as<rowvec>(RoboBayes["R"]);

  // retained run lengths
  prevUnit.RL = as<urowvec>(RoboBayes["RL"]);

  if (parms.getR) {
    prevUnit.RFull = as<mat>(RoboBayes["RFull"]);
  } else {
    prevUnit.RFull = zeros<mat>(0,0);
  }

  mat x_store, y_store; 

  if (parms.getOutliers || parms.getModels) {
    // user requested matrix of posterior distribution for all time points
    prevUnit.Rm = as<mat>(RoboBayes["Rm"]);
    y_store = as<mat>(RoboBayes["y_store"]);
    x_store = as<mat>(RoboBayes["x_store"]);

  } else {
    // if not keeping matrix of posterior distribution, create place holders
    // for these data units
    prevUnit.Rm = zeros<mat>(0,0);
    y_store.zeros(0, parms.d);
    x_store.zeros(0, parms.k);
  }

  if (parms.getOutliers) {
    // time points identified as outliers
    prevUnit.outliers = as<uvec>(RoboBayes["outliers"]);
  } else {
    prevUnit.outliers = zeros<uvec>(0);
  }

  // indices of runs retained
  prevUnit.truncRind = as<uvec>(RoboBayes["truncRind"]);

  // change points identified
  prevUnit.cpInds = as<umat>(RoboBayes["cpInds"]);

  // probability that a change occurred in the previous Lsearch time points
  prevUnit.lastLs = as<vec>(RoboBayes["lastLs"]);

  if (parms.getModels) {
    // convert models for change points to required structure
    List mods = RoboBayes["mods"];
    prevUnit.modelCP = as<uvec>(mods[0]);
    prevUnit.modelB = as<cube>(mods["B"]);
    prevUnit.modelV = as<cube>(mods["V"]);
  } else {
    // placeholders for models for change points
    prevUnit.modelCP = zeros<uvec>(0);
    prevUnit.modelB = zeros<cube>(parms.k, parms.d, 0);
    prevUnit.modelV = zeros<cube>(parms.d, parms.d, 0);
  }

  // run length dependent covariates for retained run lengths
  mat allcov = as<mat>(RoboBayes["allcov"]);

  // initial time is 1
  uword time = as<uword>(RoboBayes["time"]);

  // run analysis
  List result = roboBayes(datapts, covariates, prevUnit, model0, 
                          y_store, x_store, allcov, time, parms);

  return result;

}
