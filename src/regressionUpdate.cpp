#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "Model.h"
#include "Parameters.h"
#include "Unit.h"
#include "updateHyperParams.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Update sufficient statistics and hyperparameters

  @param datapt observations for a single time step
  @param covariate covariates for a single time step
  @param model0 initial values for hyperparameters
  @param newModel new Model object -- will be modified on return
  @param moreData covariates that depend on run length for all time steps
  @param parms analysis settings

  upon return input newUnit is updated with new hyperparameters
*/

void regressionUpdate(const arma::vec& datapt,
                      const arma::vec& covariate,
                      const Model& model0,
                      Model& newModel,
                      const arma::mat& moreData,
                      const Parameters& parms) {

  // if any observations are not finite, do not modify hyperparameters but
  // insert initial estimates at first slice
  if (!datapt.is_finite()) {
    newModel.join(model0);
    return;
  }

  // number of runs contained in previous time step
  // (okay because newModel has not yet been modified)
  int nRuns = newModel.nu.n_elem;

  mat btmp, md, val0, val1;
  mat zz(covariate.n_elem - moreData.n_rows, moreData.n_cols);

  // pad covariate matrix with 0's for covariates that do not depend
  // on run length and rever order of matrix rows
  md = join_cols(zz, moreData);
  md = reverse(md, 1);

  // calculate the additions to the sufficient statistics for each possible 
  // run length

  // {k x nRuns}
  mat xi = covariate - md.each_col(); 

  vec tcol;
  for (int i = 0; i < nRuns; ++i) {
    tcol = xi.col(i);
    newModel.XTX.slice(i) += tcol * tcol.t();
    newModel.XTY.slice(i) += tcol * datapt.t();
  }

  newModel.YTY.each_slice() += datapt * datapt.t();

  updateHyperParams(model0, newModel, parms);

  newModel.join(model0);

  return;
}
