#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "Model.h"
#include "Parameters.h"
#include "updateHyperParams.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Recover hyperparameters and sufficient statistics

  @param datapts subset of observations for current subset of Lm
  @param covariate subset of covariates for current subset of Lm
  @param tcp number of terms in subset (this may not be the same as
         the size of datapts due to previously identified outliers)
  @param model0 initial estimates for hyperparameters sand sufficient statistics
  @param modelNew current hyperparameters and sufficient statistics
  @param moreData run length dependent covariates
  @param parms analysis settings

  returns a new Model with recovered hyperparameters and sufficient statistics
*/

Model spikeUpdate(const arma::mat& datapt,
                  const arma::mat& covariate,
                  const arma::uword tcp,
                  const Model& model0,
                  const Model& modelNew,
                  const arma::mat& moreData,
                  const Parameters& parms) {

  if (!datapt.is_finite()) {
    // if data is not finite don't make any changes
    return modelNew;
  }

  // number of run lengths at current step
  uword nRuns = modelNew.nu.n_elem;

  // extract the last tcp runs
  Model starModel = modelNew;
  starModel.subset(regspace<uvec>(tcp+1, nRuns-1));

  mat md, zz(covariate.n_cols - moreData.n_rows, moreData.n_cols);

  // pad run dependent covariates to account for non-run dependent terms
  md = join_cols(zz, moreData);
  md = reverse(md, 1);

  // recover the sufficient statistics

  nRuns = starModel.nu.n_elem;

  for (uword i = tcp; i < md.n_cols; ++i) {
    mat xi = covariate.each_row() - md.col(i).t();
    starModel.XTX.slice(i-tcp) -= xi.t() * xi;
    starModel.XTY.slice(i-tcp) -= xi.t() * datapt;
  }

  starModel.YTY.each_slice() -= datapt.t() * datapt;

  updateHyperParams(model0, starModel, parms);
  starModel.nu -= 1.0 + covariate.n_rows;

  return starModel;

}
