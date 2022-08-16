#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "regressionUpdate.h"
#include "regressionProb.h"
#include <mvnorm.h>
#include <limits>
#include "Parameters.h"
#include "Model.h"
#include "Unit.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Move forward in time from change point -- used for outlier detection

  @param tcp time since change point under consideration
  @param itime current time point
  @param prevUnit previous time point results
  @param newUnit current time point results
  @param currUnit results for change point under consideration
  @param validpts vector of runs to include in step
  @param model0 initial estimates for hyperparameters
  @param moreData run length dependent covariates
  @param parms analysis settings
  @param y_store most recent Lm observations.
  @param x_store most recent Lm covariates

*/
void stepm1(const arma::uword tcp, 
            const int itime,
            const Unit& prevUnit,
            const Unit& newUnit,
            Unit& currUnit,
            const arma::uvec& validpts,
            const Model& model0,
            const arma::mat& moreData,
            const Parameters& parms,
            const arma::mat& y_store,
            const arma::mat& x_store) {

  double xmin = std::numeric_limits<double>::min();
  double lmd = log(xmin);

  double mvd = as_scalar(dmvnorm(y_store.row(tcp), 
                                 parms.outlier_mean,
                                 parms.outlier_var));

  double H = 1.0 / parms.lambda;

  uword vpat;
  uvec ids, inVP;

  vec R0(1), rFull(currUnit.RFull.n_rows), temp1, temp2;
  vec xcp(parms.k), ycp(parms.d);
  rowvec growR, templogs, temprow;

  Model tPars;

  // iterate forward from tcp run to current time step
  for (uword ti = (tcp+1); ti > 0; --ti) {

    // if run ti is a valid point (not an outlier), extract
    // the correct data point and covariate. if data is an outlier
    // set data as infinity

    inVP = find(validpts == (ti-1));

    if ((inVP.n_elem > 0) && (ti != (tcp + 1))) {
      ycp = y_store.row(ti-1).t();
      xcp = x_store.row(ti-1).t();
    } else {
      ycp.fill(datum::inf);
      xcp.fill(datum::inf);
    }

    // identify correct time points to pull from run dependent data
    vpat = moreData.n_cols - ti;

    tPars = currUnit.pars;

    // update sufficient statistics and hyperparameteres
    regressionUpdate(ycp, xcp, model0, tPars, 
                     moreData.cols(0, vpat), parms);

    // update joint distribution
    temp1 = regressionProb(ycp, currUnit.pars, tPars, parms);

    templogs = log(currUnit.jtR) + temp1.t();
    templogs.elem(find(templogs < lmd)).fill(lmd);
    templogs.elem(find_nonfinite(templogs)).fill(lmd);

    temprow = exp(templogs) * mvd;

    growR = temprow * (1.0 - H);
    R0(0) = sum(temprow * H);

    currUnit.jtR = join_rows(R0, growR);

    // update posterior distribution
    currUnit.R = currUnit.jtR / sum(currUnit.jtR);
    currUnit.R.elem(find(currUnit.R < xmin)).fill(xmin);

    // update run lengths
    currUnit.RL = join_rows(zeros<uvec>(1), currUnit.RL);
    currUnit.RL += 1;
    
    // May 2022
    int maxIndex = std::min(parms.Lsearch + parms.cp_delay - 1, 
                            (int) currUnit.R.size() - 2);
    
    vec tlastL(1);
    
    if (parms.cp_delay <= maxIndex) {
      tlastL(0) = sum(currUnit.R.subvec(parms.cp_delay, maxIndex));
    }
    
    currUnit.lastLs = join_cols(currUnit.lastLs, tlastL);

    if (parms.getR) {
      // update full matrix of posterior distributions
      rFull.zeros();
      rFull.elem(currUnit.RL(span(0, currUnit.R.n_elem-1))-1) = currUnit.R.t();

      currUnit.RFull.col(itime-ti+1) = rFull;
    }

    if (parms.getOutliers || parms.getModels) {
      // update joint probabilities
      currUnit.Rm(span(ti-1, currUnit.jtR.n_elem+ti-2),ti-1) = currUnit.jtR.t();
    }

    currUnit.pars = tPars;

  }

  return;

}
