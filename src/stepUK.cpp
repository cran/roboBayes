#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "Parameters.h"
#include "stepCP.h"
#include "Unit.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Outlier detection step

  @param cp change point
  @param newUnit current time point results
  @param prevUnit previous time point results
  @param model0 initial estimates for hyperparameters
  @param itime current time point
  @param parms analysis settings
  @param datapts current analysis observations
  @param moreData run length dependent covariates
  @param y_store most recent Lm observations.
  @param x_store most recent Lm covariates

  returns a newUnit updates with...
*/

void stepUK(const arma::uword cp, 
            Unit& newUnit,
            const Unit& prevUnit,
            const Model& model0,
            const int itime,
            const Parameters& parms,
            const arma::mat& datapts,
            const arma::mat& moreData,
            const arma::mat& y_store,
            const arma::mat& x_store) {

  bool doit = false;

  bool inOutliers = any(prevUnit.outliers == cp);

  if ((cp != 1) && 
      ((int) cp >= (itime - parms.Lm + 2)) && 
      !inOutliers) {
    doit = true;
  }

  if (!doit) {
    // if change point is 1, previously identified, or more than Lm time 
    // points back, do nothing.
    return;
  }

  uvec cpt;

  // consider a window of change points surrounding identified change point
  int mc = cp-2;
  int minc = std::max(1, mc);
  mc = cp+2;
  cpt = regspace<uvec>(minc, mc);

  // this window cannot contain the first run
  uvec tst = find(cpt == 1);
  cpt.shed_rows(tst);

  // it cannot exceed the current time point
  tst = find(cpt > (itime-1));
  cpt.shed_rows(tst);

  // and it cannot be more than Lm points back
  tst = find(cpt < (itime-parms.Lm+2));
  cpt.shed_rows(tst);

  if (cpt.n_elem == 0) {
    // if no change points in this window satisfy the conditions,
    // return without changes
    return;
  }

  // remove any change points that overlap with identified outliers
  uvec inters = intersect(prevUnit.outliers, cpt);
  if (inters.n_elem > 0) {
    for (uword i = 0; i < inters.n_elem; ++i) {
      tst = find(cpt == inters[i]);
      cpt.shed_rows(tst);
    }
  }

  if (cpt.n_elem == 0) {
    // if no change points in this window satisfy the conditions,
    // return without changes
    return;
  }

  // check for outlier
  stepCP(cpt, itime, newUnit, prevUnit, model0,
         datapts, moreData, parms,
         y_store, x_store);

  return;

}
