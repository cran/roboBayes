#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "stepm1.h"
#include "Parameters.h"
#include "spikeUpdate.h"
#include "Unit.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Outlier detection step

  @param cp vector of possible change points
  @param itime current time point
  @param newUnit current time point results
  @param prevUnit previous time point results
  @param model0 initial estimates for hyperparameters
  @param datapts current analysis observations
  @param moreData run length dependent covariates
  @param parms analysis settings
  @param y_store most recent Lm observations.
  @param x_store most recent Lm covariates

  returns an updated newUnit if outlier is detected
*/

void stepCP(arma::uvec& cp, 
            const int itime,
            Unit& newUnit,
            const Unit& prevUnit,
            const Model& model0,
            const arma::mat& datapts,
            const arma::mat& moreData,
            const Parameters& parms,
            const arma::mat& y_store,
            const arma::mat& x_store) {

  field<Unit> sCP(cp.n_elem);

  Unit st1, uniti;

  uword cpi, tcp;

  uvec eliminate, i1, validpts;

  double probY_m0 = sum(newUnit.jtR);
  vec temp1_m0 = newUnit.R.t() * probY_m0 * parms.pc;

  double sumY = sum(temp1_m0), stemp1;

  vec prob_mi(cp.n_elem), rFull(newUnit.RFull.n_rows), temp1;
  
  for (uword i = 0; i < cp.n_elem; ++i) {

    cpi = cp[i];

    // time since change point under consideration
    // note that this is guaranteed to be <  Lm by testing conditions
    // of calling function
    tcp = itime - cpi;

    // consider preceding tcp run lengths
    validpts = regspace<uvec>(tcp+1, 1);

    // if any previously identified outliers are contained in the
    // data for this iteration, remove them
    i1 = find(prevUnit.outliers >= cpi);

    if (i1.n_elem > 0) {

      eliminate = prevUnit.outliers.elem(i1) - cpi;
      eliminate = tcp - eliminate;

      validpts = reverse(validpts);
      validpts.shed_rows(eliminate);
      validpts = reverse(validpts);

    }

    // recover sufficient statistics and update hyperparameters
    // using the data and covariates from the preceding tcp runs
    uniti.pars = spikeUpdate(y_store.rows(validpts-1), 
                             x_store.rows(validpts-1), 
                             tcp, model0, newUnit.pars,
                             moreData, parms);

    // joint probability is taken as the subset of the joint
    // probabilities for the preceding tcp runs
    // note that the length of Rm is increased each iteration; it is
    // not like Rfull, which is predefined to be account for all time steps
    uniti.jtR = prevUnit.Rm(span(tcp, prevUnit.Rm.n_rows-1), tcp).t();

    // create the posterior distribution for the subset
    uniti.R = uniti.jtR / sum(uniti.jtR);
    
    // May 2022
    int maxIndex = std::min(parms.Lsearch + parms.cp_delay - 1, 
                            (int) uniti.R.size() - 2);
    
    vec tlastL(1);
    
    if (parms.cp_delay < maxIndex) {
      tlastL(0) = sum(uniti.R.subvec(parms.cp_delay, maxIndex));
    } else if (parms.cp_delay == maxIndex) {
      tlastL(0) = uniti.R(maxIndex);
    }
    
    if (cpi - parms.cp_delay - 2 >= 0) {
      uniti.lastLs = join_cols(prevUnit.lastLs.subvec(0, cpi-parms.cp_delay-2), tlastL);
    }

    // keep only the run lengths that fall in the preceding tcp runs
    // and shift to account for truncated distribution vectors
    uniti.RL = newUnit.RL.cols(find(newUnit.RL > (tcp+1))) - tcp - 1;

    if (parms.getR) {
      // reset the full posterior distribution
      rFull.zeros();
      rFull.elem(uniti.RL-1) = uniti.R.t();
      uniti.RFull = newUnit.RFull;
      uniti.RFull.col(cpi-1) = rFull;
    }

    if (parms.getOutliers || parms.getModels) {
      uniti.Rm = newUnit.Rm;
    }
    
    // reintroduce each of the preceding tcp runs
    stepm1(tcp, itime, prevUnit, newUnit, uniti, validpts-1, 
           model0, moreData, parms, y_store, x_store);
    
    // store result
    sCP[i] = uniti;

    // calculate probability
    temp1 = uniti.jtR.t() * (1.0 - parms.pc) / ((double) cp.n_elem);

    stemp1 = sum(temp1);
    sumY += stemp1;
    prob_mi[i] = stemp1;
  }

  prob_mi /= sumY;

  if (any(prob_mi > parms.alpha)) {
    // if models with outliers are more likely, remove it and replace 
    // distributions with those without outliers

    uword ind_outlier = index_max(prob_mi);

    Unit toReturn = sCP(ind_outlier);

    rowvec testVector = log(toReturn.R);
    toReturn.R(find_nonfinite(testVector)).fill(DBL_MIN);

    toReturn.RL = join_rows(zeros<urowvec>(1), prevUnit.RL);
    toReturn.RL += 1;

    // store outlier
    uvec cpOut(1, fill::value(cp(ind_outlier)));
    toReturn.outliers = join_cols(prevUnit.outliers, cpOut);

    // remove change point identified
    toReturn.cpInds = newUnit.cpInds;
    urowvec tempCP1 = newUnit.cpInds.tail_rows(1);
    tempCP1.back() = 1;
    toReturn.cpInds.tail_rows(1) = tempCP1;

    toReturn.modelCP = newUnit.modelCP;
    toReturn.modelB = newUnit.modelB;
    toReturn.modelV = newUnit.modelV;
    
    newUnit = toReturn;

    return;

  }

  return;
}
