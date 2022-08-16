#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "Model.h"
#include "Parameters.h"
#include "regressionUpdate.h"
#include "regressionProb.h"
#include "stepUK.h"
#include "spikeUpdate.h"
#include "Unit.h"
#include "updateCPs.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Primary analysis function

  @param datapts A matrix object {n x d}. The observations for n time points.
  @param covariates A matrix object {n x k}. The covariates for n time points.
  @param prevUnit A Unit object. The results of initialization or the previous
    time point analysis
  @param model0 A Model object. The initial hyperparameter estimates.
  @param y_store The most recent Lm data points.
  @param x_store The most recent Lm covariates.
  @param time An unsigned integer object. The last time point.
  @param Parameters The analysis settings.

  returns a List object containing
    R A numeric vector object. The current posterior distribution of 
      the run length.
    RL An integer vector object. The current run lengths that are 
       retained.
    truncRind An integer vector object. The indices of the run 
              lengths retained.
    jtR A numeric vector object. The current joint distribution of the 
        run length and data.
    pars A list object. Elements are run length specific summary 
         statistics and hyperparameters. 
    cpInds An integer matrix object {length(cpthresh) x nRuns}. Each 
           element contains the most recent change point.
    lastLs A numeric vector object. Each element contains the 
           probability that a change occurred in the previous Lsearch 
           time points, delayed by cp_delay points.
    outliers An integer vector object. Time points that have been 
             identified as outliers.}
    time An integer object. The current time.
    allcov A numeric matrix object. The run length dependent covariates
          for the retained run lengths. 
    model0 A list object. The initial hyperparameters.
    RFull A numeric matrix object. The posterior distribution for each time point
    Rm The joint distribution.
    modelCP An integer vector object. The identified change points.
    modelB A cube object. The coefficients of the identified change point models.
    modelV A cube object. The variance of the identified change point models
    y_store The most recent Lm data points.
    x_store The most recent Lm covariates.
*/

List roboBayes(const arma::mat& datapts,
               const arma::mat& covariates,
               Unit& prevUnit,
               Model& model0,
               arma::mat& y_store,
               arma::mat& x_store,
               arma::mat& allcov,
               arma::uword time,
               const Parameters& parms) {

  // define variables

  int iData;

  double logEPS = log(DBL_MIN);
  double H = 1.0 / parms.lambda;

  bool oflag;

  // strictly unsigned integers

  uword nNewDays = datapts.n_rows;
  uword nPriorDays = time;
  uword tcp, tempCP;

  uvec eliminate, i1, validpts;

  urowvec tempCP1;

  // doubles

  vec predProb, rFull(prevUnit.RFull.n_rows), tempv1, tempv2;

  rowvec growR, R0(1), temp1, temp2, tR;

  mat moreData, runLengthCov;

  // customs

  Unit newUnit;

  Model update;

  // if covariates depend on run length, extract these from the
  // covariates matrix.
  if (parms.kt > 0) {
    runLengthCov = covariates.cols(parms.k-parms.kt, parms.k-1);
  } else {
    runLengthCov = zeros(0,1);
  }

  for (int it = std::max((int) nPriorDays, 1); 
           it <= ((int) nNewDays + (int) nPriorDays - 1); ++it) {

    // initialize current unit to the values of the previous unit
    newUnit = prevUnit;

    // increment time step
    time += 1;

    // identify row of data corresponding to the current time step
    iData = it - nPriorDays;

    // gaussian related updates
    if (parms.kt > 0) {
      moreData = join_rows(allcov, runLengthCov.row(iData).t());
    } else {
      moreData = join_rows(allcov, runLengthCov);
    }

    // update hyperparameters and sufficient statistics
    // function returns with an updated newUnit.pars
    regressionUpdate(datapts.row(iData).t(), 
                     covariates.row(iData).t(), 
                     model0, 
                     newUnit.pars,
                     moreData, 
                     parms); 

    // update joint distribution

    predProb = regressionProb(datapts.row(iData).t(), 
                              prevUnit.pars,
                              newUnit.pars, 
                              parms);

    temp1 = log(prevUnit.R) + predProb.t();
    temp1.elem(find(temp1 < logEPS)).fill(logEPS);

    temp2 = exp(temp1);

    temp1 = temp2 * (1.0 - H);
    
    R0(0) = sum(temp2) * H;

    newUnit.jtR = join_rows(R0, temp1);

    // update posterior distribution
    newUnit.R = newUnit.jtR / sum(newUnit.jtR);
    newUnit.R.elem(find_nonfinite(log(newUnit.R))).fill(DBL_MIN);

    // update run lengths retained
    newUnit.RL = join_rows(zeros<urowvec>(1), prevUnit.RL);
    newUnit.RL += 1;

    if (parms.getR) {
      // if requested update full posterior distribution matrix
      rFull.zeros();
      rFull.elem(newUnit.RL-1) = newUnit.R.t();
      newUnit.RFull.col(it) = rFull;
    }

    // update Rm matrix
    if (parms.getOutliers || parms.getModels) {
      // shift columns up 1, first column is 0s
      newUnit.Rm.shed_col(parms.Lm-1);
      newUnit.Rm = join_cols(zeros<rowvec>(newUnit.Rm.n_cols), newUnit.Rm);
      newUnit.Rm = join_rows(newUnit.jtR.t(), newUnit.Rm);
     
      // shift data and covariate up 1, first row is for current time point
      y_store.shed_row(parms.Lm-1);
      y_store = join_cols(datapts.row(iData), y_store);

      x_store.shed_row(parms.Lm-1);
      x_store = join_cols(covariates.row(iData), x_store);
    }

    oflag = false;

    // update the change point list if the minimum time is exceeded
    if ((int)(time-1) >= parms.cptimemin) {

      // identify change points and calculated probability of change
      // occurring in previous Lsearch time points
      updateCPs(prevUnit, newUnit, time, parms);

      if (parms.getOutliers) {
        // if probabilities used for determining outlier is requested by
        // user perform outlier detection step

        // retrieve the last change point of the largest threshold
        tempCP1 = newUnit.cpInds.tail_rows(1);
        tempCP = tempCP1.back();
        if (tempCP != 1) { 
          oflag = true;
        }

        // if an outlier was detected, newUnit will be updated
        stepUK(tempCP, newUnit, prevUnit, model0, it, parms,
               datapts, moreData, y_store, x_store);

      }

    } else {

      // if the minimum time is not yet exceeded, pad with defaults

      newUnit.cpInds.ones(parms.cpthresh.size(), 1);
      newUnit.lastLs = join_cols(prevUnit.lastLs, zeros<vec>(1));

    }

    // collect most likely models every time there is a change point
    if (parms.getModels) {

      // most recent change point for largest threshold
      tempCP1 = newUnit.cpInds.tail_rows(1);
      tempCP = tempCP1.back();

      if (oflag && tempCP == 1) {
        i1 = unique(prevUnit.cpInds);
        tempCP = i1.back();
      }

      // if change point is not 1 and time since last change point is
      // less than Lm
      if ((tempCP != 1) && ((it-tempCP) < (uword) parms.Lm)) {

        // time since change point
        tcp = it - tempCP;

        // initially, include all possible time points in window since
        // last change point
        validpts = regspace<uvec>(tcp+1, 1);

        // remove point in window that correspond to outliers
        i1 = find(newUnit.outliers >= tempCP);

        if (i1.n_elem > 0) {
          // remove outlier
          eliminate = newUnit.outliers.elem(i1) - tempCP;
          eliminate = tcp - eliminate;

          validpts = reverse(validpts);
          validpts.shed_rows(eliminate);
          validpts = reverse(validpts);
        }

        // recover sufficient statistics and hyperparameters

        update = spikeUpdate(y_store.rows(validpts-1), 
                             x_store.rows(validpts-1),
                             tcp, model0, newUnit.pars,
                             moreData, parms);

        // subset joint probability to the window
        vec jtRtemp = prevUnit.Rm(span(tcp, prevUnit.Rm.n_rows-1), tcp).as_col();

        // calculate posterior distribution for this subset
        vec Rtemp = jtRtemp / sum(jtRtemp);

        // identify max posterior
        uword imax = index_max(Rtemp);

        // store model according to change point
        uvec modInd = find(tempCP == newUnit.modelCP);

        if (modInd.n_elem > 0) {
          newUnit.modelCP(modInd[0]) = tempCP;
          newUnit.modelB.slice(modInd[0]) = update.B.slice(imax);
          newUnit.modelV.slice(modInd[0]) = update.V.slice(imax) / 
            (update.nu(imax)-parms.d-1);
        } else {
          uvec ttempCP(1, fill::value(tempCP));
          newUnit.modelCP = join_cols(newUnit.modelCP, ttempCP);

          cube tb(parms.k, parms.d, 1);
          tb.slice(0) = update.B.slice(imax);
          newUnit.modelB = join_slices(newUnit.modelB, tb);

          cube tv(parms.d, parms.d, 1);
          tv.slice(0) = update.V.slice(imax) / 
                         (update.nu(imax)-parms.d-1);
          newUnit.modelV = join_slices(newUnit.modelV, tv);
        }
      }
    }

    // truncate R and pars if necessary
    // check if the R vectors need to be truncated based on threshold
    if ((int) newUnit.R.n_elem > parms.truncRmin) {

      // remove the first truncRmin runs
      tR = newUnit.R;
      tR.shed_cols(regspace<uvec>(0, parms.truncRmin-1));

      // identify any elements of the "extra" that are at or below the threshold
      uvec id1 = find(tR <= parms.truncRthresh);

      if (id1.n_elem > 0) {
         // initialize the retained indices to be all elements of R
        newUnit.truncRind = regspace<uvec>(1, newUnit.R.n_elem);

        // remove the elements that correspond to the elements
        // not satisfying the threshold
        newUnit.truncRind.shed_rows(id1 + parms.truncRmin);

        if (parms.getOutliers || parms.getModels) {

          // normalize each column of the probabilities
          mat s1 = newUnit.Rm;
          s1.each_row() /= sum(s1);

          // initialize temp indices to be all cases and
          // remove those for which the total probability
          // of the last Lm steps is below the threshold
          uvec truncRind_m = regspace<uvec>(1, newUnit.R.n_elem);

          for (int i = s1.n_rows-1; i >= 0; --i) {
            if (all(s1.row(i) <= parms.truncRthresh)) {
              truncRind_m.shed_row(i);
            }
          }

          // incorporate any kept indices into the retained index vector
          newUnit.truncRind = sort(unique(join_cols(newUnit.truncRind, truncRind_m)));
        }

        // truncate all results to include only those identified above
        // DO NOT truncate RFull

        newUnit.R = newUnit.R.cols(newUnit.truncRind-1);

        newUnit.jtR = newUnit.jtR.cols(newUnit.truncRind-1);
        newUnit.RL = newUnit.RL.cols(newUnit.truncRind-1);

        if (parms.getOutliers || parms.getModels) {
          newUnit.Rm = newUnit.Rm.rows(newUnit.truncRind-1);
        }

        newUnit.pars.subset(newUnit.truncRind-1);

        moreData = reverse(moreData,1);
        uvec ttrunc = newUnit.truncRind;
        ttrunc.shed_row(0);
        moreData = moreData.cols(ttrunc - 2);
        moreData = reverse(moreData, 1);

      } else {
        newUnit.truncRind = regspace<uvec>(1, newUnit.R.n_elem);
      }
    } else {
      newUnit.truncRind = regspace<uvec>(1, newUnit.R.n_elem);
    }

    prevUnit = newUnit;

    allcov = moreData;
  }


  List outPars = prevUnit.pars.toList();
  List outPars0 = model0.toList();

  List toReturn = List::create(Named("RFull", prevUnit.RFull), 
                               Named("R", prevUnit.R), 
                               Named("RL", prevUnit.RL), 
                               Named("truncRind", prevUnit.truncRind),
                               Named("Rm", prevUnit.Rm),
                               Named("jtR", prevUnit.jtR),
                               Named("pars", outPars), 
                               Named("cpInds", prevUnit.cpInds),
                               Named("lastLs", prevUnit.lastLs),
                               Named("outliers", prevUnit.outliers),
                               Named("modelCP", prevUnit.modelCP),
                               Named("modelB", prevUnit.modelB),
                               Named("modelV", prevUnit.modelV),
                               Named("y_store", y_store),
                               Named("x_store", x_store),
                               Named("time", time), 
                               Named("allcov", moreData),
                               Named("model0", outPars0));

  return toReturn;
}




