#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "Parameters.h"
#include "Unit.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Function to identify change points

  @param prevUnit The previous time point analysis.
  @param newUnit The current time point analysis.
  @param itime The current time point.
  @param parms The analysis settings

  returns newUnit object with updated cpInds and lastLs
*/

void updateCPs(const Unit& prevUnit,
               Unit& newUnit,
               const int itime,  
               const Parameters& parms) {

  int i1;
  uword i2;

  umat cp_ind = prevUnit.cpInds;

  // posterior current distribution ordered from largest run length to current
  rowvec cpprobs = reverse(newUnit.R);

  // number of retained run lengths
  int nRuns = newUnit.R.n_elem;

  // starting point of search window
  int startingPoint = std::max(nRuns - parms.Lwindow - parms.cp_delay, 0);

  // last point for first group in search window
  int endPoint = startingPoint + parms.Lgroup;

  // sum of the probabilities over all groups within window
  vec cpprobs_L(nRuns - parms.cp_delay - endPoint + 1);
  for (uword i = 0; i < cpprobs_L.n_elem; ++i) {
    i1 = endPoint + i;
    i2 = std::max(0, i1 - parms.Lgroup);
    cpprobs_L[i] = sum(cpprobs.subvec(i2,i1-1));
  }

  uvec bigger, cpLocs, cpWindows, itrunc, newCol(parms.cpthresh.size());
  vec cpProbs;
  rowvec tt;
  uword ij, imaxProb, maxLoc, maxProb;

  for (uword i = 0; i < parms.cpthresh.size(); ++i) {

    uword res;

    if (any(cpprobs_L >= parms.cpthresh[i])) {
      // if any of the grouped probabilities match or exceed the input
      // threshold under consideration:

      // identify which of the grouped probabilities match or exceed
      bigger = find(cpprobs_L >= parms.cpthresh[i]);

      // identify the starting point of the groups that match or exceed
      cpWindows = startingPoint + bigger + 1;

      cpLocs.zeros(cpWindows.n_elem);
      cpProbs.zeros(cpWindows.n_elem);

      for(uword j = 0; j < cpWindows.n_elem; ++j) {
        // for each group identify location and value of max probability
        ij = cpWindows[j];
        tt = cpprobs(span(ij-1,ij-1+parms.Lgroup-1));
        imaxProb = index_max(tt);
        cpLocs[j] = ij + imaxProb;
        cpProbs[j] = tt[imaxProb];
      }

      // keep the max of the group maxes
      maxProb = index_max(cpProbs);
      maxLoc = cpLocs[maxProb];
      itrunc = reverse(newUnit.truncRind);
      res = itime - itrunc(maxLoc-1);

    } else {
      // if no grouped probabilities match or exceed the input
      // threshold under consideration default to 1

      res = 1;

    }

    newCol(i) = res;

  }

  newUnit.cpInds = join_rows(cp_ind, newCol);

  // probability that a change occurred in the previous Lsearch time points

  int imax2 = std::max(nRuns - parms.Lsearch - parms.cp_delay, 1);
  vec lastL(1);
  lastL(0) = sum(cpprobs.subvec(imax2, nRuns-parms.cp_delay - 1));

  newUnit.lastLs = join_cols(prevUnit.lastLs, lastL);

  return;
}
