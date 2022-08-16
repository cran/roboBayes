#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <wishart.h>
#include "Model.h"
#include "Parameters.h"
#include "regressionUpdate.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Calculate joint probabilities

  @param datapt observations for a single time step
  @param prevModel results from the previous time step
  @param newModel results from the current time step
  @param parms Analysis settings

  @returns a vector of probabilities
*/

arma::vec regressionProb(const arma::vec& datapt,
                         const Model& prevModel,
                         const Model& newModel,
                         const Parameters& parms) {

  // number of runs in the previous time step
  uword nRuns = prevModel.nu.n_elem;

  // if data is not finite, return vector of 0s
  if (!datapt.is_finite()) {
    return zeros(nRuns);
  }

  // {nRuns}
  vec c1(nRuns), c2(nRuns), v2(nRuns), v3(nRuns), v4(nRuns), v5(nRuns);

  for (uword i = 0; i < nRuns; ++i) {
    v2[i] = log(det(prevModel.Lambda.slice(i)));
    v3[i] = log(det(newModel.Lambda.slice(i+1)));
    v4[i] = log(det(newModel.V.slice(i+1)));
    v5[i] = log(det(prevModel.V.slice(i)));
    c1[i] = lmvgamma(parms.d, (prevModel.nu(i) + 1.0) / 2.0);
    c2[i] = lmvgamma(parms.d, prevModel.nu(i) / 2.0);
  }
  vec v1 = -((double) parms.d/2.0) * log(datum::pi) + c1 - c2;

  vec predProbs = v1 + ((double) parms.d / 2.0) * (v2 - v3) -
                  ((prevModel.nu+1.0) / 2.0) % v4 + (prevModel.nu / 2.0) % v5;

  return predProbs;
  
}
