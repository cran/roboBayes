#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "Model.h"
#include "Parameters.h"
#include "Unit.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
  Update hyperparameters

  @param model0 initial values for hyperparameters
  @param newModel new Model object -- will be modified on return
  @param parms analysis settings

  upon return input newUnit is updated with new hyperparameters
*/

void updateHyperParams(const Model& model0,
                       Model& newModel,
                       const Parameters& parms) {

  // number of runs contained in previous time step
  // (okay because newModel has not yet been modified)
  int nRuns = newModel.nu.n_elem;

  mat btmp, val0, val1;

  // pull original model variables -- speed
  mat V0 = model0.V.slice(0);
  mat L0 = model0.Lambda.slice(0);
  mat B0 = model0.B.slice(0);

  mat LB = L0 * B0;

  newModel.Lambda = newModel.XTX.each_slice() + L0;

  newModel.B.set_size(parms.k, parms.d, nRuns);

  newModel.V.set_size(parms.d, parms.d, nRuns);

  for (int i = 0; i < nRuns; ++i) {

    btmp = inv(newModel.Lambda.slice(i)) * (newModel.XTY.slice(i) + LB);

    newModel.B.slice(i) = btmp;

    val0 = newModel.XTY.slice(i).t() * btmp;
    val1 = newModel.YTY.slice(i) - val0 - val0.t() + 
               btmp.t() * newModel.XTX.slice(i) * btmp;
    btmp -= B0;

    newModel.V.slice(i) = V0 + val1 + btmp.t() * L0 * btmp;
  }

  newModel.nu += 1.0;
    
  return;
}
