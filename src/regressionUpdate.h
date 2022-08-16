#include "Model.h"
#include "Parameters.h"
#include "Unit.h"

void regressionUpdate(const arma::vec& datapt,
                      const arma::vec& covariate,
                      const Model& model0,
                      Model& newModel,
                      const arma::mat& moreData,
                      const Parameters& parms);
