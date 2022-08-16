#include "Model.h"
#include "Parameters.h"

arma::vec regressionProb(const arma::vec& datapt,
                         const Model& prevModel,
                         const Model& newModel,
                         const Parameters& parms);
