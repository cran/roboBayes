#include "Parameters.h"
#include "Model.h"

Model spikeUpdate(const arma::mat& datapt,
                  const arma::mat& covariate,
                  const arma::uword tcp,
                  const Model& model0,
                  const Model& modelNew,
                  const arma::mat& moreData,
                  const Parameters& parms);
