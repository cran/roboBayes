#include "Model.h"
#include "Parameters.h"
#include "Unit.h"

void stepCP(arma::uvec& cp, 
            const int itime,
            Unit& newUnit,
            const Unit& prevUnit,
            const Model& model0,
            const arma::mat& datapts,
            const arma::mat& moreData,
            const Parameters& parms,
            const arma::mat& y_store,
            const arma::mat& x_store);
