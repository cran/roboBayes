#include "Model.h"
#include "Parameters.h"
#include "Unit.h"

void stepUK(const arma::uword cp, 
            Unit& newUnit,
            const Unit& prevUnit,
            const Model& model0,
            const int itime,
            const Parameters& parms,
            const arma::mat& datapts,
            const arma::mat& moreData,
            const arma::mat& y_store,
            const arma::mat& x_store);
