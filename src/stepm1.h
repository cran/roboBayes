#include "Model.h"
#include "Parameters.h"
#include "Unit.h"

void stepm1(const arma::uword tcp, 
            const int itime,
            const Unit& prevUnit,
            const Unit& newUnit,
            Unit& currUnit,
            const arma::uvec& validpts,
            const Model& model0,
            const arma::mat& moreData,
            const Parameters& parms,
            const arma::mat& y_store,
            const arma::mat& x_store);
