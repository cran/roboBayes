#include "Model.h"

using namespace arma;
using namespace Rcpp;

/* 
  Methods for creating and modifying object of class Model
*/

Model::Model(Parameters& parms) {

  B = zeros<cube>(parms.k, parms.d, 1);

  V = zeros<cube>(parms.d, parms.d, 1);
  mat temp(parms.d, parms.d);
  temp += 0.1;
  temp.diag().ones();
  V.slice(0) = temp;

  temp.zeros(parms.k, parms.k);
  temp.diag().fill(0.01);
  Lambda.slice(0) = temp;

  XTX = zeros<cube>(parms.k, parms.k, 1);
  XTY = zeros<cube>(parms.k, parms.d, 1);
  YTY = zeros<cube>(parms.d, parms.d, 1);

  nu = zeros(1);
  nu(0) = parms.d - 1.0 + 0.1;

}

Model::Model(List& pars) {

  B = as<cube>(pars["B"]);
  V = as<cube>(pars["V"]);
  Lambda = as<cube>(pars["Lambda"]);
  XTX = as<cube>(pars["XTX"]);
  XTY = as<cube>(pars["XTY"]);
  YTY = as<cube>(pars["YTY"]);
  nu = as<vec>(pars["nu"]);

}

List Model::toList() {

  List toReturn = List::create(Named("B", B),
                               Named("V", V),
                               Named("Lambda", Lambda),
                               Named("nu", nu),
                               Named("XTX", XTX),
                               Named("XTY", XTY),
                               Named("YTY", YTY));

  return toReturn;
}

void Model::join(const Model mod) {

  B = join_slices(mod.B.slice(0), B);
  V = join_slices(mod.V.slice(0), V);
  Lambda = join_slices(mod.Lambda.slice(0), Lambda);
  XTX = join_slices(mod.XTX.slice(0), XTX);
  XTY = join_slices(mod.XTY.slice(0), XTY);
  YTY = join_slices(mod.YTY.slice(0), YTY);
  nu = join_cols(mod.nu.row(0), nu);

  return;
}

void Model::subset(uvec idx) {

  Model toReturn;

  B = B.slices(idx);
  V = V.slices(idx);
  Lambda = Lambda.slices(idx);
  XTX = XTX.slices(idx);
  XTY = XTY.slices(idx);
  YTY = YTY.slices(idx);

  nu = nu.elem(idx);

  return;

}
