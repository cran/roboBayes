#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <RcppArmadillo.h>

using namespace Rcpp;

/*
  Class to store analysis settings
*/

class Parameters {

  public: 

  double alpha;

  int Lwindow;
  int Lsearch;
  int Lgroup;
  int Lm;

  int truncRmin;
  double truncRthresh;

  arma::vec cpthresh;
  int cp_delay;
  int cptimemin;

  bool getR;
  bool getOutliers;
  bool getModels;

  double pc;
  double lambda;

  int k;
  int kt;
  int d;

  arma::vec outlier_mean;
  arma::mat outlier_var;

  Parameters(){}

  Parameters(List params) : alpha{as<double>(params["alpha"])},
                            Lwindow{as<int>(params["Lwindow"])},
                            Lsearch{as<int>(params["Lsearch"])},
                            Lgroup{as<int>(params["Lgroup"])},
                            Lm{as<int>(params["Lm"])},
                            truncRmin{as<int>(params["truncRmin"])},
                            truncRthresh{as<double>(params["truncRthresh"])},
                            cpthresh{as<arma::vec>(params["cpthresh"])},
                            cp_delay{as<int>(params["cp_delay"])},
                            cptimemin{as<int>(params["cptimemin"])},
                            getR{as<bool>(params["getR"])},
                            getOutliers{as<bool>(params["getOutliers"])},
                            getModels{as<bool>(params["getModels"])},
                            pc{as<double>(params["pc"])},
                            lambda{as<double>(params["lambda"])},
                            k{as<int>(params["k"])},
                            kt{as<int>(params["kt"])},
                            d{as<int>(params["d"])},
                            outlier_mean{as<arma::vec>(params["outlier_mean"])},
                            outlier_var{as<arma::mat>(params["outlier_var"])} {}

};

#endif
