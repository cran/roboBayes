#' roboBayes: Robust Online Bayesian Monitoring
#'
#' The primary function for Robust Online Bayesian Monitoring (roboBayes). 
#'   Performs roboBayes on the response data regressed on the input covariates, 
#'   which should account for nonstationarity in the response between 
#'   changepoints. The algorithm outputs changepoints based on a threshold, a 
#'   roboBayes object incorporating data thus far, and algorithm settings. A 
#'   roboBayes object can be updated with new data by including it in a successive 
#'   roboBayes() function call.
#'
#' All inputs after the ellipsis are analysis settings that can only be set
#'   in the initial step of the roboBayes procedure. In the first step, input
#'   roboBayes must be NULL. All subsequent steps of the procedure must provide,
#'   through input roboBayes, the value object returned by the previous step.
#'
#' Analysis parameters must obey the following conditions:
#'   \itemize{
#'     \item Lsearch + cp_delay + Lgroup <= truncRmin + 1
#'     \item Lwindow + cp_delay + Lgroup <= truncRmin + 1
#'     \item Lgroup <= Lwindow + 1
#'     \item cptimemin > Lgroup
#'     \item cptimemin >= cp_delay - 1
#'     \item cp_delay < Lm
#'     \item Lsearch <= truncRmin
#'   }
#'
#' @param datapts A matrix object (n x d). The observations for n time points.
#'   There are methods in place to attempt to convert non-matrix input 
#'   (data.frame, numeric) to a matrix. If a numeric vector is provided, 
#'   it is assumed that d = 1.
#'
#' @param covariates A matrix object (n x k). The covariates for n time points.
#'   If any covariates depend on the run length (kt > 0), these covariates
#'   must be provided in the last kt columns of the input matrix.
#'   There are methods in place to attempt to convert non-matrix input 
#'   (data.frame, numeric) to a matrix. If a numeric vector is provided, 
#'   it is assumed that k = 1.
#'
#' @param RoboBayes A value object returned by a prior roboBayes() analysis or NULL.
#'   If NULL, the analysis is taken to be the first step, and all analysis
#'   settings (inputs after the ellipsis) can be adjusted. For all subsequent 
#'   steps, RoboBayes is the value object of the preceding step and analysis 
#'   settings are taken from that object -- any analysis setting inputs provided
#'   for continuation steps are ignored.
#'
#' @param ... Ignored. Included only to require named inputs. If RoboBayes is an
#'   object of class RoboBayes, all inputs after the ellipsis are ignored.
#' 
#' @param lambda A numeric object. The prior probability of a changepoint
#'   occurring at any time point.
#'
#' @param par_inits A list object. The initial estimates for the hyperparameters.
#'   The list can contain one or more of the following:
#'   \itemize{
#'   \item B: A (k x d) matrix. The location parameter of the matrix Normal.
#'            The default is a zero matrix.
#'
#'   \item V: A (d x d) matrix. The scale matrix of the Inverse Wishart.
#'            The default is a symmetric matrix with 1.0 on the diagonal and
#'            0.1 in the off-diagonal elements.
#'
#'   \item nu: A scalar. The degrees of freedom of the Inverse Wishart.
#'             The default is (d - 0.9).
#'
#'   \item Lambda: A (k x k) matrix. The scale parameter of the matrix Normal.
#'                 The default is a diagonal matrix of 0.01.}
#'
#' @param truncRthresh A scalar object. A probability threshold used to 
#'   limit the size of the data used at each time step. Once t > truncRmin,
#'   run lengths associated with data points from times less than t-truncRmin 
#'   for which the probability < truncRthresh are removed 
#'
#' @param truncRmin An integer object. The minimum number of time points 
#'   maintained before data truncation methods are utilized.
#'
#' @param cpthresh A numeric vector object. One or more changepoints.
#'
#' @param cptimemin An integer object. Minimum time at which updated change
#'   points can be ascertained.
#'
#' @param Lgroup An integer object. The number of runs grouped together to
#'   identify a new changepoint.
#'
#' @param Lgroup An integer object. The size of the moving window used to 
#'   aggregate run length probabilities during the initial coarse search for 
#'   the presence of a changepoint.
#'
#' @param Lsearch An integer object. The number of run length probabilities to 
#'   aggregate for calculating the probability of "recent" change.
#'
#' @param Lwindow An integer object. The number of points in the window to scan
#'   for the most recent changepoint. Changes further back than Lwindow + 
#'   cp_delay will not be recorded.
#'
#' @param Lm An integer object. The number of data points to store in memory
#'   for outlier elimination. A change recorded more than Lm points previously
#'   will not be examined for being an outlier. When getModels is TRUE, 
#'   Lm is automatically set to  Lwindow + cp_delay + 1 to ensure that model
#'   parameters associated with detected changepoints can be recovered.
#'
#' @param alpha A scalar object. The probability threshold for declaring a 
#'   point an outlier.
#'
#' @param kt An integer object. The number of covariates that depend on run
#'   length.
#'
#' @param pc A scalar object. The prior probability of the model without an
#'   outlier being the correct model, given that a candidate change/outlier has
#'   been detected.
#'
#' @param cp_delay An integer object. The minimum number of points that must
#'   be observed after the change to declare a changepoint.
#'
#' @param outlier_mean A numeric vector (d). The mean value for the 
#'   distribution of outlier values.
#'
#' @param outlier_var A matrix object (d x d). The covariance matrix for the 
#'   distribution of outlier values.
#'
#' @param getR A logical object. If TRUE, the matrix of probabilities is
#'   returned. Warning -- selecting TRUE will increase computation time, and
#'   for large analyses, may exceed available memory.
#'
#' @param getOutliers A logical object. Detect outliers. If TRUE, the matrix of
#'   probabilities used for outlier elimination are returned. Selecting TRUE
#'   will increase computation time.
#'
#' @param getModels A logical object. If TRUE, the model parameters of each
#'   changepoint are returned, as well as the matrix of probabilities used
#'   for outlier elimination. When TRUE, parameter Lm is set to 
#'   Lwindow + cp_delay + 1, regardless of that provided as input.
#'
#' @returns A RoboBayes object, which extends list and contains
#'   \itemize{
#'   \item R: A numeric vector object. The current posterior distribution of 
#'            the run length.
#'   \item RL: An integer vector object. The current run lengths that are 
#'             retained.
#'   \item truncRind: An integer vector object. The indices of the run 
#'                    lengths retained.
#'   \item jtR: A numeric vector object. The current joint distribution of the 
#'              run length and data.
#'   \item pars: A list object. Elements are run length specific summary 
#'               statistics and hyperparameters. With the exception of nu, 
#'               each hyperparameter/summary statistic is a 3 dimensional array, 
#'               where the final dimension corresponds to the run length. The 
#'               degrees of freedom, nu, are returned as a vector, of length 
#'               equal to the number of run lengths retained.
#'   \item cpInds: An integer matrix object (length(cpthresh) x nRuns). Each 
#'                 element contains the most recent changepoint.
#'   \item lastLs: A numeric vector object. Each element contains the 
#'                 probability that a change occurred in the previous Lsearch 
#'                 time points, delayed by cp_delay points.
#'   \item time: An integer object. The current time.
#'   \item allcov: A numeric matrix object. The run length dependent covariates
#'                 for the retained run lengths.
#'   \item model0: A list object. The initial hyperparameters.
#'   \item lastDataPt: A numeric vector object. The data of the last time point.
#'   \item call: The matched call.
#'   \item params: A list object. The analysis settings.
#'   }
#' Conditionally returned elements include:
#' 
#' If getR = TRUE
#'   \itemize{
#'   \item RFull: A numeric matrix object. The ith column contains the run 
#'                length posterior distribution for time i.
#'   }
#' If getOutliers = TRUE
#'   \itemize{
#'     \item Rm: A numeric matrix object. The ith column contains the joint 
#'              distribution of data and run length associated with time-i+1.
#'     \item y_store: A numeric matrix object. Contains the most recent Lm data 
#'                    points.
#'     \item x_store: A numeric matrix object. Contains the most recent Lm 
#'                    covariates.
#'     \item outliers: An integer vector object. Time points that have been 
#'                     identified as outliers.
#'   }
#' If getModels = TRUE
#'   \itemize{
#'     \item Rm: A numeric matrix object. The ith column contains the joint 
#'              distribution of data and run length associated with time-i+1.
#'     \item y_store: A numeric matrix object. Contains the most recent Lm data 
#'                    points.
#'     \item x_store: A numeric matrix object. Contains the most recent Lm 
#'                    covariates.
#'     \item mods: A list object. Each element of the list corresponds to an 
#'                 identified changepoint. For each changepoint, the expected 
#'                 model coefficients and covariance. 
#'     \item currentModel: A list object. For the current most likely run length, 
#'                    the expected model coefficients and covariance.
#'   }
#'
#'
#' @examples
#' \dontshow{
#'   RcppArmadillo::armadillo_throttle_cores(2)
#' }
#'
#'
#' nt <- 100
#'
#' ## 2 covariates each time step
#'
#' x <- cbind(rep(1,nt), rnorm(nt))
#'
#' ## 2 observations each time step
#'
#' # covariance matrix
#' sigma <- matrix(data = -0.3, nrow = 2, ncol = 2)
#' diag(sigma) <- 1.0
#'
#' # mean
#' mean1 <- c(0.5, 0.1)
#' mean2 <- c(0.2, 0.2)
#'
#' y <- rbind(mvtnorm::rmvnorm(n = nt/2, mean = mean1, sigma = 0.001*sigma),
#'            mvtnorm::rmvnorm(n = nt/2, mean = mean2, sigma = 0.001*sigma))
#'
#' # add in an outlier
#' y[30,] <- c(0.9,0.01)
#' y[y < 0.0] <- 0.0
#'
#' step1 <- roboBayes(datapts = y[1:32,], covariates = x[1:32,])
#'
#' step2 <- roboBayes(datapts = y[33:64,], covariates = x[33:64,], RoboBayes = step1)
#'
#' step3 <- roboBayes(datapts = y[65:94,], covariates = x[65:94,], RoboBayes = step2)
#'
#' step4 <- roboBayes(datapts = y[95:100,], covariates = x[95:100,], RoboBayes = step3)
#'
#' @export
#' @useDynLib roboBayes
#' @import Rcpp
#' @importFrom RcppArmadillo armadillo_throttle_cores
#' @include checkTypes.R verifyAnalysisSettings.R initRoboBayes.R regressionInit.R
#' @include RcppExports.R

roboBayes <- function(datapts,
                      covariates,
                      RoboBayes = NULL, 
                      ...,
                      lambda = 1000.0,
                      par_inits = NULL,
                      truncRthresh = 1e-4,
                      truncRmin = 100L,
                      cpthresh = 0.8,
                      cptimemin = 4L,
                      Lgroup = 3L,
                      Lsearch = 10L,
                      Lwindow = 30L,
                      Lm = 7L,
                      alpha = 0.9,
                      kt = 0L,
                      pc = 0.5,
                      cp_delay = 3L,
                      outlier_mean = rep(x = 0.5, times = ncol(x = datapts)),
                      outlier_var = diag(x = 2.0, nrow = ncol(x = datapts)),
                      getR = FALSE,
                      getOutliers = TRUE,
                      getModels = FALSE) {

  matchedCall <- match.call()

  # ensure that user did not misspell an input or provide
  # unexpected inputs
  if (length(list(...)) != 0L) {
    stop("unrecognized input provided ",
         paste(names(list(...)), collapse = ", "))
  }

  # datapts must always be provided as input
  if (is.null(x = datapts)) stop("datapts must be provided")

  # ensure that datapts is a matrix. If provided a data.frame, it is
  # converted to a data.matrix; if provided a vector, it is converted
  # to a one column matrix
  datapts <- .checkMatrix(x = datapts, name = "datapts")

  if (nrow(x = datapts) == 0L) {
    return( RoboBayes )
  }
  
  # covariates must always be provided as input
  if (is.null(x = covariates)) stop("covariates must be provided")

  # ensure that covariates is a matrix. If provided a data.frame, it is
  # converted to a data.matrix; if provided a vector, it is converted
  # to a one column matrix
  covariates <- .checkMatrix(x = covariates, name = "covariates")

  if (is.null(x = RoboBayes)) {
    # RoboBayes is null for the first step of the procedure. 

    # Ensure appropriate inputs for analysis settings.
    params <- .verifyAnalysisSettings(d = ncol(x = datapts),
                                      k = ncol(x = covariates),
                                      lambda = lambda,
                                      par_inits = par_inits,
                                      truncRthresh = truncRthresh,
                                      truncRmin = truncRmin,
                                      cpthresh = cpthresh,
                                      cptimemin = cptimemin,
                                      Lgroup = Lgroup,
                                      Lsearch = Lsearch,
                                      Lwindow = Lwindow,
                                      Lm = Lm,
                                      alpha = alpha,
                                      kt = kt,
                                      pc = pc,
                                      cp_delay = cp_delay,
                                      outlier_mean = outlier_mean,
                                      outlier_var = outlier_var,
                                      getR = getR,
                                      getOutliers = getOutliers,
                                      getModels = getModels)

    # verify or initialize initial estimates for hyper parameters
    if (is.list(x = par_inits)) {
      RoboBayes <- initRoboBayes(par_inits = par_inits, 
                                 params = params,
                                 nTimes = nrow(x = datapts))
    }

  } else if (methods::is(object = RoboBayes, class2 = "RoboBayes")) {
    # RoboBayes object indicates a continuation step. Take analysis settings
    # from previous step object
    params <- RoboBayes$params

    # ensure that data has proper number of observations for each time step
    if (ncol(x = datapts) != params$d) {
      stop("datapts is of inappropriate dimension")
    }

    # ensure that covariates has proper number of observations for each time step
    if (ncol(x = covariates) != params$k) {
      stop("covariates is of inappropriate dimension")
    }

    # remove parameters from previous object
    RoboBayes$params <- NULL

    if (params$getR) {
      # if user requested to return matrix of posterior distributions for
      # each time, extend previous result to include space for new time points

      numpts <- RoboBayes$time
  
      # number of new data points
      numNewPts <- nrow(x = datapts)

      RoboBayes$RFull <- cbind(RoboBayes$RFull, 
                               matrix(data = 0.0, 
                                      nrow = nrow(x = RoboBayes$RFull), 
                                      ncol = numNewPts))
      RoboBayes$RFull <- rbind(RoboBayes$RFull,
                               matrix(data = 0.0,
                                      nrow = numNewPts,
                                      ncol = ncol(x = RoboBayes$RFull)))
    }

  } else {
    # RoboBayes is neither NULL nor a RoboBayes object, throw error
    stop("RoboBayes is not of a recognized class")
  }

  if (is.null(x = RoboBayes)) {

    # if RoboBayes is still NULL, par_inits was not provided. Set this
    # up in the C++ code
    res <- roboBayesFirst(datapts = datapts,
                          covariates = covariates,
                          paramsList = params)

  } else if (methods::is(object = RoboBayes, class2 = "RoboBayes")) {

    # if RoboBayes is a RoboBayes object, either par_inits was provided in first
    # step or this is a continuation step

    RoboBayes <- unclass(x = RoboBayes)

    nt <- length(x = RoboBayes$mods)
    modelCP <- numeric(length = nt)
    modelB <- array(data = 0.0, dim = c(params$k, params$d, nt))
    modelV <- array(data = 0.0, dim = c(params$d, params$d, nt))
    if (nt > 0) {
      for (i in 1L:nt) {
        modelCP[i] <- RoboBayes$mods[[ i ]]$cp
        modelB[,,i] <- RoboBayes$mods[[ i ]]$B
        modelV[,,i] <- RoboBayes$mods[[ i ]]$V
      }
    }

    RoboBayes$mods <- list("cp" = modelCP, "B" = modelB, "V" = modelV)

    res <- roboBayesPrior(datapts = datapts,
                          covariates = covariates,
                          RoboBayes = RoboBayes,
                          paramsList = params)
  } else {
    # if neither condition is true, logic isn't working as expected.
    stop("unexpected RoboBayes object; contact developer")
  }

  if (params$getModels) {
    # extract the model for the most likely run length
    imax <- which.max(x = res[["R"]])
    deno <- res$pars$nu[imax] - params$d - 1.0
    res$currentModel <- list("B" = res$pars$B[,,imax],
                             "V" = res$pars$V[,,imax]/deno)

    # redefine the models for CP to conform to original definition
    if (length(res$modelCP) > 0L) {
      res$mods <- list()

      for (i in 1L:length(res$modelCP)) {
        res$mods[[ i ]] <- list("cp" = res$modelCP[i],
                                "B" = matrix(data = res$modelB[,,i],
                                             nrow = params$k,
                                             ncol = params$d),
                                "V" = matrix(data = res$modelV[,,i],
                                             nrow = params$d,
                                             ncol = params$d))
      }
    }
  }

  if (!params$getR) {
    # if posterior distribution for all time points not requested,
    # remove place holder from returned object
    res$RFull <- NULL
  }

  if (!params$getOutliers && !params$getModels) {
    # if outliers or models not requested, remove place holders from 
    # returned object
    res$Rm <- NULL
    res$y_store <- NULL
    res$x_store <- NULL
  }

  if (!params$getOutliers) {
    # if outliers not requested, remove from returned object
    res$outliers <- NULL
  }

  # remove changepoint model information from returned object
  res$modelCP <- NULL
  res$modelB <- NULL
  res$modelV <- NULL

  res$lastDataPt <- datapts[nrow(x = datapts),]
  res$call <- matchedCall
  res$params <- params
 
  class(x = res) <- c("RoboBayes", class(x = res))

  return( res )
}

