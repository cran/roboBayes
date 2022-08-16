# @param par_inits A list object containing hyperparameter information
#
# @param nChangePoints An integer object. The number of change points.
#
#' @include regressionInit.R
initRoboBayes <- function(par_inits, params, nTimes) {

  ## Initialize
  # posterior distribution of run length r_t given data y_1:t
  pars0 <- .regressionInit(pars = par_inits, params = params)
  
  result <- list("jtR" = 1.0, 
                 "R" = 1.0, 
                 "RL" = 1L, 
                 "RFull" = matrix(data = 0.0, nrow = nTimes+1, ncol = nTimes+1),
                 "Rm" = matrix(data = 1.0, nrow = 1, ncol = params$Lm),
                 "truncRind" = 1L,
                 "outliers" = numeric(0),
                 "cpInds" = matrix(data = 1L, 
                                   nrow = length(params$cpthresh), 
                                   ncol = 0),
                 "lastLs" = numeric(0),
                 "models" = array(data = 0.0, dim = c(params$k, params$d, 0)),
                 "y_store" = matrix(data = 0.0, 
                                    nrow = params$Lm, 
                                    ncol = params$d),
                 "x_store" = matrix(data = 0.0, 
                                    nrow = params$Lm, 
                                    ncol = params$k),
                 "allcov" = matrix(data = 0.0, nrow = params$kt, ncol = 0),
                 "time" = 1L, 
                 "pars" = pars0,
                 "model0" = pars0,
                 "mods" = vector(mode = "list", length = 0L))
  result$RFull[1,1] <- 1
  result$Rm[1,1] <- 1

  class(x = result) <- "RoboBayes"
  
  return(result)
}

