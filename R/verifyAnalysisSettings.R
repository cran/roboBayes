.verifyAnalysisSettings <- function(d,
                                    k,
                                    lambda,
                                    par_inits,
                                    truncRthresh,
                                    truncRmin,
                                    cpthresh,
                                    cptimemin,
                                    Lgroup,
                                    Lsearch,
                                    Lwindow,
                                    Lm,
                                    alpha,
                                    kt,
                                    pc,
                                    cp_delay,
                                    outlier_mean,
                                    outlier_var,
                                    getR,
                                    getOutliers,
                                    getModels) {

  params = list()

  # Lwindow must be integer-like and >= 1
  params$Lwindow <- .mustBeInteger(x = Lwindow, 
                                   minInt = 1L, 
                                   name = "Lwindow")

  # Lsearch must be integer-like and >= 1
  params$Lsearch <- .mustBeInteger(x = Lsearch, 
                                   minInt = 1L,
                                   name = "Lsearch")

  # Lgroup must be integer-like and >= 1
  params$Lgroup <- .mustBeInteger(x = Lgroup, 
                                  minInt = 1L,
                                  name = "Lgroup")

  # cp_delay must be integer-like and >= 0
  params$cp_delay <- .mustBeInteger(x = cp_delay, 
                                    minInt = 0L,
                                    name = "cp_delay")
  # warn if cp_delay is > 15
  if (params$cp_delay > 15L) {
    warning("cp_delay is potentially too large")
  }

  # truncRmin must be integer-like and >= 5
  params$truncRmin <- .mustBeInteger(x = truncRmin,  
                                     minInt = 5L,
                                     name = "truncRmin")

  # cptimemine must be integer-like and >= 2
  params$cptimemin <- .mustBeInteger(x = cptimemin,  
                                     minInt = 2L,
                                     name = "cptimemin")

  # truncRthresh must be numeric 0 <= truncRthresh <= 1
  params$truncRthresh <- .mustBeNumeric(x = truncRthresh, 
                                        minNum = 0.0,
                                        maxNum = 1.0,
                                        name = "truncRthresh")

  # all cpthresh must be numeric 0 <= cpthresh <= 1
  params$cpthresh <- .mustBeNumeric(x = cpthresh,  
                                    minNum = 0.0,
                                    maxNum = 1.0,
                                    name = "cpthresh")
  params$cpthresh <- sort(x = unique(x = params$cpthresh))

  # pc must be numeric 0 <= truncRthresh <= 1
  params$pc <- .mustBeNumeric(x = pc,   
                              minNum = 0.0,
                              maxNum = 1.0,
                              name = "pc")

  # lambda must be numeric
  params$lambda <- .mustBeNumeric(x = lambda, name = "lambda")

  # alpha must be numeric 0 <= alpha <= 1
  params$alpha <- .mustBeNumeric(x = alpha, 
                                 minNum = 0.0,
                                 maxNum = 1.0,
                                 name = "alpha")

  # d is the number of observations for each time point
  params$d <- d

  # k is the number of covariates for each time point
  params$k <- k

  # kt is the number of time dependent covariates for each time point
  params$kt <- .mustBeInteger(x = kt, 
                              minInt = 0,
                              maxInt = params$k,
                              name = "kt")

  params$outlier_mean <- .mustBeNumeric(x = outlier_mean, name = "outlier_mean")
  if (length(x = params$outlier_mean) != params$d) {
    stop("outlier_mean is of inappropriate length", call. = FALSE)
  }

  params$outlier_var <- .mustBeNumericMatrix(x = outlier_var, 
                                             minVal = 0.0,
                                             name = "outlier_var")
  if ((nrow(x = params$outlier_var) != params$d) ||
      (ncol(x = params$outlier_var) != params$d)) {
    stop("outlier_var is of inappropriate dimension", call. = FALSE)
  }

  # getR, getOutliers, getModels must be logical
  params$getR <- .mustBeLogical(x = getR, name = "getR")
  params$getOutliers <- .mustBeLogical(x = getOutliers, name = "getOutliers")
  params$getModels <- .mustBeLogical(x = getModels, name = "getModels")

  if (params$getModels) {
    params$Lm <- params$Lwindow + params$cp_delay + 1
  } else {
    params$Lm <- .mustBeInteger(x = Lm, name = "Lm")
  }

  # Ensure expected relationships

  tst <- params$Lsearch + params$cp_delay + params$Lgroup - 1L <= params$truncRmin
  if (!tst) {
    stop("(Lsearch + cp_delay + Lgroup - 1 <= truncRmin) is not satisfied; ",
         params$Lsearch + params$cp_delay + params$Lgroup - 1, ">",
         params$truncRmin, call. = FALSE)
  }

  tst <- params$Lwindow + params$cp_delay + params$Lgroup - 1L <= params$truncRmin
  if (!tst) {
    stop("(Lwindow + cp_delay + Lgroup -1 <= truncRmin) is not satisfied; ",
         params$Lwindow + params$cp_delay + params$Lgroup - 1, ">",
         params$truncRmin, call. = FALSE)
  }

  tst <- params$Lgroup <= params$Lwindow + 1L
  if (!tst) {
    stop("(Lgroup <= Lwindow + 1) is not satisfied; ",
         params$Lgroup, ">", params$Lwindow+1, call. = FALSE)
  }

  tst <- params$cptimemin > params$Lgroup
  if (!tst) {
    stop("(cptimemin > Lgroup) is not satisfied; ",
         params$cptimemine, "<=", params$Lgroup, call. = FALSE)
  }

  tst <- params$cptimemin > params$cp_delay-1
  if (!tst) {
    stop("(cptimemin >= cp_delay - 1) is not satisfied; ", 
         params$cptimemin, "<", params$cp_delay-1, call. = FALSE)
  }

  tst <- params$cp_delay < params$Lm
  if (!tst) {
    stop("(cp_delay < Lm) is not satisfied; ", 
         params$cp_delay, ">=", params$Lm, call. = FALSE)
  }
  
  if (params$Lsearch > params$truncRmin) {
    stop("search window too large; increase truncRmin or decrease Lsearch", 
         call. = FALSE)
  }


  return( params )
}
