# Initialize vectors for gaussian probability functions
#
# Takes in the desired initialization parameters,
# initializes the vectors needed for the gaussian probability
# function \code{gaussian_update}
#
# @param pars A list of parameters to be used for
#   initialization
# 
#
# @return A list object containing the parameters to be used in the 
#   iteratively updating algorithm of parameters describing the
#   underlying gaussian distribution of the data.
#
.regressionInit <- function(pars, params) {

  # {k x d x 1L} array of coefficients
  if (is.null(x = pars$B)) {
    pars$B <- matrix(data = 0.0, nrow = params$k, ncol = params$d)
  } else {

    if (!is.matrix(x = pars$B)) {
      stop("B must be a {k x d} matrix", call. = FALSE)
    }

    if ({nrow(x = pars$B) != params$k} || 
        {ncol(x = pars$B) != params$d}) {
      stop("B is of inappropriate dimension", call. = FALSE)
    }

  }
  pars$B <- array(data = pars$B, dim = c(params$k, params$d, 1L))

  # {d x d x 1L} covariance matrix prior
  if (is.null(x = pars$V)) {

    pars$V <- matrix(data = 0.1, nrow = params$d, ncol = params$d)

    diag(x = pars$V) <- 1.0

  } else {

    if (!is.matrix(x = pars$V)) {
      stop("V must be a {d x d} matrix", call. = FALSE)
    }

    if ({nrow(x = pars$V) != params$d} || 
        {ncol(x = pars$V) != params$d}) {
      stop("V is of inappropriate dimension", call. = FALSE)
    }

  }
  pars$V <- array(data = pars$V, dim = c(params$d, params$d, 1L))

  # {k x k x 1L} covariance prior
  if (is.null(x = pars$Lambda)) {

    pars$Lambda <- matrix(data = 0.0, nrow = params$k, ncol = params$k)
    diag(x = pars$Lambda) <- 0.01

  } else {

    if (!is.matrix(x = pars$Lambda)) {
      stop("Lambda must be a {k x k} matrix", call. = FALSE)
    }

    if ({nrow(x = pars$Lambda) != params$k} || 
        {ncol(x = pars$Lambda) != params$k}) {
      stop("Lambda is of inappropriate dimension", call. = FALSE)
    }

  }
  pars$Lambda <- array(data = pars$Lambda, dim = c(params$k, params$k, 1L))

  # scalar df prior > d-1
  if (is.null(x = pars$nu)) {
    pars$nu <- params$d - 1.0 + 0.1
  } else {
    if (pars$nu < 1e-8) {
      stop("nu cannot be <= 0", call. = FALSE)
    }
  }
  
  pars$XTX <- array(data = 0.0, dim = c(params$k, params$k, 1L))
  pars$XTY <- array(data = 0.0, dim = c(params$k, params$d, 1L))
  pars$YTY <- array(data = 0.0, dim = c(params$d, params$d, 1L))
  
  return( pars )
}
