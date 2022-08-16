# Ensure that input is a matrix object
#
# Ensures that the tested object is returned as a matrix
# non-error methods are defined for data.frame, numeric, and matrix
#
# @param x An object of class ANY, data.frame, numeric, or matrix
#
# @param name A character object. The input variable name. Used for 
#             error printing
#
#' @import methods
setGeneric(name = ".checkMatrix",
           def = function(x, ...) { standardGeneric(".checkMatrix") })

setMethod(f = ".checkMatrix",
          signature = c(x = "ANY"),
          definition = function(x, ..., name) { 
              stop(name, " is not a appropriate object", call. = FALSE) 
            })

setMethod(f = ".checkMatrix",
          signature = c(x = "data.frame"),
          definition = function(x, ...) { 
              .checkMatrix(x = data.matrix(frame = x), ...)
            })

setMethod(f = ".checkMatrix",
          signature = c(x = "numeric"),
          definition = function(x, ...) { 
              return( matrix(data = x, ncol = 1L) )
            })

setMethod(f = ".checkMatrix",
          signature = c(x = "matrix"),
          definition = function(x, ...) { 
              return( x )
            })

# Ensure that input is an integer object
#
# Ensures that the tested object is integer-like and returned as an
# integer object (cannot be a matrix with column > 1)
#
# @param x An object of class ANY, matrix, numeric, or integer
#
# @param name A character object. The input variable name. Used for 
#             error printing
#
setGeneric(name = ".mustBeInteger",
           def = function(x, ...) { standardGeneric(".mustBeInteger") })

setMethod(f = ".mustBeInteger",
          signature = c(x = "ANY"),
          definition = function(x, ..., name) { 
              stop(name, " is not an appropriate object", call. = FALSE) 
            })

setMethod(f = ".mustBeInteger",
          signature = c(x = "matrix"),
          definition = function(x, ...) { 
              if (ncol(x = x) == 1L) {
                return( .mustBeInteger(x = drop(x = x), ...) )
              } else {
                return( .mustBeInteger(x = NULL, ...) )
              }
            })

setMethod(f = ".mustBeInteger",
          signature = c(x = "numeric"),
          definition = function(x, ...) { 
              ix <- as.integer(x = round(x = x, digits = 0L))
              if (!isTRUE(x = all.equal(target = x, current = ix))) {
                return( .mustBeInteger(x = NULL, ...) )
              } else {
                return( .mustBeInteger(x = ix, ...) )
              }
            })

setMethod(f = ".mustBeInteger",
          signature = c(x = "integer"),
          definition = function(x, ..., minInt = -Inf, maxInt = Inf, name) { 
              if (any(x < minInt) || any(x > maxInt)) {
                stop(name, " is outside of allowed value range")
              }
              return( x )
            })

# Ensure that input is a numeric object
#
# Ensures that the tested object is numeric and returned as a
# numeric object (cannot be a matrix with column > 1)
#
# @param x An object of class ANY, matrix, or numeric
#
# @param name A character object. The input variable name. Used for 
#             error printing
#
setGeneric(name = ".mustBeNumeric",
           def = function(x, ...) { standardGeneric(".mustBeNumeric") })

setMethod(f = ".mustBeNumeric",
          signature = c(x = "ANY"),
          definition = function(x, ..., name) { 
              stop(name, " is not an appropriate object", call. = FALSE) 
            })

setMethod(f = ".mustBeNumeric",
          signature = c(x = "matrix"),
          definition = function(x, ...) { 
              if (ncol(x = x) == 1L) {
                return( .mustBeNumeric(x = drop(x = x), ...) )
              } else {
                return( .mustBeNumeric(x = NULL, ...) )
              }
            })

setMethod(f = ".mustBeNumeric",
          signature = c(x = "numeric"),
          definition = function(x, ..., minNum = -Inf, maxNum = Inf, name) { 
              if (any(x < minNum) || any(x > maxNum)) {
                stop(name, " is outside of allowed value range")
              }
              return( x )
            })

# Ensure that input is a numeric object
#
# Ensures that the tested object is numeric and returned as a
# numeric object (cannot be a matrix with column > 1)
#
# @param x An object of class ANY, matrix, or numeric
#
# @param name A character object. The input variable name. Used for 
#             error printing
#
setGeneric(name = ".mustBeNumericMatrix",
           def = function(x, ...) { standardGeneric(".mustBeNumericMatrix") })

setMethod(f = ".mustBeNumericMatrix",
          signature = c(x = "ANY"),
          definition = function(x, ..., name) { 
              stop(name, " is not an appropriate object", call. = FALSE) 
            })

setMethod(f = ".mustBeNumericMatrix",
          signature = c(x = "matrix"),
          definition = function(x, ..., minNum = -Inf, maxNum = Inf, name) { 
              if (any(x < minNum) || any(x > maxNum)) {
                stop(name, " is outside of allowed value range")
              }
              if (any(is.integer(x = x))) {
                x <- x*1.0
              }
              return( x )
            })

setMethod(f = ".mustBeNumericMatrix",
          signature = c(x = "numeric"),
          definition = function(x, ...) { 
              return( .mustBeNumericMatrix(x = matrix(data = x, ncol = 1L), ...) )
            })

# Ensure that input is a logical object
#
# Ensures that the tested object is logical and returned as a
# logical object
#
# @param x An object of class ANY and logical
#
# @param name A character object. The input variable name. Used for 
#             error printing
#
setGeneric(name = ".mustBeLogical",
           def = function(x, ...) { standardGeneric(".mustBeLogical") })

setMethod(f = ".mustBeLogical",
          signature = c(x = "ANY"),
          definition = function(x, ..., name) { 
              stop(name, " is not an appropriate object", call. = FALSE) 
            })

setMethod(f = ".mustBeLogical",
          signature = c(x = "logical"),
          definition = function(x, ...) { 
              return( x )
            })
