########## Utility Functions ##########

# takes a symmetric matrix x and returns lower diagonal
# note: does not check for symmetry
#
# ADM 4/18/2003 

"vech" <-
  function (x) {
    x <- as.matrix(x)
    if (dim(x)[1] != dim(x)[2]) {
      stop("Non-square matrix passed to vech().\n")
    }
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)
  }

# takes vector x and returns an nrow times nrow symmetric matrix
# this will recycle the elements of x as needed to fill the matrix
#
# ADM 4/18/2003
# ADM 11/13/2003 [bug fix]

"xpnd" <-
  function (x, nrow) {
    dim(x) <- NULL
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag=TRUE)] <- 0
    output <- output + t(hold)    
    return(output)
  }
