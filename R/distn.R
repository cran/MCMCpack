# TO DO
#  igamma
#  standardize both C++ models
#  probit model
#  verbosity in C++ models
#  clean up documentation -- add LaTeX and ASCII models and more explanations


######## Wishart distribution ######## 

# rwish delivers a pseudo-random Wishart deviate
#
# USAGE:
#
#   A <- rwish(v, S)
#
# INPUT:
#
#   v    degrees of freedom
#
#   S    Scale matrix
#
# OUTPUT:
#
#  A     a pseudo-random Wishart deviate
#
# Based on code originally posted by Bill Venables to S-news
# 6/11/1998 -- Written by KQ on 2/5/2001

rwish <- function(v, S) {
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)){
    stop(message="ERROR: S not square in rwish()\n\n")
  }
  if (v < nrow(S)){
    stop(message="ERROR: v is less than the dimension of S in rwish()\n\n")
  }
  p <- nrow(S)
  CC <- chol(S) 
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
  if(p > 1) {
    pseq <- 1:(p-1)
    Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
  }
  crossprod(Z %*% CC)
}

# dwish evaluations the Wishart pdf for positive definite matrix W.
# Uses the Gelman, et. al. parameterization
#
# USAGE:
#
#   x <- dwish(W, v, S)
#
# INPUT:
#
#   W    positive definite matrix at which to evaluate PDF
#
#   v    degrees of freedom
#
#   S    Scale matrix
#
# OUTPUT:
#
#   x    the PDF evaluated (scalar)
#
# Written by ADM on 8/16/2002

dwish <- function(W, v, S) {
   if (!is.matrix(S))
   S <- matrix(S)
   if (nrow(S) != ncol(S)){
      stop(message="ERROR: W not square in dwish()\n\n")
   }
   if (!is.matrix(W))
   S <- matrix(W)
   if (nrow(W) != ncol(W)){
      stop(message="ERROR: W not square in dwish()\n\n")
   }   
   if(nrow(S) != ncol(W)){
      stop(message="ERROR: W and X of different dimensionality in dwish()\n\n")
   }
   if (v < nrow(S)){
      stop(message="ERROR: v is less than the dimension of S in  dwish()\n\n")
  }    
  k <- nrow(S)
  
  # denomonator
  gammapart <- 1
  for(i in 1:k) {
     gammapart <- gammapart * gamma((v + 1 - i)/2)
  } 
  denom <- gammapart *  2^(v * k / 2) * pi^(k*(k-1)/4)
  
  # numerator
  detS <- det(S)
  detW <- det(W)
  hold <- solve(S) %*% W
  tracehold <- sum(hold[row(hold) == col(hold)])  
  num <- detS^(-v/2) * detW^((v - k - 1)/2) * exp(-1/2 * tracehold)

  num / denom
}


######## inverse Wishart distribution ######## 

# riwish generates a draw from the inverse Wishart distribution (using
# the Wishart generator).  
riwish <- function(v, S) {
   solve(rwish(v,S))
}

# diwish evaluations the inverse Wishart pdf for positive definite
# matrix W.  Uses the Gelman, et. al. parameterization
#
# USAGE:
#
#   x <- diwish(W, v, S)
#
# INPUT:
#
#   W    positive definite matrix at which to evaluate PDF
#
#   v    degrees of freedom
#
#   S    Scale matrix
#
# OUTPUT:
#
#   x    the PDF evaluated (scalar)
#
# Written by ADM on 8/16/2002
 
diwish <- function(W, v, S) {
   if (!is.matrix(S))
   S <- matrix(S)
   if (nrow(S) != ncol(S)){
      stop(message="ERROR: W not square in diwish()\n\n")
   }
   if (!is.matrix(W))
   S <- matrix(W)
   if (nrow(W) != ncol(W)){
      stop(message="ERROR: W not square in diwish()\n\n")
   }   
   if(nrow(S) != ncol(W)){
      stop(message="ERROR: W and X of different dimensionality in diwish()\n\n")
   }
   if (v < nrow(S)){
      stop(message="ERROR: v is less than the dimension of S in  diwish()\n\n")
  }    

  # Gelman, et. al. parameterization
  
  k <- nrow(S)   

  # denomonator
  gammapart <- 1
  for(i in 1:k) {
     gammapart <- gammapart * gamma((v + 1 - i)/2)
  } 
  denom <- gammapart *  2^(v * k / 2) * pi^(k*(k-1)/4)
  
  # numerator
  detS <- det(S)
  detW <- det(W)
  hold <- S %*% solve(W)
  tracehold <- sum(hold[row(hold) == col(hold)])  
  num <- detS^(v/2) * detW^(-(v + k + 1)/2) * exp(-1/2 * tracehold)

  num / denom

}

######## inverse Gamma distribution ######## 



# inverse gamma distribution
#digamma <- function(x, shape, rate = 1) {
#   alpha <- shape
#   beta <- rate   
#   beta^alpha / exp(lgamma(alpha)) * x^(-1*(alpha + 1)) * exp(-beta/x)
#}

#rigamma <- function(n, shape, rate = 1, scale = 1/rate) {
#  1 / rgamma(n, shape, rate, scale)
#}

