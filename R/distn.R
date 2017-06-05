##########################################################################
## Density Functions and Random Number Generators
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################



##
## Wishart
##

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
# on 6/11/1998
#
# KQ on 2/5/2001

#' The Wishart Distribution
#'
#' Density function and random generation from the Wishart distribution.
#'
#' The mean of a Wishart random variable with \code{v} degrees of freedom and
#' inverse scale matrix \code{S} is \eqn{vS}.
#'
#' @aliases dwish rwish
#'
#' @name Wishart
#'
#' @param W Positive definite matrix W \eqn{(p \times p)}.
#'
#' @param v Degrees of freedom (scalar).
#'
#' @param S Inverse scale matrix \eqn{(p \times p)}.
#'
#' @return \code{dwish} evaluates the density at positive definite matrix W.
#' \code{rwish} generates one random draw from the distribution.
#'
#' @keywords distribution
#'
#' @examples
#'
#' density <- dwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
#' draw <- rwish(3, matrix(c(1,.3,.3,1),2,2))
#'
NULL

#' @rdname Wishart
#' @export
"rwish" <-
  function(v, S) {
    if (!is.matrix(S))
      S <- matrix(S)
    if (nrow(S) != ncol(S)) {
      stop(message="S not square in rwish().\n")
    }
    if (v < nrow(S)) {
      stop(message="v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
    if(p > 1) {
      pseq <- 1:(p-1)
      Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    return(crossprod(Z %*% CC))
  }

# dwish evaluations the Wishart pdf at positive definite matrix W.
# note: uses the Gelman, et. al. parameterization.
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
# ADM 8/16/2002

#' @rdname Wishart
#' @export
"dwish" <-
  function(W, v, S) {
    if (!is.matrix(S))
      S <- matrix(S)
    if (nrow(S) != ncol(S)){
      stop(message="W not square in dwish()\n\n")
    }
    if (!is.matrix(W))
      W <- matrix(W)
    if (nrow(W) != ncol(W)){
      stop(message="W not square in dwish()\n\n")
    }
    if(nrow(S) != ncol(W)){
      stop(message="W and X of different dimensionality in dwish()\n\n")
    }
    if (v < nrow(S)){
      stop(message="v is less than the dimension of S in  dwish()\n\n")
    }
    k <- nrow(S)

    # denominator
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

    return(num / denom)
  }

##
## Inverse Wishart
##

# riwish generates a draw from the inverse Wishart distribution
# (using the Wishart generator)

#' The Inverse Wishart Distribution
#'
#' Density function and random generation from the Inverse Wishart
#' distribution.
#'
#' The mean of an inverse Wishart random variable with \code{v} degrees of
#' freedom and scale matrix \code{S} is \eqn{(v-p-1)^{-1}S}.
#'
#' @aliases diwish riwish
#'
#' @name InvWishart
#'
#' @param W Positive definite matrix W \eqn{(p \times p)}.
#'
#' @param v Degrees of freedom (scalar).
#'
#' @param S Scale matrix \eqn{(p \times p)}.
#'
#' @return \code{diwish} evaluates the density at positive definite matrix W.
#' \code{riwish} generates one random draw from the distribution.
#'
#' @keywords distribution
#'
#' @examples
#'
#' density <- diwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
#' draw <- riwish(3, matrix(c(1,.3,.3,1),2,2))
#'
NULL

#' @rdname InvWishart
#' @export
"riwish" <-
  function(v, S) {
    return(solve(rwish(v,solve(S))))
  }

# diwish evaluates the inverse Wishart pdf at positive definite
# matrix W.  note: uses the Gelman, et. al. parameterization.
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
## ADM 8/16/2002
## JHP 01/24/2017 updated
## thanks to mathew simpson (themattsimpson@gmail.com)
## for commenting and suggesting a code.

#' @rdname InvWishart
#' @export
"diwish" <-
  function(W, v, S) {
    if (!is.matrix(S))
      S <- matrix(S)
    if (nrow(S) != ncol(S)){
      stop("S not square in diwish().\n")
    }
    if (!is.matrix(W))
      W <- matrix(W)
    if (nrow(W) != ncol(W)){
      stop("W not square in diwish().\n")
    }
    if(nrow(S) != ncol(W)){
      stop("W and X of different dimensionality in diwish().\n")
    }
    if (v < nrow(S)){
      stop("v is less than the dimension of S in  diwish().\n")
    }
    p <- nrow(S)
    gammapart <- sum(lgamma((v + 1 - 1:p)/2))
    ldenom <- gammapart + 0.5*v*p*log(2) + 0.25*p*(p-1)*log(pi)
    cholS <- chol(S)
    cholW <- chol(W)
    halflogdetS <- sum(log(diag(cholS)))
    halflogdetW <- sum(log(diag(cholW)))
    invW <- chol2inv(cholW)
    exptrace <- sum(S*invW)
    lnum <- v*halflogdetS -(v + p + 1)*halflogdetW - 0.5*exptrace
    lpdf <- lnum - ldenom
    return(exp(lpdf))
}


##
## Inverse Gamma
##

## evaluate the inverse gamma density
##
## Kevin Rompala 5/6/2003
## fixed KQ 3/8/2005

#' The Inverse Gamma Distribution
#'
#' Density function and random generation from the inverse gamma distribution.
#'
#' An inverse gamma random variable with shape \eqn{a} and scale \eqn{b}
#' has mean \eqn{\frac{b}{a-1}} (assuming \eqn{a>1}) and variance
#' \eqn{\frac{b^2}{(a-1)^2(a-2)}} (assuming \eqn{a>2}).
#'
#' @aliases dinvgamma rinvgamma
#'
#' @name InvGamma
#'
#' @param x Scalar location to evaluate density.
#'
#' @param n Number of draws from the distribution.
#'
#' @param shape Scalar shape parameter.
#'
#' @param scale Scalar scale parameter (default value one).
#'
#' @return \code{dinvgamma} evaluates the density at \code{x}.
#'
#' \code{rinvgamma} takes \code{n} draws from the inverse Gamma distribution.
#' The parameterization is consistent with the Gamma Distribution in the stats
#' package.
#'
#' @seealso \code{\link[stats]{GammaDist}}
#'
#' @references Andrew Gelman, John B. Carlin, Hal S. Stern, and Donald B.
#' Rubin. 2004. \emph{Bayesian Data Analysis}. 2nd Edition. Boca Raton: Chapman
#' & Hall.
#'
#' @keywords distribution
#'
#' @examples
#'
#' density <- dinvgamma(4.2, 1.1)
#' draws <- rinvgamma(10, 3.2)
#'
NULL

#' @rdname InvGamma
#' @export
"dinvgamma" <-
  function(x, shape, scale = 1) {

    # error checking
    if(shape <= 0 | scale <=0) {
      stop("Shape or scale parameter negative in dinvgamma().\n")
    }

    alpha <- shape
    beta <- scale

    # done on log scale to allow for large alphas and betas
    log.density <- alpha * log(beta) - lgamma(alpha) -
       (alpha + 1) * log(x) - (beta/x)
    return(exp(log.density))
  }

## generate draws from the inverse gamma density (using
## the gamma simulator)
##
## Kevin Rompala 5/6/2003
## fixed KQ 3/8/2005
## shape and rate made explicit 5/25/2010 (KQ)

#' @rdname InvGamma
#' @export
"rinvgamma" <-
  function(n, shape, scale = 1) {
    return(1 / rgamma(n=n, shape=shape, rate=scale))
  }

##
## Dirichlet (Multivariate Beta)
##

# ddirichlet evaluates the density of the Dirichlet at
# vector x given shape parameter vector (or matrix)
# alpha.
#
# note: this code is taken verbatim from the R-package
# "Greg's Miscellaneous Functions" maintained by
# Gregory R. Warnes <Gregory_R_Warnes@groton.pfizer.com>
#
# Kevin Rompala 5/6/2003

#' The Dirichlet Distribution
#'
#' Density function and random generation from the Dirichlet distribution.
#'
#' The Dirichlet distribution is the multidimensional generalization of the
#' beta distribution.
#'
#' @aliases ddirichlet rdirichlet
#'
#' @name Dirichlet
#'
#' @param x A vector containing a single deviate or matrix containing one
#' random deviate per row.
#'
#' @param n Number of random vectors to generate.
#'
#' @param alpha Vector of shape parameters, or matrix of shape parameters
#' corresponding to the number of draw.
#'
#' @return \code{ddirichlet} gives the density. \code{rdirichlet} returns a
#' matrix with \code{n} rows, each containing a single Dirichlet random
#' deviate.
#'
#' @author Code is taken from Greg's Miscellaneous Functions (gregmisc).  His
#' code was based on code posted by Ben Bolker to R-News on 15 Dec 2000.
#'
#' @seealso \code{\link[stats]{Beta}}
#'
#' @keywords distribution
#'
#' @examples
#'
#'   density <- ddirichlet(c(.1,.2,.7), c(1,1,1))
#'   draws <- rdirichlet(20, c(1,1,1) )
#'
NULL

#' @rdname Dirichlet
#' @export
"ddirichlet" <-
  function(x, alpha) {

    dirichlet1 <- function(x, alpha) {
      logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
      s <- sum((alpha-1)*log(x))
      exp(sum(s)-logD)
    }

    # make sure x is a matrix
    if(!is.matrix(x))
      if(is.data.frame(x))
        x <- as.matrix(x)
      else
        x <- t(x)
    if(!is.matrix(alpha))
      alpha <- matrix( alpha, ncol=length(alpha), nrow=nrow(x), byrow=TRUE)

    if( any(dim(x) != dim(alpha)) )
      stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")

    pd <- vector(length=nrow(x))
    for(i in 1:nrow(x))
      pd[i] <- dirichlet1(x[i,],alpha[i,])

    # Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
    pd[ apply( x, 1, function(z) any( z <0 | z > 1)) ] <- 0
    pd[ apply( x, 1, function(z) all.equal(sum( z ),1) !=TRUE) ] <- 0
    return(pd)
  }


# rdirichlet generates n random draws from the Dirichlet at
# vector x given shape parameter vector (or matrix)
# alpha.
#
# note: this code is taken verbatim from the R-package
# "Greg's Miscellaneous Functions" maintained by
# Gregory R. Warnes <Gregory_R_Warnes@groton.pfizer.com>
#
# Kevin Rompala 5/6/2003

#' @rdname Dirichlet
#' @export
"rdirichlet" <-
  function(n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
    sm <- x%*%rep(1,l)
    return(x/as.vector(sm))
  }

##
## Non-Central Hypergeometric
##

# code to evaluate the noncentral hypergeometric density (at a single point
# or at all defined points).
#
# parameters:
#
#    n1, n2 -- number of subjects in group 1 and 2
#
#    Y1, Y2 -- number of subjects with positive outcome [unobserved]
#
#    psi -- odds ratio
#
#    m1 -- sum of observed values of Y1 and Y2 (Y1 + Y2)
#
# output:
#
#   pi -- Pr(Y1 = x | Y1 + Y2 = m1) x=ll,...,uu
#
#   for ll = max(0, m1-n2) and uu = min(n1, m1)
#
# if x is NA, then a matrix is returned, with the first column the possible
# values of x, and the second columns the probabilities.  if x is a scalar,
# the density is evaluated at that point.
#
# ADM on 5/8/2003
#
# note: code adapted from R code published in conjunction with:
#
# Liao, J.G. And Rosen, O. (2001) Fast and Stable Algorithms for Computing and
# Sampling from the Noncentral Hypergeometric Distribution.  The American
# Statistician 55, 366-369.
#

#' The Noncentral Hypergeometric Distribution
#'
#' Evaluates the density at a single point or all points, and generate random
#' draws from the Noncentral Hypergeometric distribution.
#'
#' The Noncentral Hypergeometric is particularly useful for conditional
#' inference for \eqn{(2 \times 2)} tables.  We use the
#' parameterization and algorithms of Liao and Rosen (2001).  The underlying R
#' code is based on their published code.  See their article for details of the
#' parameterization.
#'
#' @aliases rnoncenhypergeom dnoncenhypergeom
#'
#' @name NoncenHypergeom
#'
#' @param x The location to evaluate the density.  If \code{x} is NA, then a
#' matrix is returned with the density evaluated at all possible points.
#'
#' @param n The number of draws to make from the distribution.
#'
#' @param n1 The size of group one.
#'
#' @param n2 The size of group two.
#'
#' @param m1 The observed number of positive outcomes (in both groups).
#'
#' @param psi Odds ratio.
#'
#' @return \code{dnoncenhypergeom} evaluates the density at point \code{x}, or
#' a matrix with the first column containing the possible values of the random
#' variable, and the second column containing the probabilities.
#'
#' \code{rnoncenhypergeom} returns a list of \code{n} random draws from the
#' distribution.
#'
#' @source J. G. Liao and Ori Rosen. 2001. ``Fast and Stable Algorithms for
#' Computing and Sampling From the Noncentral Hypergeometric Distribution."
#' \emph{The American Statistician.} 55: 366-369.
#'
#' @keywords distribution
#'
#' @examples
#'
#'   density <- dnoncenhypergeom(NA, 500, 500, 500, 6.0)
#'   draws <- rnoncenhypergeom(10, 500, 500, 500, 6.0)
NULL

#' @rdname NoncenHypergeom
#' @export
"dnoncenhypergeom" <-
  function (x = NA, n1, n2, m1, psi) {

    ##
    ## AUXILIARY FUNCTIONS
    ##

    mode.compute <- function(n1, n2, m1, psi, ll, uu) {
      a <- psi - 1
      b <- -( (m1+n1+2)*psi + n2-m1 )
      c <- psi*(n1+1)*(m1+1)
      q <- b + sign(b)*sqrt(b*b-4*a*c)
      q <- -q/2

      mode <- trunc(c/q)
      if(uu>=mode && mode>=ll) return(mode)
      else return( trunc(q/a) )
    }

    r.function <- function(n1, n2, m1, psi, i) {
      (n1-i+1)*(m1-i+1)/i/(n2-m1+i)*psi
    }

    ##
    ## MAIN FUNCTION
    ##

    # upper and lower limits for density evaluation
    ll <- max(0, m1-n2)
    uu <- min(n1, m1)

    # check parameters
    if(n1 < 0 | n2 < 0) {
       stop("n1 or n2 negative in dnoncenhypergeom().\n")
    }
    if(m1 < 0 | m1 > (n1 + n2)) {
       stop("m1 out of range in dnoncenhypergeom().\n")
    }
    if(psi <=0) {
       stop("psi [odds ratio] negative in dnoncenhypergeom().\n")
    }
    if(!is.na(x) & (x < ll | x > uu)) {
       stop("x out of bounds in dnoncenhypergeom().\n")
    }
    if(!is.na(x) & length(x) > 1) {
       stop("x neither missing or scalar in dnoncenhypergeom().\n")
    }

    # evaluate density using recursion (from mode)
    mode <- mode.compute(n1, n2, m1, psi, ll, uu)
    pi <- array(1, uu-ll+1)
    shift <- 1-ll

    if(mode<uu) { # note the shift of location
      r1 <- r.function( n1, n2, m1, psi, (mode+1):uu )
      pi[(mode+1 + shift):(uu + shift)] <- cumprod(r1)
    }

    if(mode>ll) {
       r1 <- 1/r.function( n1, n2, m1, psi, mode:(ll+1) )
       pi[(mode-1 + shift):(ll + shift)] <- cumprod(r1)
    }

    pi <- pi/sum(pi)
    if(is.na(x)) return(cbind(ll:uu,pi))
    else return(pi[x + shift])
}

# code to generate random deviates from the noncentral hypergeometric density
#
# parameters:
#
#    n -- the number of draws to make
#
#    n1, n2 -- number of subjects in group 1 and 2
#
#    Y1, Y2 -- number of subjects with positive outcome [unobserved]
#
#    psi -- odds ratio
#
#    m1 -- sum of observed values of Y1 and Y2 (Y1 + Y2)
#
# output:
#
#   output -- a list of length n of random deviates
#
#
# ADM on 5/9/2003
#
# note: code adapted from R code published in conjunction with:
#
# Liao, J.G. And Rosen, O. (2001) Fast and Stable Algorithms for Computing and
# Sampling from the Noncentral Hypergeometric Distribution.  The American
# Statistician 55, 366-369.
#

#' @rdname NoncenHypergeom
#' @export
"rnoncenhypergeom" <-
  function (n, n1, n2, m1, psi) {

    ##
    ## AUXILIARY FUNCTIONS
    ##

    mode.compute <- function(n1, n2, m1, psi, ll, uu) {
      a <- psi - 1
      b <- -( (m1+n1+2)*psi + n2-m1 )
      c <- psi*(n1+1)*(m1+1)
      q <- b + sign(b)*sqrt(b*b-4*a*c)
      q <- -q/2

      mode <- trunc(c/q)
      if(uu>=mode && mode>=ll) return(mode)
      else return( trunc(q/a) )

    }

    sample.low.to.high <- function(lower.end, ran, pi, shift) {
      for(i in lower.end:uu) {
        if(ran <= pi[i+shift]) return(i)
        ran <- ran - pi[i+shift]
        }
    }

    sample.high.to.low <- function(upper.end, ran, pi, shift) {
      for(i in upper.end:ll) {
        if(ran <= pi[i+shift]) return(i)
        ran <- ran - pi[i+shift]
      }
    }

    single.draw <- function(n1, n2, m1, psi, ll, uu, mode, pi) {
      ran <- runif(1)
      shift <- 1-ll
      if(mode==ll) return(sample.low.to.high(ll, ran, pi, shift))
      if(mode==uu) return(sample.high.to.low(uu, ran, pi, shift))
      if(ran < pi[mode+shift]) return(mode)
      ran <- ran - pi[mode+shift]
      lower <- mode - 1
      upper <- mode + 1

      repeat {
        if(pi[upper + shift] >= pi[lower + shift]) {
          if(ran < pi[upper+shift]) return(upper)
          ran <- ran - pi[upper+shift]
          if(upper == uu) return( sample.high.to.low(lower, ran, pi, shift) )
          upper <- upper + 1
        }
        else {
          if(ran < pi[lower+shift]) return(lower)
          ran <- ran - pi[lower+shift]
          if(lower == ll) return( sample.low.to.high(upper, ran, pi, shift) )
          lower <- lower - 1
        }
      }
    }

    ##
    ## MAIN FUNCTION
    ##

    # upper and lower limits for density evaluation
    ll <- max(0, m1-n2)
    uu <- min(n1, m1)

    # check parameters
    if(n1 < 0 | n2 < 0) {
       stop("n1 or n2 negative in rnoncenhypergeom().\n")
    }
    if(m1 < 0 | m1 > (n1 + n2)) {
       stop("m1 out of range in rnoncenhypergeom().\n")
    }
    if(psi <=0) {
       stop("psi [odds ratio] negative in rnoncenhypergeom().\n")
    }


    # get density and other parameters
    mode <- mode.compute(n1, n2, m1, psi, ll, uu)
    pi <- dnoncenhypergeom(NA, n1, n2, m1, psi)[,2]

    output <- array(0,n)
    for(i in 1:n) output[i] <- single.draw(n1, n2, m1, psi, ll, uu, mode, pi)
    return(output)
  }
