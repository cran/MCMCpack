##########################################################################
## MCMCregress.R samples from the posterior distribution of a Gaussian
## linear regression model in R using linked C++ code in Scythe
##
## Original written by ADM and KQ 5/21/2002
## Updated with helper functions ADM 5/28/2004
## Modified to meet new developer specification 6/18/2004 KQ
## Modified for new Scythe and rngs 7/22/2004 ADM
## Modified to handle marginal likelihood calculation 1/26/2006 KQ
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for Gaussian Linear Regression
#'
#' This function generates a sample from the posterior distribution of a linear
#' regression model with Gaussian errors using Gibbs sampling (with a
#' multivariate Gaussian prior on the beta vector, and an inverse Gamma prior
#' on the conditional error variance).  The user supplies data and priors, and
#' a sample from the posterior distribution is returned as an mcmc object,
#' which can be subsequently analyzed with functions provided in the coda
#' package.
#'
#' \code{MCMCregress} simulates from the posterior distribution using standard
#' Gibbs sampling (a multivariate Normal draw for the betas, and an inverse
#' Gamma draw for the conditional error variance).  The simulation proper is
#' done in compiled C++ code to maximize efficiency.  Please consult the coda
#' documentation for a comprehensive list of functions that can be used to
#' analyze the posterior sample.
#'
#' The model takes the following form:
#'
#' \deqn{y_i = x_i ' \beta + \varepsilon_{i}}
#'
#' Where the errors are assumed to be Gaussian:
#'
#' \deqn{\varepsilon_{i} \sim \mathcal{N}(0, \sigma^2)}
#'
#' We assume standard, semi-conjugate priors:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}
#'
#' And:
#'
#' \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}
#'
#' Where
#' \eqn{\beta} and \eqn{\sigma^{-2}} are assumed \emph{a
#' priori} independent.  Note that only starting values for \eqn{\beta}
#' are allowed because simulation is done using Gibbs sampling with the
#' conditional error variance as the first block in the sampler.
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burnin.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the \eqn{\beta} vector, and the error variance are
#' printed to the screen every \code{verbose}th iteration.
#'
#' @param seed The seed for the random number generator.  If NA, the Mersenne
#' Twister generator is used with default seed 12345; if an integer is passed
#' it is used to seed the Mersenne twister.  The user can also pass a list of
#' length two to use the L'Ecuyer random number generator, which is suitable
#' for parallel computation.  The first element of the list is the L'Ecuyer
#' seed, which is a vector of length six or NA (if NA a default seed of
#' \code{rep(12345,6)} is used).  The second element of list is a positive
#' substream number. See the MCMCpack specification for more details.
#'
#' @param beta.start The starting values for the \eqn{\beta} vector.
#' This can either be a scalar or a column vector with dimension equal to the
#' number of betas.  The default value of of NA will use the OLS estimate of
#' \eqn{\beta} as the starting value.  If this is a scalar, that value
#' will serve as the starting value mean for all of the betas.
#'
#' @param b0 The prior mean of \eqn{\beta}.  This can either be a scalar
#' or a column vector with dimension equal to the number of betas. If this
#' takes a scalar value, then that value will serve as the prior mean for all
#' of the betas.
#'
#' @param B0 The prior precision of \eqn{\beta}.  This can either be a
#' scalar or a square matrix with dimensions equal to the number of betas.  If
#' this takes a scalar value, then that value times an identity matrix serves
#' as the prior precision of beta. Default value of 0 is equivalent to an
#' improper uniform prior for beta.
#'
#' @param c0 \eqn{c_0/2} is the shape parameter for the inverse Gamma
#' prior on \eqn{\sigma^2} (the variance of the disturbances). The
#' amount of information in the inverse Gamma prior is something like that from
#' \eqn{c_0} pseudo-observations.
#'
#' @param d0 \eqn{d_0/2} is the scale parameter for the inverse Gamma
#' prior on \eqn{\sigma^2} (the variance of the disturbances). In
#' constructing the inverse Gamma prior, \eqn{d_0} acts like the sum of
#' squared errors from the \eqn{c_0} pseudo-observations.
#'
#' @param sigma.mu The mean of the inverse Gamma prior on
#' \eqn{\sigma^2}.  \eqn{sigma.mu} and
#' \eqn{sigma.var} allow users to choose the inverse Gamma prior by
#' choosing its mean and variance.
#'
#' @param sigma.var The variacne of the inverse Gamma prior on
#' \eqn{\sigma^2}.  \eqn{sigma.mu} and
#' \eqn{sigma.var} allow users to choose the inverse Gamma prior by
#' choosing its mean and variance.
#'
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated, \code{Laplace} in which case the Laplace approximation (see Kass
#' and Raftery, 1995) is used, and \code{Chib95} in which case the method of
#' Chib (1995) is used.
#'
#' @param ... further arguments to be passed.
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}},
#' \code{\link[stats]{lm}}
#'
#' @references Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#' ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical
#' Software}. 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.
#'
#' Robert E. Kass and Adrian E. Raftery. 1995. ``Bayes Factors.'' \emph{Journal
#' of the American Statistical Association}. 90: 773-795.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.wustl.edu.s3-website-us-east-1.amazonaws.com/}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' line   <- list(X = c(-2,-1,0,1,2), Y = c(1,3,3,3,5))
#' posterior  <- MCMCregress(Y~X, b0=0, B0 = 0.1,
#' 	      sigma.mu = 5, sigma.var = 25, data=line, verbose=1000)
#' plot(posterior)
#' raftery.diag(posterior)
#' summary(posterior)
#' }
#'
"MCMCregress" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, sigma.mu = NA, sigma.var = NA,
           marginal.likelihood = c("none", "Laplace", "Chib95"),
           ...) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    cl <- match.call()

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    K <- ncol(X)  # number of covariates

    ## starting values and priors
    beta.start <- coef_start(beta.start, K, formula, family=gaussian, data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    if (is.na(sigma.mu)|is.na(sigma.var)) {
      check.ig.prior(c0, d0)
    }
    else {
      c0 <- 4 + 2 *(sigma.mu^2/sigma.var)
      d0 <- 2*sigma.mu *(c0/2 - 1)
    }

    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    B0.eigenvalues <- eigen(B0)$values
    if (min(B0.eigenvalues) < 0){
      stop("B0 is not positive semi-definite.\nPlease respecify and call again.\n")
    }
    if (isTRUE(all.equal(min(B0.eigenvalues), 0))){
      if (marginal.likelihood != "none"){
        warning("Cannot calculate marginal likelihood with improper prior\n")
        marginal.likelihood <- "none"
      }
    }
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }


    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, K+1)
    posterior <- NULL

    ## call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="cMCMCregress",
                     sample.nonconst=sample, Y=Y, X=X,
                     burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, c0=as.double(c0), d0=as.double(d0),
                     logmarglikeholder.nonconst=as.double(0.0),
                     chib=as.integer(chib))

    ## get marginal likelihood if Chib95
    if (marginal.likelihood == "Chib95"){
      logmarglike <- posterior$logmarglikeholder
    }
    ## marginal likelihood calculation if Laplace
    if (marginal.likelihood == "Laplace"){
      theta.start <- c(beta.start, log(0.5*var(Y)))
      optim.out <- optim(theta.start, logpost.regress, method="BFGS",
                         control=list(fnscale=-1),
                         hessian=TRUE, y=Y, X=X, b0=b0, B0=B0, c0=c0, d0=d0)

      theta.tilde <- optim.out$par
      beta.tilde <- theta.tilde[1:K]
      sigma2.tilde <- exp(theta.tilde[K+1])

      Sigma.tilde <- solve(-1*optim.out$hessian)

      logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
        log(sqrt(det(Sigma.tilde))) +
          logpost.regress(theta.tilde, Y, X, b0, B0, c0, d0)

    }

    ## pull together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior,
                               names=c(xnames, "sigma2"),
                               title="MCMCregress Posterior Sample",
                               y=Y, call=cl,
                               logmarglike=logmarglike
                               )
    return(output)
  }
