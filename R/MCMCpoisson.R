##########################################################################
## sample from the posterior distribution of a Poisson regression
## model in R using linked C++ code in Scythe
##
## ADM 1/24/2003
## KQ 3/17/2003 [bug fix]
## Modified to meet new developer specification 7/15/2004 KQ
## Modified for new Scythe and rngs 7/26/2004 KQ
## Modified to handle marginal likelihood calculation 1/27/2006 KQ
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for Poisson Regression
#'
#' This function generates a sample from the posterior distribution of a
#' Poisson regression model using a random walk Metropolis algorithm. The user
#' supplies data and priors, and a sample from the posterior distribution is
#' returned as an mcmc object, which can be subsequently analyzed with
#' functions provided in the coda package.
#'
#' \code{MCMCpoisson} simulates from the posterior distribution of a Poisson
#' regression model using a random walk Metropolis algorithm. The simulation
#' proper is done in compiled C++ code to maximize efficiency.  Please consult
#' the coda documentation for a comprehensive list of functions that can be
#' used to analyze the posterior sample.
#'
#' The model takes the following form:
#'
#' \deqn{y_i \sim \mathcal{P}oisson(\mu_i)}
#'
#' Where the inverse link
#' function:
#'
#' \deqn{\mu_i = \exp(x_i'\beta)}
#'
#' We assume a multivariate Normal prior on \eqn{\beta}:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}
#'
#' The Metropois proposal distribution is centered at the current value of
#' \eqn{\theta} and has variance-covariance \eqn{V = T (B_0 + C^{-1})^{-1} T }
#' where \eqn{T} is a the diagonal positive definite matrix formed from the
#' \code{tune}, \eqn{B_0} is the prior precision, and \eqn{C} is the
#' large sample variance-covariance matrix of the MLEs. This last calculation
#' is done via an initial call to \code{glm}.
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of Metropolis iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' mcmc iterations must be divisible by this value.
#'
#' @param tune Metropolis tuning parameter. Can be either a positive scalar or
#' a \eqn{k}-vector, where \eqn{k} is the length of
#' \eqn{\beta}.Make sure that the acceptance rate is satisfactory
#' (typically between 0.20 and 0.5) before using the posterior sample for
#' inference.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the current beta vector, and the Metropolis acceptance
#' rate are printed to the screen every \code{verbose}th iteration.
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
#' @param beta.start The starting value for the \eqn{\beta} vector.  This
#' can either be a scalar or a column vector with dimension equal to the number
#' of betas. If this takes a scalar value, then that value will serve as the
#' starting value for all of the betas.  The default value of NA will use the
#' maximum likelihood estimate of \eqn{\beta} as the starting value.
#'
#' @param b0 The prior mean of \eqn{\beta}.  This can either be a scalar
#' or a column vector with dimension equal to the number of betas. If this
#' takes a scalar value, then that value will serve as the prior mean for all
#' of the betas.
#'
#' @param B0 The prior precision of \eqn{\beta}.  This can either be a
#' scalar or a square matrix with dimensions equal to the number of betas.  If
#' this takes a scalar value, then that value times an identity matrix serves
#' as the prior precision of \eqn{\beta}. Default value of 0 is
#' equivalent to an improper uniform prior for beta.
#'
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated or \code{Laplace} in which case the Laplace approximation (see
#' Kass and Raftery, 1995) is used.
#'
#' @param ... further arguments to be passed.
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[stats]{glm}}
#'
#' @references Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#' ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical
#' Software}. 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
#'    counts <- c(18,17,15,20,10,20,25,13,12)
#'    outcome <- gl(3,1,9)
#'    treatment <- gl(3,3)
#'    posterior <- MCMCpoisson(counts ~ outcome + treatment)
#'    plot(posterior)
#'    summary(posterior)
#'    }
#'
"MCMCpoisson" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin=1, tune=1.1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0,
           marginal.likelihood = c("none", "Laplace"),...) {

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
    beta.start <- coef.start(beta.start, K, formula, family=poisson, data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]


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


    ## form the tuning parameter
    tune <- vector.tune(tune, K)
    V <- vcov(glm(formula=formula, data=data, family=poisson))

    ## test y non-negative
    if (sum(Y < 0) > 0) {
      cat("\n Elements of Y negative. ")
      stop("\n Check data and call MCMCpoisson() again. \n")
    }

    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, dim(X)[2] )



    ## marginal likelihood calculation if Laplace
    if (marginal.likelihood == "Laplace"){
      theta.start <- beta.start
      optim.out <- optim(theta.start, logpost.poisson, method="BFGS",
                         control=list(fnscale=-1),
                         hessian=TRUE, y=Y, X=X, b0=b0, B0=B0)

      theta.tilde <- optim.out$par
      beta.tilde <- theta.tilde[1:K]

      Sigma.tilde <- solve(-1*optim.out$hessian)

      logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
        log(sqrt(det(Sigma.tilde))) +
          logpost.poisson(theta.tilde, Y, X, b0, B0)

    }
    posterior <- NULL

    ## call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="cMCMCpoisson",
                     sample.nonconst=sample, Y=Y, X=X,
                     burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     tune=tune, lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, V=V)

    ## put together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior, names=xnames,
                               title="MCMCpoisson Posterior Sample",
                               y=Y, call=cl, logmarglike=logmarglike)
    return(output)

  }
