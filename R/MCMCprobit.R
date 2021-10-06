##########################################################################
## sample from the posterior distribution of a probit
## model in R using linked C++ code in Scythe
##
## ADM and KQ 5/21/2002
## Modified to meet new developer specification 7/26/2004 KQ
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

#' Markov Chain Monte Carlo for Probit Regression
#'
#' This function generates a sample from the posterior distribution of a probit
#' regression model using the data augmentation approach of Albert and Chib
#' (1993). The user supplies data and priors, and a sample from the posterior
#' distribution is returned as an mcmc object, which can be subsequently
#' analyzed with functions provided in the coda package.
#'
#' \code{MCMCprobit} simulates from the posterior distribution of a probit
#' regression model using data augmentation. The simulation proper is done in
#' compiled C++ code to maximize efficiency.  Please consult the coda
#' documentation for a comprehensive list of functions that can be used to
#' analyze the posterior sample.
#'
#' The model takes the following form:
#'
#' \deqn{y_i \sim \mathcal{B}ernoulli(\pi_i)}
#'
#' Where the inverse link function:
#'
#' \deqn{\pi_i = \Phi(x_i'\beta)}
#'
#' We assume a multivariate Normal prior on \eqn{\beta}:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}
#'
#' See Albert and Chib (1993)
#' for estimation details.
#'
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of Gibbs iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' Gibbs iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number and the betas are printed to the screen every
#' \code{verbose}th iteration.
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
#' of betas.  If this takes a scalar value, then that value will serve as the
#' starting value for all of the betas. The default value of NA will use the
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
#' equivalent to an improper uniform prior on \eqn{\beta}.
#'
#' @param bayes.resid Should latent Bayesian residuals (Albert and Chib, 1995)
#' be returned? Default is FALSE meaning no residuals should be returned.
#' Alternatively, the user can specify an array of integers giving the
#' observation numbers for which latent residuals should be calculated and
#' returned. TRUE will return draws of latent residuals for all observations.
#'
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated, \code{Laplace} in which case the Laplace approximation (see Kass
#' and Raftery, 1995) is used, or \code{Chib95} in which case Chib (1995)
#' method is used.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[stats]{glm}}
#'
#' @references Albert, J. H. and S. Chib. 1993. ``Bayesian Analysis of Binary
#' and Polychotomous Response Data.'' \emph{J. Amer. Statist. Assoc.} 88,
#' 669-679
#'
#' Albert, J. H. and S. Chib. 1995. ``Bayesian Residual Analysis for Binary
#' Response Regression Models.'' \emph{Biometrika.} 82, 747-759.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.
#' <doi: 10.1080/01621459.1995.10476635>
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
#'    data(birthwt)
#'    out1 <- MCMCprobit(low~as.factor(race)+smoke, data=birthwt,
#'    	b0 = 0, B0 = 10, marginal.likelihood="Chib95")
#'    out2 <- MCMCprobit(low~age+as.factor(race), data=birthwt,
#'    	b0 = 0, B0 = 10,  marginal.likelihood="Chib95")
#'    out3 <- MCMCprobit(low~age+as.factor(race)+smoke, data=birthwt,
#'    	b0 = 0, B0 = 10,  marginal.likelihood="Chib95")
#'    BayesFactor(out1, out2, out3)
#'    plot(out3)
#'    summary(out3)
#'    }
#'
"MCMCprobit" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin = 1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, bayes.resid=FALSE,
           marginal.likelihood = c("none", "Laplace", "Chib95"), ...) {

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
    beta.start <- coef.start(beta.start, K, formula,
                             family=binomial(link="probit"), data)
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


    ## residuals setup
    resvec <- NULL
    if (is.logical(bayes.resid) && bayes.resid==TRUE){
      resvec <- matrix(1:length(Y), length(Y), 1)
    }
    else if (!is.logical(bayes.resid)){
      resvec <- matrix(bayes.resid, length(bayes.resid), 1)
      if (min(resvec %in% 1:length(Y)) == 0){
        cat("Elements of bayes.resid are not valid row numbers.\n")
        stop("Check data and call MCMCprobit() again.\n")
      }
    }

    ## y \in {0, 1} error checking
    if (sum(Y!=0 & Y!=1) > 0) {
      cat("Elements of Y equal to something other than 0 or 1.\n")
      stop("Check data and call MCMCprobit() again.\n")
    }

    ## if Chib95 is true
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    posterior <- NULL

    if (is.null(resvec)){
      ## define holder for posterior density sample
      sample <- matrix(data=0, mcmc/thin, dim(X)[2] )

      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="cMCMCprobit",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose), betastart=beta.start,
                       b0=b0, B0=B0, logmarglikeholder.nonconst = as.double(0.0),
                       chib = as.integer(chib))

      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
      }
      ## marginal likelihood calculation if Laplace
      if (marginal.likelihood == "Laplace"){
        theta.start <- beta.start
        optim.out <- optim(theta.start, logpost.probit, method="BFGS",
                           control=list(fnscale=-1),
                           hessian=TRUE, y=Y, X=X, b0=b0, B0=B0)

        theta.tilde <- optim.out$par
        beta.tilde <- theta.tilde[1:K]

        Sigma.tilde <- solve(-1*optim.out$hessian)

        logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
          log(sqrt(det(Sigma.tilde))) +
            logpost.probit(theta.tilde, Y, X, b0, B0)

      }

      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCprobit Posterior Sample",
                                 y=Y, call=cl,
                                 logmarglike=logmarglike)

    }
    else{
      # define holder for posterior density sample
      sample <- matrix(data=0, mcmc/thin, dim(X)[2]+length(resvec) )

      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCprobitres",
                       sample.nonconst=sample, Y=Y, X=X, resvec=resvec,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose), betastart=beta.start,
                       b0=b0, B0=B0,  logmarglikeholder.nonconst= as.double(0.0),
                       chib = as.integer(chib))

      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
      }
      ## marginal likelihood calculation if Laplace
      if (marginal.likelihood == "Laplace"){
        theta.start <- beta.start
        optim.out <- optim(theta.start, logpost.probit, method="BFGS",
                           control=list(fnscale=-1),
                           hessian=TRUE, y=Y, X=X, b0=b0, B0=B0)

        theta.tilde <- optim.out$par
        beta.tilde <- theta.tilde[1:K]

        Sigma.tilde <- solve(-1*optim.out$hessian)

        logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
          log(sqrt(det(Sigma.tilde))) +
            logpost.probit(theta.tilde, Y, X, b0, B0)

      }

      ## put together matrix and build MCMC object to return
      xnames <- c(xnames, paste("epsilonstar", as.character(resvec), sep="") )

      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCprobit Posterior Sample",
                                 y=Y, call=cl, logmarglike=logmarglike)

    }
    return(output)

  }
