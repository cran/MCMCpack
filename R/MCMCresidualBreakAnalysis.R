#########################################################
## residual break analysis
## JHP 07/01/2007
## JHP 03/03/2009
#########################################################

#' Break Analysis of Univariate Time Series using Markov Chain Monte Carlo
#'
#' This function performs a break analysis for univariate time series data
#' using a linear Gaussian changepoint model. The code is written mainly for an
#' internal use in \code{testpanelSubjectBreak}.
#'
#' \code{MCMCresidualBreakAnalysis} simulates from the posterior distribution
#' using standard Gibbs sampling (a multivariate Normal draw for the betas, and
#' an inverse Gamma draw for the conditional error variance).  The simulation
#' proper is done in compiled C++ code to maximize efficiency.  Please consult
#' the coda documentation for a comprehensive list of functions that can be
#' used to analyze the posterior sample.
#'
#' The model takes the following form:
#'
#' \deqn{y_{i} \sim \mathcal{N}(\beta_{m}, \sigma^2_{m}) \;\; m = 1, \ldots, M}
#'
#' We assume standard, semi-conjugate priors:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}
#'
#' And: \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}
#'
#' Where \eqn{\beta} and \eqn{\sigma^{-2}} are
#' assumed \emph{a priori} independent.
#'
#' And:
#'
#' \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}
#'
#' Where \eqn{M} is the number of states.
#'
#' @param resid Univariate time series
#'
#' @param m The number of breaks.
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
#' @param a \eqn{a} is the shape1 beta prior for transition probabilities.
#' By default, the expected duration is computed and corresponding a and b
#' values are assigned. The expected duration is the sample period divided by
#' the number of states.
#'
#' @param b \eqn{b} is the shape2 beta prior for transition probabilities.
#' By default, the expected duration is computed and corresponding a and b
#' values are assigned. The expected duration is the sample period divided by
#' the number of states.
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
#' @param P.start The starting values for the transition matrix.  A user should
#' provide a square matrix with dimension equal to the number of states.  By
#' default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper
#' transition matrix for each raw except the last raw.
#'
#' @param random.perturb If TRUE, randomly sample hidden states whenever
#' regularly sampled hidden states have at least one single observation state.
#' It's one method to avoid overfitting in a non-ergodic hidden Markov models.
#' See Park and Sohn (2017).
#'
#' @param WAIC Compute the Widely Applicable Information Criterion (Watanabe
#' 2010).
#'
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated, and \code{Chib95} in which case the method of Chib (1995) is
#' used.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}},
#' \code{\link[stats]{lm}}
#'
#' @references Jong Hee Park and Yunkyu Sohn. 2017. "Detecting Structural
#' Changes in Network Data: An Application to Changes in Military Alliance
#' Networks, 1816-2012".  Working Paper.
#'
#' Jong Hee Park, 2012. ``Unified Method for Dynamic and Cross-Sectional
#' Heterogeneity: Introducing Hidden Markov Panel Models.'' \emph{American
#' Journal of Political Science}.56: 1040-1054. <doi: 10.1111/j.1540-5907.2012.00590.x>
#'
#' Sumio Watanabe. 2010. "Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory"
#' \emph{Journal of Machine Learning Research}. 11: 3571-3594.
#'
#' Siddhartha Chib. 1995. "Marginal Likelihood from the Gibbs Output."
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.
#' <doi: 10.1016/S0304-4076(97)00115-2>
#'
#' Siddhartha Chib. 1998. "Estimation and comparison of multiple change-point
#' models."  \emph{Journal of Econometrics}. 86: 221-241.
#'  <doi: 10.1080/01621459.1995.10476635>
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' line   <- list(X = c(-2,-1,0,1,2), Y = c(1,3,3,3,5))
#' ols <- lm(Y~X)
#' residual <-   rstandard(ols)
#' posterior  <- MCMCresidualBreakAnalysis(residual, m = 1, data=line, mcmc=1000, verbose=200)
#' plotState(posterior)
#' summary(posterior)
#' }
#'
"MCMCresidualBreakAnalysis"<-
  function(resid, m = 1,
           b0 = 0, B0 = 0.001, c0 = 0.1, d0 = 0.1, a = NULL, b = NULL,
           mcmc = 1000, burnin = 1000,  thin = 1, verbose = 0,
           seed = NA, beta.start = NA, P.start = NA, random.perturb = FALSE,
           WAIC = FALSE, marginal.likelihood = c("none", "Chib95"), ...){

    ## form response and model matrices
    y <- as.matrix(resid, ,1)
    n <- nrow(y)
    ns <- m + 1 # number of states

    ## check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    cl <- match.call()
    nstore <- mcmc/thin

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    if(!is.na(seed)) set.seed(seed)

   ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)

    ## following MCMCregress, set marginalrun as binary
    logmarglike <- loglike <- NULL
    marginalrun <- ifelse(marginal.likelihood == "Chib95", 1, 0)
    sos <- ifelse(random.perturb, 1, 0)

    ## initial values
    if(m == 0){
      output <- MCMCregress(y~1, mcmc=mcmc, burnin=burnin, verbose=verbose, thin=thin,
                            b0=b0, B0=B0, c0=c0, d0=d0, seed = seed,
                            beta.start=beta.start, marginal.likelihood = marginal.likelihood)
      attr(output, "y") <- y
    }
    else{
      ## prior
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)
      Pstart  <-  check.P(P.start, m, a=a, b=b)
      betastart  <- rep(mean(y), ns)
      Sigmastart <- rep(var(y), ns)
      statestart <- sort(sample(1:ns, n, replace=T))

      ## call C++ code to draw sample
      posterior <- .C("cMCMCresidualBreakAnalysis",
                      betaout = as.double(rep(0.0, nstore*ns)),
                      Sigmaout = as.double(rep(0.0, nstore*ns)),
                      ## Pout = as.double(rep(0.0, nstore*ns*ns)),
                      psout = as.double(rep(0.0, n*ns)),
                      sout = as.double(rep(0.0, nstore*n)),
                      yloglike = as.double(rep(0.0, nstore*n)),

                      Ydata = as.double(y),
                      Yrow = as.integer(nrow(y)),
                      Ycol = as.integer(ncol(y)),

                      m = as.integer(m),
                      burnin = as.integer(burnin),
                      mcmc = as.integer(mcmc),
                      thin = as.integer(thin),
                      verbose = as.integer(verbose),

                      lecuyer=as.integer(lecuyer),
                      seedarray=as.integer(seed.array),
                      lecuyerstream=as.integer(lecuyer.stream),

                      betastart = as.double(betastart),
                      Sigmastart = as.double(Sigmastart),
                      Pstart = as.double(Pstart),
                      statestart = as.integer(statestart),

                      a = as.double(a),
                      b = as.double(b),
                      b0data = as.double(b0),
                      B0data = as.double(B0),
                      c0 = as.double(c0),
                      d0 = as.double(d0),
                      A0data = as.double(A0),
                      logmarglikeholder = as.double(0.0),
                      loglikeholder = as.double(0.0),
                      marginalrun = as.integer(marginalrun),
                      sos = as.integer(sos))

      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
        loglike <- posterior$loglikeholder
      }
      Waic.out <- NA
      if(WAIC){
          Y.loglike.mat <- matrix(posterior$yloglike, nrow = nstore, ncol = n)
          Waic.out <- waic(Y.loglike.mat)$total
          rm(Y.loglike.mat)
          cat("    Waic: ", Waic.out[1], "\n")
      }

      ## pull together matrix and build MCMC object to return
      beta.holder <- matrix(posterior$betaout, nstore, ns)
      Sigma.holder <- matrix(posterior$Sigmaout, nstore, ns)
      ## P.holder    <- matrix(posterior$Pout, nstore, )
      ps.holder   <- matrix(posterior$psout, n, )
      s.holder    <- matrix(posterior$sout, nstore, )

      output1 <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output1)  <- sapply(c(1:ns),
                                   function(i){
                                     paste("mu_regime", i, sep = "")
                                   })
      output2 <- mcmc(data=Sigma.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output2)  <- sapply(c(1:ns),
                                   function(i){
                                   paste("sigma2_regime", i, sep = "")
                                 })
      output <- as.mcmc(cbind(output1, output2))

      attr(output, "title") <- "MCMCresidualBreakAnalysis Posterior Sample"
      attr(output, "formula") <- formula
      attr(output, "y")       <- y
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      ## attr(output, "P.store") <- P.holder
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "s.store") <- s.holder
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglike") <- loglike
    }
    return(output)

}## end of MCMC function
