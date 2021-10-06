## sample from the posterior distribution
## of a binary model with multiple changepoints
## using linked C++ code in Scythe
##
## Written by JHP 07/01/2007
## Revised by JHP 07/16/2009

#' Markov Chain Monte Carlo for a Binary Multiple Changepoint Model
#'
#' This function generates a sample from the posterior distribution of
#' a binary model with multiple changepoints. The function uses the
#' Markov chain Monte Carlo method of Chib (1998).  The user supplies
#' data and priors, and a sample from the posterior distribution is
#' returned as an mcmc object, which can be subsequently analyzed with
#' functions provided in the coda package.
#'
#' \code{MCMCbinaryChange} simulates from the posterior distribution
#' of a binary model with multiple changepoints.
#'
#' The model takes the following form: \deqn{Y_t \sim
#' \mathcal{B}ernoulli(\phi_i),\;\; i = 1, \ldots, k} Where \eqn{k}
#' is the number of states.
#'
#' We assume Beta priors for \eqn{\phi_{i}} and for transition
#' probabilities: \deqn{\phi_i \sim \mathcal{B}eta(c_0, d_0)} And:
#' \deqn{p_{mm} \sim \mathcal{B}eta{a}{b},\;\; m = 1, \ldots, k} Where
#' \eqn{M} is the number of states.
#'
#' @param data The data.
#'
#' @param m The number of changepoints.
#'
#' @param c0 \eqn{c_0} is the shape1 parameter for Beta prior on
#'   \eqn{\phi} (the mean).
#'
#' @param d0 \eqn{d_0} is the shape2 parameter for Beta prior on
#'   \eqn{\phi} (the mean).
#'
#' @param a \eqn{a} is the shape1 beta prior for transition
#'   probabilities.  By default, the expected duration is computed and
#'   corresponding a and b values are assigned. The expected duration
#'   is the sample period divided by the number of states.
#'
#' @param b \eqn{b} is the shape2 beta prior for transition
#'   probabilities.  By default, the expected duration is computed and
#'   corresponding a and b values are assigned. The expected duration
#'   is the sample period divided by the number of states.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burn-in.
#'
#' @param thin The thinning interval used in the simulation.  The
#'   number of MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the
#'   progress of the sampler is printed to the screen.  If
#'   \code{verbose} is greater than 0, the iteration number and the
#'   posterior density samples are printed to the screen every
#'   \code{verbose}th iteration.
#'
#' @param seed The seed for the random number generator.  If NA,
#'   current R system seed is used.
#'
#' @param phi.start The starting values for the mean. The default
#'   value of NA will use draws from the Uniform distribution.
#'
#' @param P.start The starting values for the transition matrix. A
#'   user should provide a square matrix with dimension equal to the
#'   number of states. By default, draws from the \code{Beta(0.9,
#'   0.1)} are used to construct a proper transition matrix for each
#'   raw except the last raw.
#'
#' @param marginal.likelihood How should the marginal likelihood be
#'   calculated?  Options are: \code{none} in which case the marginal
#'   likelihood will not be calculated, and \code{Chib95} in which
#'   case the method of Chib (1995) is used.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This
#'   object can be summarized by functions provided by the coda
#'   package.  The object contains an attribute \code{prob.state}
#'   storage matrix that contains the probability of \eqn{state_i} for
#'   each period, and the log-marginal likelihood of the model
#'   (\code{logmarglike}).
#'
#' @export
#'
#' @seealso \code{\link{MCMCpoissonChange}},\code{\link{plotState}},
#'   \code{\link{plotChangepoint}}
#'
#' @references Jong Hee Park. 2011. ``Changepoint Analysis of Binary
#'   and Ordinal Probit Models: An Application to Bank Rate Policy
#'   Under the Interwar Gold Standard."
#'   \emph{Political Analysis}. 19: 188-204.
#'   <doi:10.1093/pan/mpr007>
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#' ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of
#' Statistical Software}.  42(9): 1-21.
#' \doi{10.18637/jss.v042.i09}.
#'
#' Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs
#' Output.''  \emph{Journal of the American Statistical
#' Association}. 90: 1313-1321. <doi: 10.1080/01621459.1995.10476635>
#'
#' @keywords models
#'
#' @examples
#'
#'     \dontrun{
#'     set.seed(19173)
#'     true.phi<- c(0.5, 0.8, 0.4)
#'
#'     ## two breaks at c(80, 180)
#'     y1 <- rbinom(80, 1,  true.phi[1])
#'     y2 <- rbinom(100, 1, true.phi[2])
#'     y3 <- rbinom(120, 1, true.phi[3])
#'     y  <- as.ts(c(y1, y2, y3))
#'
#'     model0 <- MCMCbinaryChange(y, m=0, c0=2, d0=2, mcmc=100, burnin=100, verbose=50,
#'     	      marginal.likelihood = "Chib95")
#'     model1 <- MCMCbinaryChange(y, m=1, c0=2, d0=2, mcmc=100, burnin=100, verbose=50,
#'     	      marginal.likelihood = "Chib95")
#'     model2 <- MCMCbinaryChange(y, m=2, c0=2, d0=2, mcmc=100, burnin=100, verbose=50,
#'     	      marginal.likelihood = "Chib95")
#'     model3 <- MCMCbinaryChange(y, m=3, c0=2, d0=2, mcmc=100, burnin=100, verbose=50,
#'     	      marginal.likelihood = "Chib95")
#'     model4 <- MCMCbinaryChange(y, m=4, c0=2, d0=2, mcmc=100, burnin=100, verbose=50,
#'     	      marginal.likelihood = "Chib95")
#'     model5 <- MCMCbinaryChange(y, m=5, c0=2, d0=2, mcmc=100, burnin=100, verbose=50,
#'     	      marginal.likelihood = "Chib95")
#'
#'     print(BayesFactor(model0, model1, model2, model3, model4, model5))
#'
#'     ## plot two plots in one screen
#'     par(mfrow=c(attr(model2, "m") + 1, 1), mai=c(0.4, 0.6, 0.3, 0.05))
#'     plotState(model2, legend.control = c(1, 0.6))
#'     plotChangepoint(model2, verbose = TRUE, ylab="Density", start=1, overlay=TRUE)
#'
#'     }
#'
"MCMCbinaryChange"<-
  function(data,  m = 1, c0 = 1,  d0 = 1,  a = NULL, b = NULL,
           burnin = 10000, mcmc = 10000, thin = 1, verbose = 0,
           seed = NA, phi.start = NA, P.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...) {

    ## check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    cl <- match.call()

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    if(!is.na(seed)) set.seed(seed)

    ## sample size
    y <- as.matrix(data)
    n <- nrow(y)
    ns <- m+1


    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)

    ## following MCMCregress, set chib as binary
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    nstore <- mcmc/thin
    if (m == 0){
      b0 <- c0/(c0 + d0)
      B0 <- c0*d0/(c0 + d0)^2*(c0 + d0 + 1)
      output <- MCMCprobit(y~1, burnin = burnin, mcmc = mcmc,
                           thin = thin, verbose = verbose, b0 = b0, B0 = B0,
                           marginal.likelihood = marginal.likelihood)
      attr(output, "y") <- y
    }
    else {
      ## prior for transition matrix
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)
      Pstart <- check.P(P.start, m=m, a=a, b=b)
      phistart <- check.theta(phi.start, ns, y, min=0, max=1)


      ## call C++ code to draw sample
      posterior <- .C("cMCMCbinaryChange",
                      phiout = as.double(rep(0.0, nstore*ns)),
                      Pout = as.double(rep(0.0, nstore*ns*ns)),
                      psout = as.double(rep(0.0, n*ns)),
                      sout = as.double(rep(0.0, nstore*n)),
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
                      phistart = as.double(phistart),
                      Pstart = as.double(Pstart),
                      a = as.double(a),
                      b = as.double(b),
                      c0 = as.double(c0),
                      d0 = as.double(d0),
                      A0data = as.double(A0),
                      logmarglikeholder = as.double(0.0),
                      chib = as.integer(chib))

      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
      }

      ## pull together matrix and build MCMC object to return
      phi.holder <- matrix(posterior$phiout, nstore, )
      P.holder    <- matrix(posterior$Pout,  nstore, )
      s.holder    <- matrix(posterior$sout,  nstore, )
      ps.holder   <- matrix(posterior$psout, n, )

      output <- mcmc(data=phi.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output)  <- paste("phi.", 1:ns, sep = "")
      attr(output,"title") <- "MCMCbinaryChange Posterior Sample"
      attr(output, "y")       <- y
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      attr(output, "logmarglike") <- logmarglike
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "s.store") <- s.holder
    }
    return(output)
  }## end of MCMC function
