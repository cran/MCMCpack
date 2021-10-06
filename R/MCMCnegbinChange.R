################################
## Negative Binomial Changepoint Model
##
##
## 05/20/2017 Matthew Blackwell
################################

#' Markov Chain Monte Carlo for Negative Binomial Regression
#' Changepoint Model
#'
#' This function generates a sample from the posterior distribution of
#' a Negative Binomial regression model with multiple changepoints.
#' For the changepoints, the sampler uses the Markov Chain Monte Carlo
#' method of Chib (1998). The user supplies data and priors, and a
#' sample from the posterior distribution is returned as an mcmc
#' object, which can be subsequently analyzed with functions provided
#' in the coda package.
#'
#' \code{MCMCnegbinChange}simulates from the posterior distribution of a
#' Negative Binomial regression model with multiple changepoints using the methods of
#' Chib (1998) and Fruehwirth-Schnatter et al (2009).  The details of the
#' model are discussed in Blackwell (2017).
#'
#'  The model takes the following form:
#'
#' \deqn{y_t \sim \mathcal{P}oisson(\nu_t\mu_t)}
#'
#' \deqn{\mu_t = x_t ' \beta_m,\;\; m = 1, \ldots, M}
#'
#' \deqn{\nu_t \sim \mathcal{G}amma(\rho_m, \rho_m)}
#'
#'
#' Where
#' \eqn{M} is the number of states and \eqn{\beta_m} and \eqn{\rho_m}
#' are parameters when a state is \eqn{m} at \eqn{t}.
#'
#' We assume Gaussian distribution for prior of \eqn{\beta}:
#'
#' \deqn{\beta_m \sim \mathcal{N}(b_0,B_0^{-1}),\;\; m = 1, \ldots, M}
#'
#' And:
#'
#' \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M} Where \eqn{M} is the number of states.
#'
#' The overdispersion parameters have a prior with the following form:
#'
#' \deqn{f(\rho_m|e,f,g) \propto \rho^{e-1}(\rho + g)^{-(e+f)}}
#'
#' The model is simulated via blocked Gibbs conditonal on the states.
#' The \eqn{\beta} being simulated via the auxiliary mixture sampling
#' method of Fuerhwirth-Schanetter et al. (2009). The \eqn{\rho} is
#' updated via slice sampling. The \eqn{\nu_i} are updated their
#' (conjugate) full conditional, which is also Gamma. The states are
#' updated as in Chib (1998)
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param m The number of changepoints.
#'
#' @param fixed.m A logical indicator for whether or not the number of
#'   changepoints in the sampler should be exactly equal to \code{m}
#'   or if that is simply an upper bound. Setting \code{fixed.m} to
#'   \code{FALSE} is equivalent to assuming a weak-limit approximation
#'   to a Dirichlet process mixture.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of Metropolis iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' mcmc iterations must be divisible by this value.
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
#' @param e The hyperprior for the distribution \eqn{\rho} See details.
#'
#' @param f The hyperprior for the distribution \eqn{\rho}. See details.
#'
#' @param g The hyperprior for the distribution \eqn{\rho}. See details.
#'
#' @param beta.start The starting value for the \eqn{\beta} vector.  This
#' can either be a scalar or a column vector with dimension equal to the number
#' of betas. If this takes a scalar value, then that value will serve as the
#' starting value for all of the betas.  The default value of NA will use the
#' maximum likelihood estimate of \eqn{\beta} as the starting value
#'   for all regimes.
#'
#' @param P.start The starting values for the transition matrix. A user should
#' provide a square matrix with dimension equal to the number of states. By
#' default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper
#' transition matrix for each raw except the last raw.
#'
#' @param rho.start The starting value for the \eqn{\rho} variable.
#'   This can either be a scalar or a column vector with dimension
#'   equal to the number of regimes. If the value is scalar, it will
#'   be used for all regimes. The default value is a vector of ones.
#'
#' @param nu.start The starting values for the random effect,
#'     \eqn{\nu}. The default value is a vector of ones.
#'
#' @param rho.step Tuning parameter for the slice sampling approach to
#'     sampling \eqn{rho}. Determines the size of the step-out used to
#'     find the correct slice to draw from. Lower values are more
#'     accurate, but will take longer (up to a fixed searching limit).
#'     Default is 0.1.
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
#' @seealso \code{\link{MCMCpoissonChange}}, \code{\link{plotState}},
#' \code{\link{plotChangepoint}}
#'
#' @references Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#' ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical
#' Software}. 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Sylvia Fruehwirth-Schnatter, Rudolf Fruehwirth, Leonhard Held, and
#'     Havard Rue. 2009. ``Improved auxiliary mixture sampling for
#'     hierarchical models of non-Gaussian data'', \emph{Statistics
#'     and Computing} 19(4): 479-492.
#'     <doi:10.1007/s11222-008-9109-4>
#'
#' Matthew Blackwell. 2017. ``Game Changers: Detecting Shifts in
#'   Overdispersed Count Data,'' \emph{Political Analysis}
#'   26(2), 230-239. <doi:10.1017/pan.2017.42>
#'
#' @keywords models
#'
#' @examples
#'
#'  \dontrun{
#'    n <- 150
#'    reg <- 3
#'    true.s <- gl(reg, n/reg, n)
#'    rho.true <- c(1.5, 0.5, 3)
#'    b0.true <- c(1, 3, 1)
#'    b1.true <- c(1, -2, 2)
#'    x1 <- runif(n, 0, 2)
#'    nu.true <- rgamma(n, rho.true[true.s], rho.true[true.s])
#'    mu <- nu.true * exp(b0.true[true.s] + x1 * b1.true[true.s])
#'    y <- rpois(n, mu)
#'
#'    posterior <- MCMCnegbinChange(y ~ x1, m = 2, verbose = 1000,
#'                           marginal.likelihood = "Chib95",
#'                           e = 2, f = 2, g = 10,
#'                           b0 = rep(0, 2), B0 = (1/9) * diag(2),
#'                           rho.step = rep(0.75, times = 3),
#'                           seed = list(NA, 2))
#'
#'    par(mfrow=c(attr(posterior, "m") + 1, 1), mai=c(0.4, 0.6, 0.3, 0.05))
#'    plotState(posterior, legend.control = c(1, 0.6))
#'    plotChangepoint(posterior, verbose = TRUE, ylab="Density",
#'   start=1, overlay=TRUE)
#'
#'
#'    open.ended <- MCMCnegbinChange(y ~ x1, m = 10, verbose = 1000,
#'                           fixed.m = FALSE, mcmc = 2000, burnin = 2000,
#'                           e = 2, f = 2, g = 10,
#'                           a = 100, b = 1,
#'                           b0 = rep(0, 2), B0 = (1/9) * diag(2),
#'                           rho.step = 0.75,
#'                           seed = list(NA, 2))
#'
#'    plotState(open.ended, legend.control = c(1, 0.6))
#'    }
#'
"MCMCnegbinChange"<-
  function(formula, data = parent.frame(), m = 1, fixed.m = TRUE,
           b0 = 0, B0 = 1, a = NULL, b = NULL,
           e = 2, f = 2, g = 10,
           burnin = 1000, mcmc = 1000, thin = 1, verbose = 0,
           seed = NA, beta.start = NA, P.start = NA,
           rho.start = NA, rho.step, nu.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...) {

    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    k <- ncol(X)
    n <- length(y)
    n.arrival<-  y + 1
    NT      <-  max(n.arrival)
    tot.comp <-  n + sum(y)
    ns      <- m + 1

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
    if(!any(is.na(seed)) & (length(seed) == 1)) set.seed(seed)


    mvn.prior <- form.mvn.prior(b0, B0, k)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }
    if (m == 0){
      output <- MCMCnegbin(formula, data = data, burnin = burnin, mcmc = mcmc,
                           thin = thin, verbose = verbose, b0 = b0, B0 = B0,
                           e = e, f = f, g = g, rho.step = rho.step,
                           marginal.likelihood = marginal.likelihood)
      return(output)
    }

    if (!fixed.m & chib == 1) {
      chib <- 0
      warning("Marginal likelihood not implemented for open-ended changepoint model.")
    }

    ## prior
    if (m == 0) {
      A0 <- 0
    } else {
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)
    }
    ## get initial values of tau from observed y
    if (m > 0) {
      Pstart <- check.P(P.start, m, a = a, b = b)
    } else {
      Pstart <- 1
    }
    betastart  <- beta.negbin.start(beta.start, ns, k, formula, data)
    if (is.na(rho.start[1])) {
      rho.start <- rep(1, times = m+1)
    }
    if (!(length(rho.step) %in% c(1, m+1))) {
      stop("rho.step has the wrong dimensions")
    }
    if (length(rho.step) == 1) {
      rho.step <- rep(rho.step, times = m+1)
    }
    if (is.na(nu.start[1])) {
      nu.start <- rep(1, times = n) #rgamma(n, rho.start[1], rho.start[1])
    }
    ## new taus/
    taustart <- tau.negbin.initial(y)
    component1start  <-  round(runif(tot.comp, 1, 10))

    ## call C++ code to draw sample
    posterior <- .C("cMCMCnegbinChange",
                    betaout = as.double(rep(0.0, nstore*ns*k)),
                    Pout = as.double(rep(0.0, nstore*ns*ns)),
                    psout = as.double(rep(0.0, n*ns)),
                    sout = as.double(rep(0.0, nstore*n)),
                    nuout = as.double(rep(0.0, nstore*n)),
                    rhoout = as.double(rep(0.0, nstore*ns)),
                    tau1out = as.double(rep(0.0, nstore*n)),
                    tau2out = as.double(rep(0.0, nstore*n)),
                    comp1out = as.integer(rep(0, nstore*n)),
                    comp2out = as.integer(rep(0, nstore*n)),
                    sr1out = as.double(rep(0, nstore*n)),
                    sr2out = as.double(rep(0, nstore*n)),
                    mr1out = as.double(rep(0, nstore*n)),
                    mr2out = as.double(rep(0, nstore*n)),
                    rhosizes = as.double(rep(0.0, ns)),
                    Ydata = as.double(y),
                    Yrow = as.integer(nrow(y)),
                    Ycol = as.integer(ncol(y)),
                    Xdata = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    m = as.integer(m),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    verbose = as.integer(verbose),
                    betastart = as.double(betastart),
                    Pstart = as.double(Pstart),
                    nustart = as.double(nu.start),
                    rhostart = as.double(rho.start),
                    tau1start = as.double(taustart$tau1),
                    tau2start = as.double(taustart$tau2),
                    component1start = as.double(component1start),
                    a = as.double(a),
                    b = as.double(b),
                    e = as.double(e),
                    f = as.double(f),
                    g = as.double(g),
                    rhostepdata = as.double(rho.step),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    b0data = as.double(b0),
                    B0data = as.double(B0),
                    A0data = as.double(A0),
                    fixed_m = as.integer(fixed.m),
                    logmarglikeholder = as.double(0.0),
                    loglikeholder = as.double(0.0),
                    chib = as.integer(chib),
                    PACKAGE="MCMCpack")



    ## pull together matrix and build MCMC object to return
    beta.holder <- matrix(posterior$betaout, nstore, )
    P.holder    <- matrix(posterior$Pout, nstore, )
    s.holder    <- matrix(posterior$sout, nstore, )
    ps.holder   <- matrix(posterior$psout, n, )
    nu.holder   <- matrix(posterior$nuout, nstore, )
    rho.holder  <- matrix(posterior$rhoout, nstore, )
    tau1.holder <- matrix(posterior$tau1out, nstore, )
    tau2.holder <- matrix(posterior$tau2out, nstore, )
    comp1.holder <- matrix(posterior$comp1out, nstore, )
    comp2.holder <- matrix(posterior$comp2out, nstore, )
    sr1.holder <- matrix(posterior$sr1out, nstore, )
    sr2.holder <- matrix(posterior$sr2out, nstore, )
    mr1.holder <- matrix(posterior$mr1out, nstore, )
    mr2.holder <- matrix(posterior$mr2out, nstore, )

    pr.chpt <- rowMeans(rbind(0, diff(t(s.holder)) != 0))
    output <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
    if (m > 0) {
      varnames(output)  <- sapply(c(1:ns),
                                  function(i) {
                                    paste(xnames, "_regime", i, sep = "")
                                  })    }
    attr(output, "title") <- "MCMCnegbinChange Posterior Sample"
    attr(output, "formula") <- formula
    attr(output, "y")       <- y
    attr(output, "X")       <- X
    attr(output, "m")       <- m
    attr(output, "call")    <- cl
    attr(output, "prob.state") <- ps.holder/nstore
    attr(output, "pr.chpt") <- pr.chpt
    attr(output, "s.store") <- s.holder
    attr(output, "P.store") <- P.holder
    attr(output, "rho.store") <- rho.holder
    attr(output, "nu.store") <- nu.holder
    ## below are auxiliary variables, sometimes useful
    ## for debugging
    ## attr(output, "tau1.store") <- tau1.holder
    ## attr(output, "tau2.store") <- tau2.holder
    ## attr(output, "comp1.store") <- comp1.holder
    ## attr(output, "comp2.store") <- comp2.holder
    ## attr(output, "sr1.store") <- sr1.holder
    ## attr(output, "sr2.store") <- sr2.holder
    ## attr(output, "mr1.store") <- mr1.holder
    ## attr(output, "mr2.store") <- mr2.holder
    attr(output, "rho.step") <- posterior$rhosizes

    ## get marginal likelihood if Chib95
    if (marginal.likelihood == "Chib95" & fixed.m){
      logmarglike <- posterior$logmarglikeholder
      loglike <- posterior$loglikeholder
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglike") <- loglike
    }
    return(output)
  }
