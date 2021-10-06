##########################################################################
## sample from the posterior distribution of a Negative Binomial regression
## model in R using linked C++ code in Scythe
##
## 05/20/2017 Matthew Blackwell
##########################################################################

#' Markov Chain Monte Carlo for Negative Binomial Regression
#'
#' This function generates a sample from the posterior distribution of a
#' Negative Binomial regression model via auxiliary mixture sampling. The user
#' supplies data and priors, and a sample from the posterior distribution is
#' returned as an mcmc object, which can be subsequently analyzed with
#' functions provided in the coda package.
#'
#' \code{MCMCnegbin} simulates from the posterior distribution of a
#' Negative Binomial regression model using a combination of auxiliary
#' mixture sampling and slice sampling. The simulation proper is done
#' in compiled C++ code to maximize efficiency. Please consult the
#' coda documentation for a comprehensive list of functions that can
#' be used to analyze the posterior sample.
#'
#' The model takes the following form:
#'
#' \deqn{y_i \sim \mathcal{P}oisson(\nu_i\mu_i)}
#'
#' Where the inverse link function:
#'
#' \deqn{\mu_i = \exp(x_i'\beta)}
#'
#' We assume a multivariate Normal prior on \eqn{\beta}:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}
#'
#' The unit-level random effect that handles overdispersion is assumed
#' to be distributed Gamma:
#'
#' \deqn{\nu_i \sim \mathcal{G}amma(\rho, \rho)}
#'
#' The overdispersion parameter has a prior with the following form:
#'
#' \deqn{f(\rho|e,f,g) \propto \rho^{e-1}(\rho + g)^{-(e+f)}}
#'
#' The model is simulated via blocked Gibbs, with the the \eqn{\beta}
#' being simulated via the auxiliary mixture sampling method of
#' Fuerhwirth-Schanetter et al. (2009). The \eqn{\rho} is updated via
#' slice sampling. The \eqn{\nu_i} are updated their (conjugate) full
#' conditional, which is also Gamma.
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
#' @param rho.start The starting value for the \eqn{\rho} variable.
#'     The default value is 1.
#'
#' @param nu.start The starting values for the random effect,
#'     \eqn{\nu}. The default value is a vector of ones.
#'
#' @param e The hyperprior for the distribution \eqn{\rho}. See details.
#'
#' @param f The hyperprior for the distribution \eqn{\rho}. See details.
#'
#' @param g The hyperprior for the distribution \eqn{\rho}. See details.
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
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MASS]{glm.nb}}
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
#' Sylvia Fruehwirth-Schnatter, Rudolf Fruehwirth, Leonhard Held, and
#'     Havard Rue. 2009. ``Improved auxiliary mixture sampling for
#'     hierarchical models of non-Gaussian data'', \emph{Statistics
#'     and Computing} 19(4): 479-492.
#'     <doi:10.1007/s11222-008-9109-4>
#'
#' @keywords models
#'
#' @examples
#'
#'  \dontrun{
#'    n <- 150
#'    mcmcs <- 5000
#'    burnin <- 5000
#'    thin <- 5
#'    x1 <- runif(n, 0, 2)
#'    rho.true <- 1.5
#'    nu.true <- rgamma(n, rho.true, rho.true)
#'    mu <- nu.true * exp(1 + x1 * 1)
#'    y <- rpois(n, mu)
#'    posterior <- MCMCnegbin(y ~ x1)
#'    plot(posterior)
#'    summary(posterior)
#'    }
#'
"MCMCnegbin"<-
  function(formula, data = parent.frame(),
           b0 = 0, B0 = 1,
           e = 2, f = 2, g = 10,
           burnin = 1000, mcmc = 1000, thin = 1, verbose = 0,
           seed = NA, beta.start = NA,
           rho.start = NA, rho.step = 0.1, nu.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...) {

    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    k <- ncol(X)
    n <- length(y)
    n.arrival<-  y + 1
    tot.comp <-  n + sum(y)
    NT      <-  max(n.arrival)

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

    betastart  <- beta.negbin.start(beta.start, 1, k, formula, data)
    if (is.na(rho.start[1])) {
      rho.start <- 1
    }
    if (is.na(nu.start[1])) {
      nu.start <- rep(1, times = n)
    }
    ## new taus/
    taustart <- tau.negbin.initial(y)
    component1start  <-  round(runif(tot.comp, 1, 10))

    ## call C++ code to draw sample
    posterior <- .C("cMCMCnegbin",
                    betaout = as.double(rep(0.0, nstore*k)),
                    nuout = as.double(rep(0.0, nstore*n)),
                    rhoout = as.double(rep(0.0, nstore)),
                    tau1out = as.double(rep(0.0, nstore*n)),
                    tau2out = as.double(rep(0.0, nstore*n)),
                    comp1out = as.integer(rep(0, nstore*n)),
                    comp2out = as.integer(rep(0, nstore*n)),
                    sr1out = as.double(rep(0, nstore*n)),
                    sr2out = as.double(rep(0, nstore*n)),
                    mr1out = as.double(rep(0, nstore*n)),
                    mr2out = as.double(rep(0, nstore*n)),
                    rhosizes = as.double(0.0),
                    Ydata = as.double(y),
                    Yrow = as.integer(nrow(y)),
                    Ycol = as.integer(ncol(y)),
                    Xdata = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    verbose = as.integer(verbose),
                    betastart = as.double(betastart),
                    nustart = as.double(nu.start),
                    rhostart = as.double(rho.start),
                    tau1start = as.double(taustart$tau1),
                    tau2start = as.double(taustart$tau2),
                    component1start = as.double(component1start),
                    e = as.double(e),
                    f = as.double(f),
                    g = as.double(g),
                    rhostepdata = as.double(rho.step),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    b0data = as.double(b0),
                    B0data = as.double(B0),
                    logmarglikeholder = as.double(0.0),
                    loglikeholder = as.double(0.0),
                    chib = as.integer(chib),
                    PACKAGE="MCMCpack")



    ## pull together matrix and build MCMC object to return
    beta.holder <- matrix(posterior$betaout, nstore, )
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

    output <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
    varnames(output) <- xnames
    attr(output, "title") <- "MCMCnegbin Posterior Sample"
    attr(output, "formula") <- formula
    attr(output, "y")       <- y
    attr(output, "X")       <- X
    attr(output, "call")    <- cl
    attr(output, "rho.store") <- rho.holder
    attr(output, "nu.store") <- nu.holder
    attr(output, "tau1.store") <- tau1.holder
    attr(output, "tau2.store") <- tau2.holder
    attr(output, "comp1.store") <- comp1.holder
    attr(output, "comp2.store") <- comp2.holder
    attr(output, "sr1.store") <- sr1.holder
    attr(output, "sr2.store") <- sr2.holder
    attr(output, "mr1.store") <- mr1.holder
    attr(output, "mr2.store") <- mr2.holder
    attr(output, "rho.step") <- posterior$rhosizes

    ## get marginal likelihood if Chib95
    if (marginal.likelihood == "Chib95"){
      logmarglike <- posterior$logmarglikeholder
      loglike <- posterior$loglikeholder
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglike") <- loglike
    }

    return(output)
  }


## initial values of tau in negbin samplers
"tau.negbin.initial" <- function(y) {
  tau1 <- rep(NA, length(y))
  tau2 <- rep(NA, length(y))
  lambda.t        <-  exp(5)
  for (t in 1:length(y)){
    nt      <-  y[t]
    if (nt==0) {
      tau1[t]    <-  1 + rexp(1, lambda.t)
      tau2[t] <- 0
    }
    else{
      xi <- rexp(1, lambda.t)
      tau2[t] <- rbeta(1, nt, 1)
      tau1[t] <- 1 - tau2[t] + xi
    }
  }
  return(list(tau1 = tau1, tau2 = tau2))
}

## beta starting values in negbin samplers
"beta.negbin.start"<- function (beta.start, ns, k, formula, data){
  ## if a user does not specify beta.start, use a coefficient vector from mle
  if (is.na(beta.start[1])) {
    b0 <- coef(glm.nb(formula, data = data))
    beta.start  <-  matrix(rep(b0, ns), ns, k, byrow=TRUE)
  }
  ## if beta.start is scalar or k by 1 vector, repeat this
  else if (is.null(dim(beta.start))&&length(beta.start)<=k) {
    beta.start <- beta.start * matrix(1, ns, k)
    ## this alternates beta.start if beta.start is not a scalar
  }
  ## if the length of beta.start is same to ns*k, make this as a matrix
  else if (is.null(dim(beta.start))&&length(beta.start)==ns*k) {
    beta.start <- matrix(beta.start, ns, k)
  }
  else if (is.null(dim(beta.start))&&length(beta.start)>=k) {
    cat("Error: Starting value for beta not conformable.\n")
    stop("Please respecify and call ", calling.function(),
         " again.\n", call. = FALSE)
  }
  ## else, report an error message and stop
  else if (!all(dim(beta.start) == c(ns, k))) {
    cat("Error: Starting value for beta not conformable.\n")
    stop("Please respecify and call ", calling.function(),
         " again.\n", call. = FALSE)
  }
  return(beta.start)
}
