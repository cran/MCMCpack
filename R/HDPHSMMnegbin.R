################################
## HDP-HMM Negative Binomial Model
##
## Matthew Blackwell 5/22/2017
################################

#' Markov Chain Monte Carlo for HDP-HSMM with a Negative
#' Binomial outcome distribution
#'
#' This function generates a sample from the posterior distribution of
#' a Hidden Semi-Markov Model with a Heirarchical Dirichlet Process
#' and a Negative Binomial outcome distribution
#' (Johnson and Willsky, 2013).  The user supplies data and priors, and a
#' sample from the posterior distribution is returned as an mcmc
#' object, which can be subsequently analyzed with functions provided
#' in the coda package.
#'
#' \code{HDPHSMMnegbin} simulates from the posterior distribution of a
#' HDP-HSMM with a Negative Binomial outcome distribution,
#' allowing for multiple, arbitrary changepoints in the model. The details of the
#' model are discussed in Johnson & Willsky (2013). The implementation here is
#' based on a weak-limit approximation, where there is a large, though
#' finite number of regimes that can be switched between. Unlike other
#' changepoint models in \code{MCMCpack}, the HDP-HSMM approach allows
#' for the state sequence to return to previous visited states.
#'
#'  The model takes the following form, where we show the fixed-limit version:
#'
#' \deqn{y_t \sim \mathcal{P}oisson(\nu_t\mu_t)}
#'
#' \deqn{\mu_t = x_t ' \beta_k,\;\; k = 1, \ldots, K}
#'
#' \deqn{\nu_t \sim \mathcal{G}amma(\rho_k, \rho_k)}
#'
#' Where \eqn{K} is an upper bound on the number of states and
#' \eqn{\beta_k} and \eqn{\rho_k} are parameters when a state is
#' \eqn{k} at \eqn{t}.
#'
#' In the HDP-HSMM, there is a super-state sequence that, for a given
#' observation, is drawn from the transition distribution and then a
#' duration is drawn from a duration distribution to determin how long
#' that state will stay active. After that duration, a new super-state
#' is drawn from the transition distribution, where self-transitions
#' are  disallowed. The transition probabilities between states are
#' assumed to follow a heirarchical Dirichlet process:
#'
#' \deqn{\pi_k \sim \mathcal{D}irichlet(\alpha\delta_1, \ldots ,
#' \alpha\delta_K)}
#'
#' \deqn{\delta \sim \mathcal{D}irichlet(\gamma/K, \ldots, \gamma/K)}
#'
#' In the algorithm itself, these \eqn{\pi} vectors are modified to
#' remove self-transitions as discussed above. There is a unique
#' duration distribution for each regime with the following
#' parameters:
#'
#' \deqn{D_k \sim \mathcal{N}egBin(r, \omega_k)}
#'
#' \deqn{\omega_k \sim \mathcal{B}eta(a_{\omega,k}, b_{\omega, k})}
#'
#' We assume Gaussian distribution for prior of \eqn{\beta}:
#'
#' \deqn{\beta_k \sim \mathcal{N}(b_0,B_0^{-1}),\;\; m = 1, \ldots, K}
#'
#' The overdispersion parameters have a prior with the following form:
#'
#' \deqn{f(\rho_k|e,f,g) \propto \rho^{e-1}(\rho + g)^{-(e+f)}}
#'
#' The model is simulated via blocked Gibbs conditonal on the states.
#' The \eqn{\beta} being simulated via the auxiliary mixture sampling
#' method of Fuerhwirth-Schanetter et al. (2009). The \eqn{\rho} is
#' updated via slice sampling. The \eqn{\nu_t} are updated their
#' (conjugate) full conditional, which is also Gamma. The states and
#' their durations are drawn as in Johnson & Willsky (2013).
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param K The number of regimes under consideration. This should be
#'   larger than the hypothesized number of regimes in the data. Note
#'   that the sampler will likely visit fewer than \code{K} regimes.
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
#' @param r Parameter of the Negative Binomial prior for regime
#'   durations. It is the target number of successful trials. Must be
#'   strictly positive. Higher values increase the variance of the
#'   duration distributions.
#'
#' @param a.omega,b.omega Paramaters for the Beta prior on
#'   \eqn{\omega}, which determines the regime length distribution,
#'   which is Negative Binomial, with parameters \code{r} and \code{omega}.

#' @param beta.start The starting value for the \eqn{\beta} vector.  This
#' can either be a scalar or a column vector with dimension equal to the number
#' of betas. If this takes a scalar value, then that value will serve as the
#' starting value for all of the betas.  The default value of NA will use the
#' maximum likelihood estimate of \eqn{\beta} as the starting value
#'   for all regimes.
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
#'
#' @param a.alpha,b.alpha Shape and scale parameters for the Gamma
#'   distribution on \eqn{\alpha}.
#'
#' @param a.gamma,b.gamma Shape and scale parameters for the Gamma
#'   distribution on \eqn{\gamma}.
#'
#' @param alpha.start,gamma.start Scalar starting values for the
#'   \eqn{\alpha}, and \eqn{\gamma} parameters.
#'
#' @param omega.start A vector of starting values for the probability
#'   of success parameter in the Negative Binomial distribution that
#'   governs the duration distributions.
#'
#' @param P.start Initial transition matrix between regimes. Should be
#'   a \code{K} by \code{K} matrix. If not provided, the default value
#'   will be uniform transition distributions.
#'
#' @param ... further arguments to be passed.
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link{MCMCnegbinChange}},
#'   \code{\link{HDPHMMnegbin}},
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
#' Matthew J. Johnson and Alan S. Willsky. 2013. ``Bayesian Nonparametric Hidden Semi-Markov Models.'' \emph{Journal of Machine Learning Research}, 14(Feb), 673-701.
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
#'    b1.true <- c(1, -2, 2)
#'    x1 <- runif(n, 0, 2)
#'    nu.true <- rgamma(n, rho.true[true.s], rho.true[true.s])
#'    mu <- nu.true * exp(1 + x1 * b1.true[true.s])
#'    y <- rpois(n, mu)
#'
#'    posterior <- HDPHSMMnegbin(y ~ x1, K = 10, verbose = 1000,
#'                           e = 2, f = 2, g = 10,
#'                           b0 = 0, B0 = 1/9,
#'                           a.omega = 1, b.omega = 100, r = 1,
#'                           rho.step = rep(0.75, times = 10),
#'                           seed = list(NA, 2),
#'                           omega.start = 0.05, gamma.start = 10,
#'                           alpha.start = 5)
#'
#'    plotHDPChangepoint(posterior, ylab="Density", start=1)
#'    }
#'



"HDPHSMMnegbin"<-
  function(formula, data = parent.frame(), K = 10,
           b0 = 0, B0 = 1,
           a.alpha = 1, b.alpha = 0.1, a.gamma = 1, b.gamma = 0.1,
           a.omega, b.omega, e = 2, f = 2, g = 10, r = 1,
           burnin = 1000, mcmc = 1000, thin = 1, verbose = 0,
           seed = NA, beta.start = NA, P.start = NA,
           rho.start = NA, rho.step, nu.start = NA,
           omega.start = NA, gamma.start = 0.5, alpha.start = 100, ...) {

    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    J <- ncol(X)
    n <- length(y)
    n.arrival<-  y + 1
    NT      <-  max(n.arrival)
    tot.comp <-  n + sum(y)

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


    mvn.prior <- form.mvn.prior(b0, B0, J)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    if (K == 1) {
      stop("HDP-HSMM only works with K > 1")
    }

    ## get initial values of tau from observed y
    Pstart <- check.hdp.P(P.start, K, n, 1/K)

    betastart  <- beta.negbin.start(beta.start, K, J, formula, data)
    if (is.na(rho.start[1])) {
      rho.start <- rep(1, times = K)
    }
    if (!(length(rho.step) %in% c(1, K))) {
      stop("rho.step has the wrong dimensions")
    }
    if (length(rho.step) == 1) {
      rho.step <- rep(rho.step, times = K)
    }
    if (is.na(nu.start[1])) {
      nu.start <- rep(1, times = n)
    }
    if (!(length(omega.start) %in% c(1, K))) {
      stop("omega.start has the wrong dimensions")
    }
    if (is.na(omega.start[1])) {
      omega.start <- rbeta(K, a.omega, b.omega)
    } else {
      if (length(omega.start) == 1) {
        omega.start <- rep(omega.start, times = K)
      }
    }
    ## new taus
    taustart <- tau.negbin.initial(y)
    component1start  <-  round(runif(tot.comp, 1, 10))

    ## call C++ code to draw sample
    posterior <- .C("cHDPHSMMnegbin",
                    betaout = as.double(rep(0.0, nstore * K * J)),
                    Pout = as.double(rep(0.0, nstore * K * K)),
                    omegaout = as.double(rep(0.0, nstore * K)),
                    sout = as.double(rep(0.0, nstore * n)),
                    nuout = as.double(rep(0.0, nstore * n)),
                    rhoout = as.double(rep(0.0, nstore * K)),
                    tau1out = as.double(rep(0.0, nstore * n)),
                    tau2out = as.double(rep(0.0, nstore * n)),
                    comp1out = as.integer(rep(0, nstore * n)),
                    comp2out = as.integer(rep(0, nstore * n)),
                    sr1out = as.double(rep(0, nstore * n)),
                    sr2out = as.double(rep(0, nstore * n)),
                    mr1out = as.double(rep(0, nstore*n)),
                    mr2out = as.double(rep(0, nstore*n)),
                    gammaout = as.double(rep(0, nstore)),
                    alphaout = as.double(rep(0, nstore)),
                    rhosizes = as.double(rep(0.0, K)),
                    Ydata = as.double(y),
                    Yrow = as.integer(nrow(y)),
                    Ycol = as.integer(ncol(y)),
                    Xdata = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    K = as.integer(K),
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
                    alphastart = as.double(alpha.start),
                    gammastart = as.double(gamma.start),
                    omegastart = as.double(omega.start),
                    a_alpha = as.double(a.alpha),
                    b_alpha = as.double(b.alpha),
                    a_gamma = as.double(a.gamma),
                    b_gamma = as.double(b.gamma),
                    a_omega = as.double(a.omega),
                    b_omega = as.double(b.omega),
                    e = as.double(e),
                    f = as.double(f),
                    g = as.double(g),
                    r = as.double(r),
                    rhostepdata = as.double(rho.step),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    b0data = as.double(b0),
                    B0data = as.double(B0),
                    PACKAGE="MCMCpack")



    ## pull together matrix and build MCMC object to return
    beta.holder <- matrix(posterior$betaout, nstore, )
    P.holder    <- matrix(posterior$Pout, nstore, )
    s.holder    <- matrix(posterior$sout, nstore, )
    nu.holder   <- matrix(posterior$nuout, nstore, )
    rho.holder  <- matrix(posterior$rhoout, nstore, )
    omega.holder  <- matrix(posterior$omegaout, nstore, )
    tau1.holder <- matrix(posterior$tau1out, nstore, )
    tau2.holder <- matrix(posterior$tau2out, nstore, )
    comp1.holder <- matrix(posterior$comp1out, nstore, )
    comp2.holder <- matrix(posterior$comp2out, nstore, )
    sr1.holder <- matrix(posterior$sr1out, nstore, )
    sr2.holder <- matrix(posterior$sr2out, nstore, )
    mr1.holder <- matrix(posterior$mr1out, nstore, )
    mr2.holder <- matrix(posterior$mr2out, nstore, )

    nsts <- apply(s.holder, 1, function(x) length(unique(x)))
    mode.mod <- as.numeric(names(which.max(table(nsts))))
    pr.chpt <- rowMeans(rbind(0, diff(t(s.holder)) != 0))
    output <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output)  <- sapply(c(1:K), function(i) { paste(xnames,
    "_regime", i, sep = "")})
    attr(output, "title") <- "HDPHSMMnegbin Posterior Sample"
    attr(output, "formula") <- formula
    attr(output, "y")       <- y
    attr(output, "X")       <- X
    attr(output, "K")       <- K
    attr(output, "call")    <- cl
    ## attr(output, "prob.state") <- ps.holder/nstore
    attr(output, "pr.chpt") <- pr.chpt
    attr(output, "s.store") <- s.holder
    attr(output, "P.store") <- P.holder
    attr(output, "rho.store") <- rho.holder
    attr(output, "nu.store") <- nu.holder
    ## attr(output, "tau1.store") <- tau1.holder
    ## attr(output, "tau2.store") <- tau2.holder
    ## attr(output, "comp1.store") <- comp1.holder
    ## attr(output, "comp2.store") <- comp2.holder
    ## attr(output, "sr1.store") <- sr1.holder
    ## attr(output, "sr2.store") <- sr2.holder
    ## attr(output, "mr1.store") <- mr1.holder
    ## attr(output, "mr2.store") <- mr2.holder
    attr(output, "rho.step") <- posterior$rhosizes
    attr(output, "num.regimes") <- nsts
    attr(output, "gamma.store") <- posterior$gammaout
    attr(output, "omega.store") <- omega.holder
    attr(output, "alpha.store") <- posterior$alphaout

    return(output)
  }


## initial values of tau in MCMCnegbinChange
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
