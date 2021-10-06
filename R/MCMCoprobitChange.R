#########################################################
##
## sample from the posterior distribution
## of ordinal probit changepoint regression model
## using a linear Gaussian approximation
##
## JHP 07/01/2007
## JHP 03/03/2009
## JHP 09/08/2010
#########################################################

#' Markov Chain Monte Carlo for Ordered Probit Changepoint Regression Model
#'
#' This function generates a sample from the posterior distribution of an
#' ordered probit regression model with multiple parameter breaks. The function
#' uses the Markov chain Monte Carlo method of Chib (1998).  The user supplies
#' data and priors, and a sample from the posterior distribution is returned as
#' an mcmc object, which can be subsequently analyzed with functions provided
#' in the coda package.
#'
#' \code{MCMCoprobitChange} simulates from the posterior distribution of an
#' ordinal probit regression model with multiple parameter breaks. The
#' simulation of latent states is based on the linear approximation method
#' discussed in Park (2011).
#'
#' The model takes the following form:
#'
#' \deqn{\Pr(y_t = 1) = \Phi(\gamma_{c, m} - x_i'\beta_m) - \Phi(\gamma_{c-1, m} - x_i'\beta_m)\;\; m = 1, \ldots, M}
#'
#' Where \eqn{M} is the number of states, and \eqn{\gamma_{c, m}} and
#' \eqn{\beta_m} are paramters when a state is \eqn{m} at \eqn{t}.
#'
#' We assume Gaussian distribution for prior of \eqn{\beta}:
#'
#' \deqn{\beta_m \sim \mathcal{N}(b_0,B_0^{-1}),\;\; m = 1, \ldots, M}
#'
#' And:
#'
#' \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}
#'
#' Where \eqn{M} is the number of states.
#'
#' Note that when the fitted changepoint model has very few observations in any
#' of states, the marginal likelihood outcome can be ``nan," which indicates
#' that too many breaks are assumed given the model and data.
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param m The number of changepoints.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burnin.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param tune The tuning parameter for the Metropolis-Hastings step. Default
#' of NA corresponds to a choice of 0.05 divided by the number of categories in
#' the response variable.
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
#' number of betas.  The default value of of NA will use the MLE estimate of
#' \eqn{\beta} as the starting value.  If this is a scalar, that value
#' will serve as the starting value mean for all of the betas.
#'
#' @param gamma.start The starting values for the \eqn{\gamma} vector.
#' This can either be a scalar or a column vector with dimension equal to the
#' number of gammas.  The default value of of NA will use the MLE estimate of
#' \eqn{\gamma} as the starting value.  If this is a scalar, that value
#' will serve as the starting value mean for all of the gammas.
#'
#' @param P.start The starting values for the transition matrix.  A user should
#' provide a square matrix with dimension equal to the number of states.  By
#' default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper
#' transition matrix for each raw except the last raw.
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
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated, and \code{Chib95} in which case the method of Chib (1995) is
#' used.
#'
#' @param gamma.fixed 1 if users want to constrain \eqn{\gamma} values
#' to be constant. By default, \eqn{\gamma} values are allowed to vary
#' across regimes.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.  The object
#' contains an attribute \code{prob.state} storage matrix that contains the
#' probability of \eqn{state_i} for each period, the log-likelihood of
#' the model (\code{loglike}), and the log-marginal likelihood of the model
#' (\code{logmarglike}).
#'
#' @export
#'
#' @seealso \code{\link{plotState}}, \code{\link{plotChangepoint}}
#'
#' @references Jong Hee Park. 2011. ``Changepoint Analysis of Binary and
#' Ordinal Probit Models: An Application to Bank Rate Policy Under the Interwar
#' Gold Standard."  \emph{Political Analysis}. 19: 188-204. <doi:10.1093/pan/mpr007>
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point
#' models.'' \emph{Journal of Econometrics}. 86: 221-241.
#'
#' @keywords models
#'
#' @examples
#'
#' set.seed(1909)
#' N <- 200
#' x1 <- rnorm(N, 1, .5);
#'
#' ## set a true break at 100
#' z1 <- 1 + x1[1:100] + rnorm(100);
#' z2 <- 1 -0.2*x1[101:200] + rnorm(100);
#' z <- c(z1,  z2);
#' y <- z
#'
#' ## generate y
#' y[z < 1] <- 1;
#' y[z >= 1 & z < 2] <- 2;
#' y[z >= 2] <- 3;
#'
#' ## inputs
#' formula <- y ~ x1
#'
#' ## fit multiple models with a varying number of breaks
#' out1 <- MCMCoprobitChange(formula, m=1,
#'       	mcmc=100, burnin=100, thin=1, tune=c(.5, .5), verbose=100,
#'      	b0=0, B0=0.1, marginal.likelihood = "Chib95")
#' out2 <- MCMCoprobitChange(formula, m=2,
#'       	mcmc=100, burnin=100, thin=1, tune=c(.5, .5, .5), verbose=100,
#'      	b0=0, B0=0.1, marginal.likelihood = "Chib95")
#'
#' ## Do model comparison
#' ## NOTE: the chain should be run longer than this example!
#' BayesFactor(out1, out2)
#'
#' ## draw plots using the "right" model
#' plotState(out1)
#' plotChangepoint(out1)
#'
"MCMCoprobitChange"<-
  function(formula, data=parent.frame(),  m = 1,
           burnin = 1000, mcmc = 1000, thin = 1, tune = NA, verbose = 0,
           seed = NA, beta.start = NA, gamma.start = NA, P.start = NA,
           b0 = NULL, B0 = NULL, a = NULL, b = NULL,
           marginal.likelihood = c("none", "Chib95"),
           gamma.fixed=0, ...){

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    cl <- match.call()
    nstore <- mcmc/thin

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    totiter <- mcmc+burnin
    holder <- parse.formula(formula, data=data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    K <- ncol(X)
    Y <- factor(y, ordered = TRUE)
    ncat <- nlevels(Y)
    cat <- levels(Y)
    ns <- m + 1
    N <- nrow(X)
    gk <- ncat + 1

    if(sum(is.na(tune))==1) {
      stop("Please specify a tune parameter and call MCMCoprobitChange() again.\n")
    }
    else if (length(tune)==1){
      tune <- rep(tune, ns)
    }
    else if(length(tune)>1&length(tune)<ns){
      tune <- rep(tune[1], ns)
      cat("The first element of tune is repeated to make it conformable to the number of states.\n")
    }
    else{

    }

    xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if (xint > 0) {
      new.X <- X[, -xint, drop = FALSE]
    }
    else
      warning("An intercept is needed and assumed in MCMCoprobitChange()\n.")
    if (ncol(new.X) == 0) {
      polr.out <- polr(ordered(Y) ~ 1)
    }
    else {
      polr.out <- polr(ordered(Y) ~ new.X)
    }

    ## prior for transition matrix
    A0 <- trans.mat.prior(m=m, n=N, a=a, b=b)

    ## prior for beta error checking
    if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,K,1)
    }
    if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
      cat("N(b0,B0) prior b0 not conformable.\n")
      stop("Please respecify and call MCMCoprobitChange() again.\n")
    }
    if(is.null(dim(B0))) {
      B0 <- B0 * diag(K)
    }
    if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
      cat("N(b0,B0) prior B0 not conformable.\n")
      stop("Please respecify and call MCMCoprobitChange() again.\n")
    }
    marginal.likelihood  <- match.arg(marginal.likelihood)
    B0.eigenvalues <- eigen(B0)$values
     if (isTRUE(all.equal(min(B0.eigenvalues), 0))){
      if (marginal.likelihood != "none"){
        warning("Cannot calculate marginal likelihood with improper prior\n")
        marginal.likelihood <- "none"
      }
    }
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    ## to save time
    B0inv <- solve(B0)
    gamma.start <- matrix(NA, ncat + 1, 1)
    gamma.start[1] <- -300
    gamma.start[2] <- 0
    gamma.start[3:ncat] <- (polr.out$zeta[2:(ncat - 1)] - polr.out$zeta[1]) * 0.588
    gamma.start[ncat + 1] <- 300

    ## initial values
    mle <- polr(Y ~ X[,-1])
    beta <- matrix(rep(c(mle$zeta[1], coef(mle)), ns), ns, , byrow=TRUE)
    ols <- lm(as.double(Y) ~ X-1)
    betalinearstart <- matrix(rep(coef(ols), ns), ns, , byrow=TRUE)
    P <-  trans.mat.prior(m=m, n=N, a=0.9, b=0.1)
    Sigmastart <- summary(ols)$sigma
    if (gamma.fixed==1){
      gamma <- gamma.start
      gamma.storage <-rep(0.0, nstore*gk)
    }
    else {
      gamma <- matrix(rep(gamma.start, ns), ns, ,  byrow=T)
      gamma.storage <- rep(0.0, nstore*ns*gk)
    }

    ## call C++ code to draw sample
    posterior <- .C("cMCMCoprobitChange",
                    betaout = as.double(rep(0.0, nstore*ns*K)),
                    betalinearout = as.double(rep(0.0, nstore*ns*K)),
                    gammaout = as.double(gamma.storage),
                    Pout = as.double(rep(0.0, nstore*ns*ns)),
                    psout = as.double(rep(0.0, N*ns)),
                    sout = as.double(rep(0.0, nstore*N)),

                    Ydata = as.double(Y),
                    Xdata = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),

                    m = as.integer(m),
                    ncat = as.integer(ncat),

                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    verbose = as.integer(verbose),

                    tunedata = as.double(tune),
                    lecuyer=as.integer(lecuyer),
                    seedarray=as.integer(seed.array),
                    lecuyerstream=as.integer(lecuyer.stream),

                    betastart = as.double(beta),
                    betalinearstart = as.double(betalinearstart),
                    gammastart = as.double(gamma),
                    Pstart = as.double(P),
                    sigmastart = as.double(Sigmastart),

                    a = as.double(a),
                    b = as.double(b),
                    b0data = as.double(b0),
                    B0data = as.double(B0),
                    A0data = as.double(A0),
                    logmarglikeholder = as.double(0.0),
                    loglikeholder = as.double(0.0),
                    chib = as.integer(chib),
                    gammafixed= as.integer(gamma.fixed))

    ## get marginal likelihood if Chib95
    if (chib==1){
      logmarglike <- posterior$logmarglikeholder
      loglike <- posterior$loglikeholder
    }
    else{
      logmarglike <- loglike <- 0
    }

    ## pull together matrix and build MCMC object to return
    beta.holder <- mcmc(matrix(posterior$betaout, nstore, ns*K))
    if (gamma.fixed==1){
      gamma.holder <- mcmc(matrix(posterior$gammaout, nstore, gk))
    }
    else {
      gamma.holder <- mcmc(matrix(posterior$gammaout, nstore, ns*gk))
    }
    P.holder    <- matrix(posterior$Pout, nstore, )
    s.holder    <- matrix(posterior$sout, nstore, )
    ps.holder   <- matrix(posterior$psout, N, )

    varnames(beta.holder)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c(xnames), "_regime", i, sep = "")
                                     })
    ## betalinear
    betalinear.holder <- mcmc(matrix(posterior$betalinearout, nstore, ns*K))
    varnames(betalinear.holder)  <- sapply(c(1:ns),
                                           function(i){
                                             paste(c(xnames), "_regime", i, sep = "")
                                           })
    gamma.holder <- gamma.holder[, as.vector(sapply(1:ns, function(i){gk*(i-1) + (3:(gk-1))}))]
    gamma.names <- paste("gamma", 3:(gk-1), sep="")
    varnames(gamma.holder)  <- sapply(c(1:ns),
                                      function(i){
                                        paste(gamma.names, "_regime", i, sep = "")
                                      })

    output <- mcmc(cbind(beta.holder, gamma.holder))
    attr(output, "title") <- "MCMCoprobitChange Posterior Sample"
    ## attr(output, "betalinear") <- mcmc(betalinear.holder)
    attr(output, "formula") <- formula
    attr(output, "y")       <- Y
    attr(output, "X")       <- X
    attr(output, "m")       <- m
    attr(output, "call")    <- cl
    attr(output, "logmarglike") <- logmarglike
    attr(output, "loglike") <- loglike
    attr(output, "prob.state") <- ps.holder/nstore
    attr(output, "s.store") <- s.holder
    return(output)

  }## end of MCMC function
