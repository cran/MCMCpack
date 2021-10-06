##########################################################################
## MCMCSVDreg.R samples from the posterior distribution of a Gaussian
## linear regression model in which the X matrix has been decomposed
## with an SVD. Useful for prediction when number of columns of X
## is (possibly much) greater than the number of rows of X.
##
## See West, Mike. 2000. "Bayesian Regression in the 'Large p, Small n'
##      Paradigm." Duke ISDS Discussion Paper 2000-22.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 9/9/2005
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

parse.formula.SVDreg <- function(formula, data, intercept){

  ## extract Y, X, and variable names for model formula and frame
  mt <- terms(formula, data=data)
  if(missing(data)) data <- sys.frame(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- NULL
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass ## for specialty handling of missing data

  mf[[1]] <- as.name("model.frame")

   mf <- eval(mf, sys.frame(sys.parent()))
  if (!intercept){
    attributes(mt)$intercept <- 0
  }

  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)
  X <- as.matrix(X)         # X matrix
  Y <- as.matrix(model.response(mf, "numeric")) # Y matrix


  ## delete obs that are missing in X but potentially keep obs that are
  ## missing Y
  ## These are used to for the full SVD of X'
  keep.indic <- apply(is.na(X), 1, sum) == 0
  Y.full <- as.matrix(Y[keep.indic,])
  X.full <- X[keep.indic,]
  xvars.full <- dimnames(X.full)[[2]] # X variable names
  xobs.full  <- dimnames(X.full)[[1]] # X observation names

  return(list(Y.full, X.full, xvars.full, xobs.full))

}

#' Markov Chain Monte Carlo for SVD Regression
#'
#' This function generates a sample from the posterior distribution of
#' a linear regression model with Gaussian errors in which the design
#' matrix has been decomposed with singular value decomposition.The
#' sampling is done via the Gibbs sampling algorithm.  The user
#' supplies data and priors, and a sample from the posterior
#' distribution is returned as an mcmc object, which can be
#' subsequently analyzed with functions provided in the coda package.
#'
#' The model takes the following form: \deqn{y = X \beta +
#' \varepsilon} Where the errors are assumed to be iid Gaussian:
#' \deqn{\varepsilon_{i} \sim \mathcal{N}(0, \sigma^2)}
#'
#' Let \eqn{N} denote the number of rows of \eqn{X} and \eqn{P} the
#' number of columns of \eqn{X}. Unlike the standard regression setup
#' where \eqn{N >> P} here it is the case that \eqn{P >> N}.
#'
#' To deal with this problem a singular value decomposition of
#' \eqn{X'} is performed: \eqn{X' = ADF} and the regression model
#' becomes
#'
#' \deqn{y = F'D \gamma + \varepsilon}
#'
#' where \eqn{\gamma = A' \beta}
#'
#' We assume the following priors:
#'
#' \deqn{\sigma^{-2} \sim \mathcal{G}amma(a_0/2, b_0/2)}
#'
#' \deqn{\tau^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}
#'
#' \deqn{\gamma_i \sim w0_i \delta_0 + (1-w0_i) \mathcal{N}(g0_i,
#' \sigma^2 \tau_i^2/ d_i^2)}
#'
#' where \eqn{\delta_0} is a unit point mass at 0 and \eqn{d_i} is the
#' \eqn{i}th diagonal element of \eqn{D}.
#'
#' @param formula Model formula. Predictions are returned for elements
#'   of y that are coded as NA.
#'
#' @param data Data frame.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burnin.
#'
#' @param thin The thinning interval used in the simulation.  The
#'   number of MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the
#'   progress of the sampler is printed to the screen.  If
#'   \code{verbose} is greater than 0 the iteration number, the
#'   \eqn{\beta} vector, and the error variance are printed to the
#'   screen every \code{verbose}th iteration.
#'
#' @param seed The seed for the random number generator.  If NA, the
#'   Mersenne Twister generator is used with default seed 12345; if an
#'   integer is passed it is used to seed the Mersenne twister.  The
#'   user can also pass a list of length two to use the L'Ecuyer
#'   random number generator, which is suitable for parallel
#'   computation.  The first element of the list is the L'Ecuyer seed,
#'   which is a vector of length six or NA (if NA a default seed of
#'   \code{rep(12345,6)} is used).  The second element of list is a
#'   positive substream number. See the MCMCpack specification for
#'   more details.
#'
#' @param tau2.start The starting values for the \eqn{\tau^2} vector.
#'   Can be either a scalar or a vector. If a scalar is passed then
#'   that value will be the starting value for all elements of
#'   \eqn{\tau^2}.
#'
#' @param g0 The prior mean of \eqn{\gamma}.  This can either be a
#'   scalar or a column vector with dimension equal to the number of
#'   gammas. If this takes a scalar value, then that value will serve
#'   as the prior mean for all of the betas.
#'
#' @param a0 \eqn{a_0/2} is the shape parameter for the inverse Gamma
#'   prior on \eqn{\sigma^2} (the variance of the disturbances). The
#'   amount of information in the inverse Gamma prior is something
#'   like that from \eqn{a_0} pseudo-observations.
#'
#' @param b0 \eqn{b_0/2} is the scale parameter for the inverse Gamma
#'   prior on \eqn{\sigma^2} (the variance of the disturbances). In
#'   constructing the inverse Gamma prior, \eqn{b_0} acts like the sum
#'   of squared errors from the \eqn{a_0} pseudo-observations.
#'
#' @param c0 \eqn{c_0/2} is the shape parameter for the inverse Gamma
#'   prior on \eqn{\tau_i^2}.
#'
#' @param d0 \eqn{d_0/2} is the scale parameter for the inverse Gamma
#'   prior on \eqn{\tau_i^2}.
#'
#' @param w0 The prior probability that \eqn{\gamma_i = 0}.  Can be
#'   either a scalar or an \eqn{N} vector where \eqn{N} is the number
#'   of observations.
#'
#' @param beta.samp Logical indicating whether the sampled elements of
#'   beta should be stored and returned.
#'
#' @param intercept Logical indicating whether the original design
#'   matrix should include a constant term.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This
#'   object can be summarized by functions provided by the coda
#'   package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},
#'   \code{\link[coda]{summary.mcmc}}, \code{\link[stats]{lm}}
#'
#' @references Mike West, Josheph Nevins, Jeffrey Marks, Rainer Spang,
#'   and Harry Zuzan. 2000. ``DNA Microarray Data Analysis and
#'   Regression Modeling for Genetic Expression
#'   Profiling." Duke ISDS working paper.
#'
#' Gottardo, Raphael, and Adrian Raftery. 2004. ``Markov chain Monte
#' Carlo with mixtures of singular distributions.'' Statistics
#' Department, University of Washington, Technical Report 470.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#' ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of
#' Statistical Software}.  42(9): 1-21.
#' \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.
#' \emph{Scythe Statistical Library 1.0.}
#' \url{http://scythe.lsa.umich.edu}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.
#' ``Output Analysis and Diagnostics for MCMC (CODA)'', \emph{R
#' News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' @keywords models
"MCMCSVDreg" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA, tau2.start = 1,
           g0 = 0, a0 = 0.001, b0 = 0.001, c0=2, d0=2, w0=1,
           beta.samp=FALSE, intercept=TRUE, ...) {

    # checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    # seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula.SVDreg(formula, data, intercept)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    obsnames <- holder[[4]]
    K <- ncol(X)  # number of covariates in unidentified model
    N <- nrow(X)  # number of obs (including obs for which predictions
                  # are required) N is also the length of gamma
    Y.miss.indic <- as.numeric(is.na(Y))
    n.miss.Y <- sum(Y.miss.indic)
    Y[is.na(Y)] <- mean(Y, na.rm=TRUE)


    ## create SVD representation of t(X)
    svd.out <- svd(t(X)) # t(X) = A %*% D %*% F
    A <- svd.out$u
    D <- diag(svd.out$d)
    F <- t(svd.out$v)

    ## starting values and priors
    if (length(tau2.start) < N){
      tau2.start <- rep(tau2.start, length.out=N)
    }
    tau2.start <- matrix(tau2.start, N, 1)
    mvn.prior <- form.mvn.prior(g0, 0, N)
    g0 <- mvn.prior[[1]]
    c0 <- rep(c0, length.out=N)
    d0 <- rep(d0, length.out=N)
    check.ig.prior(a0, b0)
    for (i in 1:N){
      check.ig.prior(c0[i], d0[i])
    }
    w0 <- rep(w0, length.out=N)
    if (min(w0) < 0 | max(w0) > 1){
      cat("Element(s) of w0 not in [0, 1].\n")
      stop("Please respecify and call MCMCSVDreg again.\n")
    }

    ## define holder for posterior sample
    if (beta.samp){
      sample <- matrix(data=0, mcmc/thin, n.miss.Y + 2*N + 1 + K)
    }
    else{
      sample <- matrix(data=0, mcmc/thin, n.miss.Y + 2*N + 1)
    }


    ## call C++ code to draw sample
    posterior <- .C("cMCMCSVDreg",
                    sampledata = as.double(sample),
                    samplerow = as.integer(nrow(sample)),
                    samplecol = as.integer(ncol(sample)),
                    Y = as.double(Y),
                    Yrow = as.integer(nrow(Y)),
                    Ycol = as.integer(ncol(Y)),
                    Ymiss = as.integer(Y.miss.indic),
                    A = as.double(A),
                    Arow = as.integer(nrow(A)),
                    Acol = as.integer(ncol(A)),
                    D = as.double(D),
                    Drow = as.integer(nrow(D)),
                    Dcol = as.integer(ncol(D)),
                    F = as.double(F),
                    Frow = as.integer(nrow(F)),
                    Fcol = as.integer(ncol(F)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    tau2.start = as.double(tau2.start),
                    tau2row = as.integer(nrow(tau2.start)),
                    tau2col = as.integer(ncol(tau2.start)),
                    g0 = as.double(g0),
                    g0row = as.integer(nrow(g0)),
                    g0col = as.integer(ncol(g0)),
                    a0 = as.double(a0),
                    b0 = as.double(b0),
                    c0 = as.double(c0),
                    d0 = as.double(d0),
                    w0 = as.double(w0),
                    betasamp = as.integer(beta.samp),
                    PACKAGE="MCMCpack"
                    )

    ## pull together matrix and build MCMC object to return
    Y.miss.names <- NULL
    if (sum(Y.miss.indic) > 0){
      Y.miss.names <- paste("y", obsnames[Y.miss.indic==1], sep=".")
    }
    gamma.names <- paste("gamma", 1:N, sep=".")
    tau2.names <- paste("tau^2", 1:N, sep=".")
    beta.names <- paste("beta", xnames, sep=".")
    if (beta.samp){
      output <- form.mcmc.object(posterior,
                                 names=c(Y.miss.names, gamma.names, tau2.names,
                                   "sigma^2", beta.names),
                                 title="MCMCSVDreg Posterior Sample")
    }
    else{
      output <- form.mcmc.object(posterior,
                                 names=c(Y.miss.names, gamma.names, tau2.names,
                                   "sigma^2"),
                                 title="MCMCSVDreg Posterior Sample")
    }
    return(output)
  }
