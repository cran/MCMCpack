##########################################################################
## sample from the posterior distribution of an ordered probit model
## via the data augmentation approach of Cowles (1996)
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 1/25/2003
## Modified to meet new developer specification 7/26/2004 KQ
## Modified for new Scythe and rngs 7/26/2004 KQ
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for Ordered Probit Regression
#'
#' This function generates a sample from the posterior distribution of
#' an ordered probit regression model using the data augmentation approach of
#' Albert and Chib (1993), with cut-points sampled according to Cowles (1996)
#' or Albert and Chib (2001). The user supplies data and priors, and a sample from the
#' posterior distribution is returned as an mcmc object, which can be
#' subsequently analyzed with functions provided in the coda package.
#'
#' \code{MCMCoprobit} simulates from the posterior distribution of a ordered
#' probit regression model using data augmentation. The simulation proper is
#' done in compiled C++ code to maximize efficiency.  Please consult the coda
#' documentation for a comprehensive list of functions that can be used to
#' analyze the posterior sample.
#'
#' The observed variable \eqn{y_i} is ordinal with a total of \eqn{C}
#' categories, with distribution governed by a latent variable: \deqn{z_i =
#' x_i'\beta + \varepsilon_i} The errors are
#' assumed to be from a standard Normal distribution.  The probabilities of
#' observing each outcome is governed by this latent variable and
#' \eqn{C-1} estimable cutpoints, which are denoted
#' \eqn{\gamma_c}.  The probability that individual \eqn{i} is in
#' category \eqn{c} is computed by:
#'
#' \deqn{\pi_{ic} = \Phi(\gamma_c - x_i'\beta) - \Phi(\gamma_{c-1} -
#' x_i'\beta)}
#'
#' These probabilities are used to form the multinomial distribution
#' that defines the likelihoods.
#'
#' \code{MCMCoprobit} provides two ways to sample the cutpoints. Cowles (1996)
#' proposes a sampling scheme that groups sampling of a latent variable with
#' cutpoints.  In this case, for identification the first element
#' \eqn{\gamma_1} is normalized to zero. Albert and Chib (2001) show
#' that we can sample cutpoints indirectly without constraints by transforming
#' cutpoints into real-valued parameters (\eqn{\alpha}).
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' Gibbs iterations must be divisible by this value.
#'
#' @param tune The tuning parameter for the Metropolis-Hastings step. Default
#' of NA corresponds to a choice of 0.05 divided by the number of categories in
#' the response variable.
#'
#' @param tdf Degrees of freedom for the multivariate-t proposal distribution
#' when \code{mcmc.method} is set to "IndMH". Must be positive.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the beta vector, and the Metropolis-Hastings acceptance
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
#' starting value for all of the betas. The default value of NA will use
#' rescaled estimates from an ordered logit model.
#'
#' @param b0 The prior mean of \eqn{\beta}.  This can either be a scalar
#' or a column vector with dimension equal to the number of betas. If this
#' takes a scalar value, then that value will serve as the prior mean for all
#' of the betas.
#'
#' @param B0 The prior precision of \eqn{\beta}.  This can either be a
#' scalar or a square matrix with dimensions equal to the number of betas.  If
#' this takes a scalar value, then that value times an identity matrix serves
#' as the prior precision of \eqn{\beta}.  Default value of 0 is
#' equivalent to an improper uniform prior on \eqn{\beta}.
#'
#' @param a0 The prior mean of \eqn{\gamma}.  This can either be a
#' scalar or a column vector with dimension equal to the number of betas. If
#' this takes a scalar value, then that value will serve as the prior mean for
#' all of the betas.
#'
#' @param A0 The prior precision of \eqn{\gamma}.  This can either be a
#' scalar or a square matrix with dimensions equal to the number of betas.  If
#' this takes a scalar value, then that value times an identity matrix serves
#' as the prior precision of \eqn{\gamma}.  Default value of 0 is
#' equivalent to an improper uniform prior on \eqn{\gamma}.
#'
#' @param mcmc.method Can be set to either "Cowles" (default) or "AC" to
#' perform posterior sampling of cutpoints based on Cowles (1996) or Albert and
#' Chib (2001) respectively.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}
#'
#' @references Albert, J. H. and S. Chib. 1993. ``Bayesian Analysis of Binary
#' and Polychotomous Response Data.'' \emph{J. Amer. Statist. Assoc.} 88,
#' 669-679
#'
#' M. K. Cowles. 1996. ``Accelerating Monte Carlo Markov Chain Convergence for
#' Cumulative-link Generalized Linear Models." \emph{Statistics and Computing.}
#' 6: 101-110.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Valen E. Johnson and James H. Albert. 1999. \emph{Ordinal Data Modeling}.
#' Springer: New York.
#'
#' Albert, James and Siddhartha Chib. 2001. ``Sequential Ordinal Modeling with
#' Applications to Survival Data." \emph{Biometrics.} 57: 829-836.
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
#'    x1 <- rnorm(100); x2 <- rnorm(100);
#'    z <- 1.0 + x1*0.1 - x2*0.5 + rnorm(100);
#'    y <- z; y[z < 0] <- 0; y[z >= 0 & z < 1] <- 1;
#'    y[z >= 1 & z < 1.5] <- 2; y[z >= 1.5] <- 3;
#'    out1 <- MCMCoprobit(y ~ x1 + x2, tune=0.3)
#'    out2 <- MCMCoprobit(y ~ x1 + x2, tune=0.3, tdf=3, verbose=1000, mcmc.method="AC")
#'    summary(out1)
#'    summary(out2)
#'    plot(out1)
#'    plot(out2)
#'    }
#'
"MCMCoprobit" <-
  function(formula, data = parent.frame(), burnin = 1000, mcmc = 10000,
           thin = 1, tune = NA, tdf = 1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, a0 = 0, A0 = 0, mcmc.method = c("Cowles", "AC"), ...) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## extract X, Y, and variable names from the model formula and frame
    call <- match.call()
    mt <- terms(formula, data=data)
    if(missing(data)) data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$burnin <- mf$mcmc <- mf$b0 <- mf$B0 <- mf$a0 <- mf$A0 <- NULL
    mf$thin <- mf$... <- mf$tune <- mf$tdf <- mf$verbose <- mf$seed <- NULL
    mf$beta.start <- mf$mcmc.method <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    vars <- as.character(attr(mt, "variables"))[-1] # y varname and x varnames

    ## null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf)# else NULL
    X.names <- dimnames(X)[[2]]
    Y    <- model.response(mf, "numeric")
    Y    <- factor(Y, ordered=TRUE)
    ncat <- nlevels(Y)             # number of categories of y
    cat  <- levels(Y)              # values of categories of y
    N <- nrow(X)	               # number of observations
    K <- ncol(X)	               # number of covariates
    if (length(Y) != N){
      cat("X and Y do not have same number of rows.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }

    ## convert data to matrices to be passed
    Y <- as.matrix(as.integer(Y))
    X <- as.matrix(X)

    ## check tuning parameter
    if (is.na(tune)){
      tune <- 0.05/ncat
    }

    xint <- match("(Intercept)", colnames(X), nomatch=0)
    if (xint > 0){
      new.X <- X[, -xint, drop=FALSE]
    }
    else warning("An intercept is needed and assumed in MCMCoprobit()\n.")
    if (ncol(new.X) == 0){
      polr.out <- polr(ordered(Y)~1)
    }
    else {
      polr.out <- polr(ordered(Y)~new.X)
    }

    ## starting values for beta error checking
    if (is.na(beta.start)){
      beta.start <- matrix(0, K, 1)
      beta.start[1] <- -.588 * polr.out$zeta[1]
      if( ncol(new.X) > 0){
        beta.start[2:K] <- .588 * coef(polr.out)
      }
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,K,1)
    }
    else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
      cat("Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }

    ## prior for beta error checking
    if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,K,1)
    }
    if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
      cat("N(b0,B0) prior b0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }
    if(is.null(dim(B0))) {
      B0 <- B0 * diag(K)
    }
    if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
      cat("N(b0,B0) prior B0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }

    ## prior for alpha error checking
    if(is.null(dim(a0))) {
      a0 <- a0 * matrix(1, ncat-1, 1)
    }
    if((dim(a0)[1] != ncat-1) || (dim(a0)[2] != 1)) {
      cat("N(a0,A0) prior a0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }
    if(is.null(dim(A0))) {
      A0 <- A0 + diag(ncat - 1)
    }
    if((dim(A0)[1] != ncat - 1) || (dim(A0)[2] != ncat - 1)) {
      cat("N(a0, A0) prior A0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }

    ## form gamma starting values (note: not changeable)
    gamma <- matrix(NA,ncat+1,1)
    gamma[1] <- -300
    gamma[2] <- 0
    gamma[3:ncat] <- (polr.out$zeta[2:(ncat-1)] - polr.out$zeta[1])*.588
    gamma[ncat+1] <- 300

    ## posterior sample
    sample <- matrix(data=0, mcmc/thin, K + ncat + 1)

    ## call C++ code to draw sample
    nY <- as.matrix(as.numeric(Y))

    ## mcmc.method
    cowles <- as.integer(1)
    if(mcmc.method[1]=="AC") {cowles <- as.integer(0)}

    ## form the tuning parameter
    tune <- vector.tune(tune, ncat-1)
    posterior <- NULL

    auto.Scythe.call(output.object="posterior", cc.fun.name="cMCMCoprobit",
                     sample.nonconst=sample, Y=as.integer(Y), nY=nY, X=X,
                     burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     tune=tune, tdf=as.double(tdf),
                     lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), beta=beta.start,
                     gamma=gamma, b0=b0, B0=B0, a0=a0, A0=A0,
                     cowles=as.integer(cowles))

    ## put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol, byrow=FALSE)
    if(mcmc.method[1]=="AC"){
      sample[ , 1] <- sample[, 1] - sample[, K+2] ## post-MCMC normalization
      sample[ , (K+2):(K+ncat)] <- sample[ , (K+2):(K+ncat)] - sample[, K+2] ## post-MCMC normalization
    }
    sample <- sample[, c(1:K, (K+3):(K+ncat))]

    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    xnames <- c(X.names, paste("gamma", 2:(ncat-1), sep=""))
    varnames(output) <- xnames
    attr(output, "title") <- "MCMCoprobit Posterior Sample"

    return(output)
  }
