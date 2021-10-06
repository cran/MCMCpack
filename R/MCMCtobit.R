##########################################################################
## tobit regression model
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


#' Markov Chain Monte Carlo for Gaussian Linear Regression with a Censored
#' Dependent Variable
#'
#' This function generates a sample from the posterior distribution of a linear
#' regression model with Gaussian errors using Gibbs sampling (with a
#' multivariate Gaussian prior on the beta vector, and an inverse Gamma prior
#' on the conditional error variance).  The dependent variable may be censored
#' from below, from above, or both. The user supplies data and priors, and a
#' sample from the posterior distribution is returned as an mcmc object, which
#' can be subsequently analyzed with functions provided in the coda package.
#'
#' \code{MCMCtobit} simulates from the posterior distribution using standard
#' Gibbs sampling (a multivariate Normal draw for the betas, and an inverse
#' Gamma draw for the conditional error variance). \code{MCMCtobit} differs
#' from \code{MCMCregress} in that the dependent variable may be censored from
#' below, from above, or both. The simulation proper is done in compiled C++
#' code to maximize efficiency.  Please consult the coda documentation for a
#' comprehensive list of functions that can be used to analyze the posterior
#' sample.
#'
#' The model takes the following form:
#'
#' \deqn{y_i = x_i ' \beta + \varepsilon_{i},}
#'
#' where the errors are assumed
#' to be Gaussian:
#'
#' \deqn{\varepsilon_{i} \sim \mathcal{N}(0, \sigma^2).}
#'
#' Let \eqn{c_1} and \eqn{c_2} be the two censoring points, and let
#' \eqn{y_i^\ast} be the partially observed dependent variable. Then,
#'
#' \deqn{y_i = y_i^{\ast} \texttt{ if } c_1 < y_i^{\ast} < c_2,}
#'
#' \deqn{y_i = c_1 \texttt{ if } c_1 \geq y_i^{\ast},}
#'
#' \deqn{y_i = c_2 \texttt{ if } c_2 \leq y_i^{\ast}.}
#'
#' We assume standard, semi-conjugate priors:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1}),}
#'
#' and:
#'
#' \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2),}
#'
#' where \eqn{\beta} and \eqn{\sigma^{-2}} are
#' assumed \emph{a priori} independent.  Note that only starting
#' values for \eqn{\beta} are allowed because simulation is done
#' using Gibbs sampling with the conditional error variance as the
#' first block in the sampler.
#'
#' @param formula A model formula.
#'
#' @param data A dataframe.
#'
#' @param below The point at which the dependent variable is censored from
#' below. The default is zero. To censor from above only, specify that below =
#' -Inf.
#' @param above The point at which the dependent variable is censored from
#' above. To censor from below only, use the default value of Inf.
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
#' iteration number, the \eqn{\beta} vector, and the error variance is
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
#' number of betas. The default value of of NA will use the OLS estimate of
#' \eqn{\beta} as the starting value.  If this is a scalar, that value
#' will serve as the starting value mean for all of the betas.
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
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @author Ben Goodrich, \email{goodrich.ben@gmail.com},
#' \url{http://www.columbia.edu/~bg2382/}
#'
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}},
#' \code{\link[survival]{survreg}}, \code{\link[MCMCpack]{MCMCregress}}
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
#' Siddhartha Chib. 1992. ``Bayes inference in the Tobit censored regression
#' model."  \emph{Journal of Econometrics}. 51:79-99.
#'
#' James Tobin. 1958. ``Estimation of relationships for limited dependent
#' variables." \emph{Econometrica.} 26:24-36.
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' library(survival)
#' example(tobin)
#' summary(tfit)
#' tfit.mcmc <- MCMCtobit(durable ~ age + quant, data=tobin, mcmc=30000,
#'                         verbose=1000)
#' plot(tfit.mcmc)
#' raftery.diag(tfit.mcmc)
#' summary(tfit.mcmc)
#' }
#'
"MCMCtobit" <-
function(formula, data=NULL, below = 0.0, above = Inf,
           burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, ...) {

    # checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    if (!is.numeric(below) | !is.numeric(above)) {
        cat("Error: Censoring points must be numeric, which includes +-Inf.\n")
        stop("Please respecify and call ", calling.function(), " again.",
        call.=FALSE)
    }
    if (below >= above) {
        cat("Error: Truncation points are logically inconsistent.\n")
        stop("Please respecify and call ", calling.function(), " again.",
        call.=FALSE)
    }

    # convert infinite values to finite approximations
    if(is.infinite(below)) below <- .Machine$double.xmax*-1
    if(is.infinite(above)) above <- .Machine$double.xmax

    # seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    # form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    K <- ncol(X)  # number of covariates

    # starting values and priors
    beta.start <- coef.start(beta.start, K, formula, family=gaussian, data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    check.ig.prior(c0, d0)

    # define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, K+1)
    posterior <- NULL

    # call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="cMCMCtobit",
                     sample.nonconst=sample, Y=Y, X=X, below=as.double(below),
                     above=as.double(above),
                     burnin=as.integer(burnin), mcmc=as.integer(mcmc),
                     thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, c0=as.double(c0), d0=as.double(d0))

    # pull together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior,
                               names=c(xnames, "sigma2"),
                               title="MCMCtobit Posterior Sample")
    return(output)
  }
