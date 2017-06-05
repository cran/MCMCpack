##########################################################################
## simple instructional models using Monte Carlo simulation
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

## Monte Carlo simulation from the likelihood of a
## binomial distribution with a Beta(alpha, beta) prior
## ADM 1/25/2006

#' Monte Carlo Simulation from a Binomial Likelihood with a Beta Prior
#'
#' This function generates a sample from the posterior distribution of a
#' binomial likelihood with a Beta prior.
#'
#' \code{MCbinomialbeta} directly simulates from the posterior distribution.
#' This model is designed primarily for instructional use.  \eqn{\pi} is
#' the probability of success for each independent Bernoulli trial.  We assume
#' a conjugate Beta prior:
#'
#' \deqn{\pi \sim \mathcal{B}eta(\alpha, \beta)}
#'
#' \eqn{y} is the number of successes in \eqn{n} trials.  By
#' default, a uniform prior is used.
#'
#' @param y The number of successes in the independent Bernoulli trials.
#'
#' @param n The number of independent Bernoulli trials.
#'
#' @param alpha Beta prior distribution alpha parameter.
#'
#' @param beta Beta prior distribution beta parameter.
#'
#' @param mc The number of Monte Carlo draws to make.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' posterior <- MCbinomialbeta(3,12,mc=5000)
#' summary(posterior)
#' plot(posterior)
#' grid <- seq(0,1,0.01)
#' plot(grid, dbeta(grid, 1, 1), type="l", col="red", lwd=3, ylim=c(0,3.6),
#'   xlab="pi", ylab="density")
#' lines(density(posterior), col="blue", lwd=3)
#' legend(.75, 3.6, c("prior", "posterior"), lwd=3, col=c("red", "blue"))
#' }
#'
MCbinomialbeta <- function(y, n, alpha=1, beta=1, mc=1000, ...) {

   # check data
   if (any(y < 0)) {
        cat("Error: Number of successes negative.\n")
        stop("Please respecify and call function again.")
    }
    if (any(n < 0)) {
        cat("Error: Number of trials negative.\n")
        stop("Please respecify and call function again.")
    }
    if (any(y > n)) {
        cat("Error: Number of successes greater than number of trials.\n")
        stop("Please respecify and call function again.")
    }

   # check other parameters
   check.beta.prior(alpha, beta)
   check.mc.parameter(mc)

   # draw sample and return
   output <- mcmc(matrix(rbeta(mc, alpha+y, beta+n-y),mc,1))
   varnames(output) <- as.list("pi")
   attr(output,"title") <- "MCbinomialbeta Posterior Sample"

   return(output)
}

## Monte Carlo simulation from the likelihood of a
## Poisson distribution with a Gamma(alpha, beta) prior
## ADM 1/25/2006
#' Monte Carlo Simulation from a Poisson Likelihood with a Gamma Prior
#'
#' This function generates a sample from the posterior distribution of a
#' Poisson likelihood with a Gamma prior.
#'
#' \code{MCpoissongamma} directly simulates from the posterior distribution.
#' This model is designed primarily for instructional use.
#' \eqn{\lambda} is the parameter of interest of the Poisson
#' distribution.  We assume a conjugate Gamma prior:
#'
#' \deqn{\lambda \sim \mathcal{G}amma(\alpha, \beta)}
#'
#' \eqn{y} is a vector of counts.
#'
#' @param y A vector of counts (must be non-negative).
#'
#' @param alpha Gamma prior distribution shape parameter.
#'
#' @param beta Gamma prior distribution scale parameter.
#'
#' @param mc The number of Monte Carlo draws to make.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' data(quine)
#' posterior <- MCpoissongamma(quine$Days, 15, 1, 5000)
#' summary(posterior)
#' plot(posterior)
#' grid <- seq(14,18,0.01)
#' plot(grid, dgamma(grid, 15, 1), type="l", col="red", lwd=3, ylim=c(0,1.3),
#'   xlab="lambda", ylab="density")
#' lines(density(posterior), col="blue", lwd=3)
#' legend(17, 1.3, c("prior", "posterior"), lwd=3, col=c("red", "blue"))
#' }
#'
MCpoissongamma <- function(y, alpha, beta, mc=1000, ...) {

   # check data
   if(any(y < 0)) {
      cat("Error: Some counts negative in y.\n")
      stop("Please respecify and call function again.")
    }
   n <- length(y)

   # check other parameters
   check.gamma.prior(alpha, beta)
   check.mc.parameter(mc)

   # draw sample and return
   output <- mcmc(matrix(rgamma(mc, alpha+sum(y), beta+n),mc,1))
   varnames(output) <- as.list("lambda")
   attr(output,"title") <- "MCpoissongamma Posterior Sample"

   return(output)
}

## Monte Carlo simulation from the likelihood of a
## Normal distribution with a Normal(mu0, tau20) prior
## the variance sigma2 is known
## ADM 1/26/2006
#' Monte Carlo Simulation from a Normal Likelihood (with known variance) with a
#' Normal Prior
#'
#' This function generates a sample from the posterior distribution of a Normal
#' likelihood (with known variance) with a Normal prior.
#'
#' \code{MCnormalnormal} directly simulates from the posterior distribution.
#' This model is designed primarily for instructional use.  \eqn{\mu} is
#' the parameter of interest of the Normal distribution.  We assume a conjugate
#' normal prior:
#'
#' \deqn{\mu \sim \mathcal{N}(\mu_0, \tau^2_0)}
#'
#' \eqn{y} is a vector of observed data.
#'
#' @param y The data.
#'
#' @param sigma2 The known variance of y.
#'
#' @param mu0 The prior mean of mu.
#'
#' @param tau20 The prior variance of mu.
#'
#' @param mc The number of Monte Carlo draws to make.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' y <- c(2.65, 1.80, 2.29, 2.11, 2.27, 2.61, 2.49, 0.96, 1.72, 2.40)
#' posterior <- MCMCpack:::MCnormalnormal(y, 1, 0, 1, 5000)
#' summary(posterior)
#' plot(posterior)
#' grid <- seq(-3,3,0.01)
#' plot(grid, dnorm(grid, 0, 1), type="l", col="red", lwd=3, ylim=c(0,1.4),
#'    xlab="mu", ylab="density")
#' lines(density(posterior), col="blue", lwd=3)
#' legend(-3, 1.4, c("prior", "posterior"), lwd=3, col=c("red", "blue"))
#' }
#'
MCnormalnormal <- function(y, sigma2, mu0, tau20, mc=1000, ...) {

   n <- length(y)
   if(sigma2 <= 0) {
      cat("Error: Known variance sigma2 is less than or equal to zero.\n")
      stop("Please respecify and call function again.")
    }

   # check other parameters
   check.normal.prior(mu0, tau20)
   check.mc.parameter(mc)

   # draw sample and return
   mu1 = (1/tau20 * mu0 + n/sigma2 * mean(y)) / (1/tau20 + n/sigma2)
   tau21 = 1 / (1/tau20 + n/sigma2)
   output <- mcmc(matrix(rnorm(mc, mu1, sqrt(tau21)),mc,1))
   varnames(output) <- as.list("mu")
   attr(output,"title") <- "MCnormalnormal Posterior Sample"

   return(output)
}

## Monte Carlo simulation from the likelihood of a
## multinomal distribution with a Dirichlet(alpha) prior
#' Monte Carlo Simulation from a Multinomial Likelihood with a Dirichlet Prior
#'
#' This function generates a sample from the posterior distribution of a
#' multinomial likelihood with a Dirichlet prior.
#'
#' \code{MCmultinomdirichlet} directly simulates from the posterior
#' distribution.  This model is designed primarily for instructional use.
#' \eqn{\pi} is the parameter of interest of the multinomial distribution.
#' It is of dimension \eqn{(d \times 1)}. We assume a conjugate
#' Dirichlet prior:
#'
#' \deqn{\pi \sim \mathcal{D}irichlet(\alpha_0)}
#'
#' \eqn{y} is a \eqn{(d \times 1)} vector of
#' observed data.
#'
#' @param y A vector of data (number of successes for each category).
#'
#' @param alpha0 The vector of parameters of the Dirichlet prior.
#'
#' @param mc The number of Monte Carlo draws to make.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' ## Example from Gelman, et. al. (1995, p. 78)
#' posterior <- MCmultinomdirichlet(c(727,583,137), c(1,1,1), mc=10000)
#' bush.dukakis.diff <- posterior[,1] - posterior[,2]
#' cat("Pr(Bush > Dukakis): ",
#'    sum(bush.dukakis.diff > 0) / length(bush.dukakis.diff), "\n")
#' hist(bush.dukakis.diff)
#' }
#'
MCmultinomdirichlet <- function(y, alpha0, mc=1000, ...) {

   # check data
   d <- length(y)
   if(any(y < 0)) {
      cat("Error: Some counts negative in y.\n")
      stop("Please respecify and call function again.")
    }

   # check alpha
   if(length(alpha0) != d) {
      cat("Error: Dimension of alpha and y do not match.\n")
      stop("Please respecify and call function again.")
   }
   if(any(alpha0 <= 0)) {
      cat("Error: At least one alpha in Dirichlet prior less than or equal to zero.\n")
      stop("Please respecify and call function again.")
   }

   # draw sample and return
   output <- mcmc(rdirichlet(mc,y + alpha0))
   varnames(output) <- paste("pi.", 1:d, sep="")
   attr(output,"title") <- "MCmultinomdirichlet Posterior Sample"

   return(output)
}
