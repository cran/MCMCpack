\name{MCMCtobit}
\alias{MCMCtobit}
\title{Markov Chain Monte Carlo for Gaussian Linear Regression with a Censored Dependent Variable}
\description{
  This function generates a posterior density sample
  from a linear regression model with Gaussian errors using
  Gibbs sampling (with a multivariate Gaussian prior on the
  beta vector, and an inverse Gamma prior on the conditional
  error variance).  The dependent variable may be censored
  from below, from above, or both. The user supplies data 
  and priors, and a sample from the posterior density is 
  returned as an mcmc object, which can be subsequently 
  analyzed with functions provided in the coda package.
  }
\usage{
MCMCtobit(formula, data = parent.frame(), below = 0, above = Inf, 
  burnin = 1000, mcmc = 10000, thin = 1, verbose = FALSE, seed = NA, 
  beta.start = NA, b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, ...)
}
\arguments{
    \item{formula}{A model formula.}
  
    \item{data}{A dataframe.}
  
    \item{below}{The point at which the dependent variable is 
      censored from below. The default is zero. To censor from 
      above only, specify that below = -Inf.}
    
    \item{above}{The point at which the dependent variable is 
      censored from above. To censor from below only, use the
      default value of Inf.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burnin.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number, the
    \eqn{\beta}{beta} vector, and the conditional error variance is printed to 
    the screen every 1000 iterations.}

    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack specification 
    for more details.}

    \item{beta.start}{The starting values for the \eqn{\beta}{beta} vector.
    This can either be a scalar or a column vector with dimension equal to the 
    number of betas. The default value of of NA will use the OLS estimate of 
    \eqn{\beta}{beta} as the starting value.  If this is a scalar, that value 
    will serve as the starting value mean for all of the betas.}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a column vector with dimension equal to the number of betas. If 
    this takes a scalar  value, then that value will serve as the prior
    mean for all of the betas.} 
    
    \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either be a 
    scalar or a square matrix with dimensions equal to the number of
    betas.
    If this takes a scalar value, then that value times an
    identity matrix serves as the prior precision of beta. Default value
    of 0 is equivalent to an improper uniform prior for beta.}
    
    \item{c0}{\eqn{c_0/2}{c0/2} is the shape parameter for the inverse
    Gamma prior on \eqn{\sigma^2}{sigma^2} (the variance of the
    disturbances). The amount of information in the inverse Gamma prior
    is something like that from \eqn{c_0}{c0} pseudo-observations.} 
    
    \item{d0}{\eqn{d_0/2}{d0/2} is the scale parameter for the
    inverse Gamma prior on \eqn{\sigma^2}{sigma^2} (the variance of the
    disturbances). In constructing the inverse Gamma prior,
    \eqn{d_0}{d0} acts like the sum of squared errors from the
    \eqn{c_0}{c0} pseudo-observations.}
    
    \item{...}{further arguments to be passed}       
    
}

\details{
  \code{MCMCtobit} simulates from the posterior density using standard
  Gibbs sampling (a multivariate Normal draw for the betas, and an
  inverse Gamma draw for the conditional error variance). \code{MCMCtobit}
  differs from \code{MCMCregress} in that the dependent variable may be
  censored from below, from above, or both. The simulation
  proper is done in compiled C++ code to maximize efficiency.  Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior density sample.
  
  The model takes the following form:
  \deqn{y_i = x_i ' \beta + \varepsilon_{i},}{y_i = x_i'beta + epsilon_i,}
  where the errors are assumed to be Gaussian:
  \deqn{\varepsilon_{i} \sim \mathcal{N}(0, \sigma^2).}{epsilon_i ~ N(0,
    sigma^2).}
  Let \eqn{c_1} and \eqn{c_2} be the two censoring points, and let
  \eqn{y_i^\ast}{y_i^star} be the partially observed dependent variable. Then,
  \deqn{y_i = y_i^{\ast} \texttt{ if } c_1 < y_i^{\ast} < c_2,}{y_i = y_i^star if c_1 < y_i^star < c_2,}
  \deqn{y_i = c_1 \texttt{ if } c_1 \geq y_i^{\ast},}{y_i = c_1 if c_1 >= y_i^star,}
  \deqn{y_i = c_2 \texttt{ if } c_2 \leq y_i^{\ast}.}{y_i = c_2 if c_1 <= y_i^star.}
  We assume standard, semi-conjugate priors:
  \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1}),}{beta ~ N(b0,B0^(-1)),}
  and:
  \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2),}{sigma^(-2) ~
    Gamma(c0/2, d0/2),}
  where \eqn{\beta}{beta} and \eqn{\sigma^{-2}}{sigma^(-2)} are assumed 
  \emph{a priori} independent.  Note that only starting values for
  \eqn{\beta}{beta} are allowed because simulation is done using
  Gibbs sampling with the conditional error variance
  as the first block in the sampler.
  }

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\references{
   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  
   \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
   
   James Tobin. 1958. ``Estimation of relationships for limited dependent 
   variables." \emph{Econometrica.}, 26:24-36.
}

\author{Ben Goodrich}

\seealso{
  \code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}, 
  \code{\link[survival]{survreg}},
  \code{\link[MCMCpack]{MCMCregress}}}

\examples{
library(survival)
example(tobin)
summary(tfit)
tfit.mcmc <- MCMCtobit(durable ~ age + quant, data=tobin, mcmc=30000, verbose=TRUE)
plot(tfit.mcmc)
raftery.diag(tfit.mcmc)
summary(tfit.mcmc)
}

\keyword{models}
