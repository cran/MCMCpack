\name{MCMCregress}
\alias{MCMCregress}
\title{Markov chain Monte Carlo for Gaussian Linear Regression}
\description{
  This function generates a posterior density sample
  from a linear regression model with Gaussian errors using
  Gibbs sampling (with a multivariate Gaussian prior on the
  beta vector, and an inverse-Gamma prior on the conditional
  error variance).  The user supplies data and priors, and 
  a sample from the posterior density is returned as an \code{mcmc}
  object, which can be subsequently analyzed with functions 
  provided in the coda package.
  }
  
\usage{
MCMCregress(formula, data = list(), burnin = 1000, mcmc = 10000,
   thin = 5, verbose = FALSE, seed = 0, beta.start = NA, sigma2.start = NA,
   b0 = 0, B0 = 0, nu = 0.001, delta = 0.001, ...) }

\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Gibbs iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number, the
    \eqn{\beta}{beta} vector, and the conditional error variance is printed to 
    the screen 
    every 500 iterations.}

    \item{seed}{The seed for the random number generator.  The code uses the
    Mersenne Twister, which requires an integer as an input.  If nothing
    is provided, the Scythe default seed is used.}

    \item{beta.start}{The starting value for the \eqn{\beta}{beta} vector.  
    This can either 
    be a scalar or a column vector with dimension equal to the number of 
    betas. If this takes a scalar value, then that value will serve as the 
    starting value for all of the betas. The default value of NA will
    use the maximum likelihood estimate of \eqn{\beta}{beta} as the starting 
    value.}

    \item{sigma2.start}{The starting value for the conditional error
      variance. The default value of of NA will use the maximum
      likelihood estimate of \eqn{\sigma^2}{sigma2} as the starting value.}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
     scalar or a
      column vector with dimension equal to the number of betas. If this
      takes a scalar  value, then that value will serve as the prior
      mean for all of the betas.} 
    
    \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either be a 
    scalar or a square matrix with dimensions equal to the number of betas.  
    If this
    takes a scalar value, then that value times an identity matrix serves
    as the prior precision of beta. Default value of 0 is equivalent to
    an improper uniform prior for beta.}
    
    \item{nu}{\eqn{\nu/2}{nu/2} is the shape parameter for inverse-Gamma prior 
     on the conditional error  variance.}
    
    \item{delta}{\eqn{\delta/2}{delta/2} is the scale parameter for 
    inverse-Gamma prior on the conditional error variance.}
    
    \item{...}{further arguments to be passed}       
}

\value{
   An \code{mcmc} object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCMCregress} simulates from the posterior density using 
  standard Gibbs sampling (a multivariate Normal draw for the betas, and an
  inverse-Gamma draw for the conditional error variance).  The simulation
  proper is done in compiled C++ code to maximize efficiency.  Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior density sample.
  
  The model takes the following form:
  \deqn{y_i = x_i ' \beta + \varepsilon_{i}}{y_i = x_i'beta + epsilon_i}
  Where the errors are assumed to be Gaussian:
  \deqn{\varepsilon_{i} \sim \mathcal{N}(0, \sigma^2)}{epsilon_i ~ N(0,
    sigma^2)} 
  We assume standard, conjugate priors:
  \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}{beta ~ N(b0,B0^(-1))}
  And:
  \deqn{\sigma^{-2} \sim \mathcal{G}amma(\nu/2, \delta/2)}{sigma^(-2) ~
    Gamma(nu/2, delta/2)}
  Where \eqn{\beta}{beta} and \eqn{\sigma^{-2}}{sigma^(-2)} are assumed 
  \emph{a priori} independent.
  }
  
\references{
   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2002.  
   \emph{Scythe Statistical Library 0.3.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}


\examples{
line   <- list(X = c(1,2,3,4,5), Y = c(1,3,3,3,5))
line$X <- line$X - mean(line$X)
post2  <- MCMCregress(Y~X, data=line)
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}, \code{\link[base]{lm}}}

