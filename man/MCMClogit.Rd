\name{MCMClogit}
\alias{MCMClogit}
\title{Markov chain Monte Carlo for Logistic Regression}
\description{
  This function generates a posterior density sample
  from a logistic regression model using a random walk Metropolis
  algorithm. The user supplies data and priors,
  and a sample from the posterior density is returned as an mcmc
  object, which can be subsequently analyzed with functions 
  provided in the coda package.
  }
  
\usage{
MCMClogit(formula, data = list(), burnin = 1000, mcmc = 10000,
   thin=5, tune=1.1, verbose = FALSE, seed = 0,  beta.start = NA,
   b0 = 0, B0 = 0.001, ...) }

\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Metropolis iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    mcmc iterations must be divisible by this value.}

    \item{tune}{Metropolis tuning parameter.  Make sure that the
    acceptance rate is satisfactory before using the posterior density sample
    for inference.}
    
    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number and
    the betas are printed to the screen every 500 iterations.}

    \item{seed}{The seed for the random number generator.  The code uses the
    Mersenne Twister, which requires an integer as an input.  If nothing
    is provided, the Scythe default seed is used.}

    \item{beta.start}{The starting value for the \eqn{\beta}{beta} vector.
    This can either 
    be a scalar or a column vector with dimension equal to the number of 
    betas. If this takes a scalar value, then that value will serve as the 
    starting value for all of the betas.  The default value of NA will
    use the maximum likelihood estimate of \eqn{\beta}{beta} as the starting 
    value.}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a column      
    vector with dimension equal to the number of betas. If this takes a scalar
    value, then that value will serve as the prior mean for all of the
    betas.}
    
    \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either be a
    scalar
    or a square matrix with dimensions equal to the number of betas.  If this
    takes a scalar value, then that value times an identity matrix serves
    as the prior precision of \eqn{\beta}{beta}.}
    
    \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{\code{MCMClogit} simulates from the posterior density of a logistic
  regression model using a random walk Metropolis algorithm. The simulation
  proper is done in compiled C++ code to maximize efficiency.  Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior density sample.
  
  The model takes the following form:
  \deqn{y_i \sim \mathcal{B}ernoulli(\pi_i)}{y_i ~ Bernoulli(pi_i)}
  Where the inverse link function:
  \deqn{\pi_i = \frac{\exp(x_i'\beta)}{1 + \exp(x_i'\beta)}}{pi_i = 
    exp(x_i'beta) / [1 + exp(x_i'beta)]}
  We assume a multivariate Normal prior on \eqn{\beta}{beta}:
    \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}{beta ~ N(b0,B0^(-1))}
  The candidate generating density is a multivariate Normal density centered
  at the current value of \eqn{\beta}{beta} with variance-covariance matrix
  that is an approximation of the posterior based on the maximum likelihood 
  estimates and the prior precision multiplied by the tuning parameter  
   squared.
  }
  
\references{
      
  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2003.  
   \emph{Scythe Statistical Library 0.4.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}


\examples{
   \dontrun{
   data(birthwt)
   posterior <- MCMClogit(low~age+as.factor(race)+smoke, data=birthwt)
   plot(posterior)
   summary(posterior)
   }
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}, \code{\link[base]{glm}}}

