\name{MCMCoprobit}
\alias{MCMCoprobit}
\title{Markov chain Monte Carlo for Ordered Probit Regression}
\description{
  This function generates a posterior density sample
  from an ordered probit regression model using the data augmentation 
  approach of Cowles (1996). The user supplies data and priors,
  and a sample from the posterior density is returned as an \code{mcmc}
  object, which can be subsequently analyzed with functions 
  provided in the coda package.
  }
  
\usage{
MCMCoprobit(formula, data = list(), burnin = 1000, mcmc = 10000,
   thin=5, tune = NA, verbose = FALSE, seed = 0, beta.start = NA,
   b0 = 0, B0 = 0.001, ...) }

\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value.}

    \item{tune}{The tuning parameter for the Metropolis-Hastings
      step. Default of NA corresponds to a choice of 0.05 divided by the
      number of categories in the response variable.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number and
    the betas are printed to the screen every 500 iterations.}

    \item{seed}{The seed for the random number generator.  The code uses the
    Mersenne Twister, which requires an integer as an input.  If nothing
    is provided, the Scythe default seed is used.}
  
    \item{beta.start}{The starting value for the \eqn{\beta}{beta} vector.   
    This can either be a scalar or a column vector with dimension equal
    to the number of betas. If this takes a scalar value, then that
    value will serve as the starting value for all of the betas. The
    default value of NA will use rescaled estimates from an ordered
    logit model.}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a column vector with dimension equal to the number of
    betas. If this takes a scalar value, then that value will serve as
    the prior mean for all of the betas.}
    
    \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either be a 
    scalar or a square matrix with dimensions equal to the number of
    betas.  If this takes a scalar value, then that value times an
    identity matrix serves as the prior precision of \eqn{\beta}{beta}.}
    
    \item{...}{further arguments to be passed}       
}

\value{
   An \code{mcmc} object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
\code{MCMCoprobit} simulates from the posterior density of a ordered probit
  regression model using data augmentation. The simulation proper is
  done in compiled C++ code to maximize efficiency.  Please consult the
  coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior density sample. 
  
  The observed variable \eqn{y_i}{y_i} is ordinal with a total of \eqn{C}{C} 
  categories, with distribution
  governed by a latent variable:
  \deqn{z_i = x_i'\beta + \varepsilon_i}{z_i = x_i'beta + epsilon_i}
  The errors are assumed to be from a standard Normal distribution.  The 
  probabilities of observing each outcome is governed by this latent
  variable and \eqn{C-1}{C-1} estimable cutpoints, which are denoted
  \eqn{\gamma_c}{gamma_c}.  The probability that individual \eqn{i}{i}
  is in category \eqn{c}{c} is computed by:
  \deqn{
   \pi_{ic} = \Phi(\gamma_c - x_i'\beta) - \Phi(\gamma_{c-1} - x_i'\beta)
   }{
   pi_ic = Phi(gamma_c - x_i'beta) - Phi(gamma_(c-1) - x_i'beta)
   }
   These probabilities are used to form the multinomial distribution
   that defines the likelihoods.
   
   The algorithm employed is discussed in depth by Cowles (1996).  Note that 
   the model does include a constant in the data matrix.  Thus, the first
   element  \eqn{\gamma_1}{gamma_1} is normalized to zero, and thus not 
   returned in the \code{mcmc} object.
}
  
\references{
  M. K. Cowles. 1996. ``Accelerating Monte Carlo Markov Chain Convergence for
  Cumulative-link Generalized Linear Models." \emph{Statistics and Computing.}
  6: 101-110.
     
  Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling." 
  Springer: New York.
      
  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2002.  
   \emph{Scythe Statistical Library 0.3.} \url{http://scythe.wustl.edu}.
   
  Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}
}


\examples{\dontrun{
x1 <- rnorm(100); x2 <- rnorm(100);
z <- 1.0 + x1*0.1 - x2*0.5 + rnorm(100);
y <- z; y[z < 0] <- 0; y[z >= 0 & z < 1] <- 1;
y[z >= 1 & z < 1.5] <- 2; y[z >= 1.5] <- 3;
posterior <- MCMCoprobit(y ~ x1 + x2, tune=0.3, mcmc=20000)
}}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}
