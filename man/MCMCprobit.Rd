\name{MCMCprobit}
\alias{MCMCprobit}
\title{Markov chain Monte Carlo for Probit Regression}
\description{
  This function generates a posterior density sample
  from a probit regression model using the data augmentation
  approach of Albert and Chib (1993). The user supplies data and priors,
  and a sample from the posterior density is returned as an \code{mcmc}
  object, which can be subsequently analyzed with functions 
  provided in the coda package.
  }
  
\usage{
MCMCprobit(formula, data = list(), burnin = 1000, mcmc = 10000,
   thin = 5, verbose = FALSE, seed = 0, beta.start = NA,
   b0 = 0, B0 = 0, bayes.resid = FALSE, ...) }

\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Gibbs iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number and
    the betas are printed to the screen.}

    \item{seed}{The seed for the random number generator.  The code uses the
    Mersenne Twister, which requires an integer as an input.  If nothing
    is provided, the Scythe default seed is used.}
  
    \item{beta.start}{The starting value for the \eqn{\beta}{beta} vector.    
    This can either 
    be a scalar or a column vector with dimension equal to the number of 
    betas.  If  If this takes a scalar value, then that value will serve as 
    the 
    starting value for all of the betas. The default value of NA will
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
    as the prior precision of \eqn{\beta}{beta}. Default value of 0 is 
    equivalent to
    an improper uniform prior on \eqn{\beta}{beta} .}

   \item{bayes.resid}{Should latent Bayesian residuals (Albert and Chib,
    1995) be returned? Default is FALSE meaning no residuals should be
    returned. Alternatively, the user can specify an array of integers
    giving the observation numbers for which latent residuals should be
    calculated and returned. TRUE will return draws of
    latent residuals for all observations.}
    
    \item{...}{further arguments to be passed}       
}

\value{
   An \code{mcmc} object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
\code{MCMCprobit} simulates from the posterior density of a probit
  regression model using data augmentation. The simulation
  proper is done in compiled C++ code to maximize efficiency.  Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior density sample.

  The model takes the following form:
  \deqn{y_i \sim \mathcal{B}ernoulli(\pi_i)}{y_i ~ Bernoulli(pi_i)}
  Where the inverse link function:
  \deqn{\pi_i = \Phi(x_i'\beta)}{pi_i = Phi(x_i'beta)}
  We assume a multivariate Normal prior on \eqn{\beta}{beta}:
    \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}{beta ~ N(b0,B0^(-1))}
  See Albert and Chib (1993) for estimation details.
  
 }
  
\references{
  Albert, J. H. and S. Chib. 1993. ``Bayesian Analysis of Binary and
  Polychotomous Response Data.'' \emph{J. Amer. Statist. Assoc.} 88, 669-679

  Albert, J. H. and S. Chib. 1995. ``Bayesian Residual Analysis for
  Binary Response Regression Models.'' \emph{Biometrika.} 82, 747-759.
      
  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2002.  
   \emph{Scythe Statistical Library 0.3.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}


\examples{\dontrun{data(birthwt)
posterior <- MCMCprobit(low~age+as.factor(race)+smoke, data=birthwt)}}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}, \code{\link[base]{glm}}}

