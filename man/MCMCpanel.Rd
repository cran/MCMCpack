\name{MCMCpanel}
\alias{MCMCpanel}
\title{Markov Chain Monte Carlo for the General Linear Panel Model}
\description{
  MCMCpanel generates a posterior density sample from a General
  Linear Panel Model using Algorithm 2 of Chib and Carlin (1999).
  This model uses a multivariate Normal prior for the fixed
  effects parameters, a Wishart prior on the random effects
  precision matrix, and a Gamma prior on the conditional error
  precision. The user supplies data and priors, and a sample from
  the posterior density is returned as an mcmc object,
  which can be subsequently analyzed with functions provided in
  the coda package.
  }
  
\usage{
MCMCpanel(obs, Y, X, W, burnin = 1000, mcmc = 10000, thin = 5, 
    verbose = FALSE, seed = NA, sigma2.start = NA,
    D.start = NA, b0 = 0, B0 = 1, eta0, R0, nu0 = 0.001,
    delta0 = 0.001, ...)
   }

\arguments{
    \item{obs}{An \eqn{(nk \times 1)}{(nk * 1)} vector that contains unique
    observation numbers for each subject.}

    \item{Y}{An \eqn{(nk \times 1)}{(nk * 1)} vector of response variables, 
    stacked across all subjects.}

    \item{X}{An \eqn{(nk \times p)}{(nk * p)} matrix of fixed effects
    covariates, stacked across all subjects.}

    \item{W}{An \eqn{(nk \times q)}{(nk * q)} matrix of random effects
    covariates, stacked across all subjects.}

    \item{burnin}{The number of burnin iterations for the sampler.}
    
    \item{mcmc}{The number of Gibbs iterations for the sampler.}
    
    \item{thin}{The thinning interval used in the simulation.  The number of
    mcmc iterations must be divisible by this value.}

    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details.}


    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number
    and parameters are printed to the screen.}

    \item{sigma2.start}{The starting value for the conditional error
    variance. Default value of NA uses the least squares estimates.}

    \item{D.start}{The starting value for precision matrix of the random 
    effects.  This can either be a scalar or square matrix with dimension
    equal to the number of random effects. If this takes a scalar value, then
    that value multiplied by an identity matrix will be the starting
    value. Default value of NA uses an identity matrix multiplied by 0.5
    the OLS \eqn{\sigma^2}{sigma2} estimate. }
    
    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a
    column vector with dimension equal to the number of betas. If this takes
    a scalar value, then that value will serve as the prior mean for all
    of the betas.} 
  
    \item{B0}{The prior precision of \eqn{\beta}{beta}.
    This can either be a scalar
    or a square matrix with dimensions equal to the number of betas.  If
    this takes  a scalar value, then that value times an identity matrix 
    serves as the prior precision of beta. Default value of 0 is equivalent
    to an improper uniform prior for beta.} 
    
    \item{eta0}{The shape parameter for the Wishart prior on
    precision matrix for the random effects.}
    
    \item{R0}{The scale matrix for the Wishart prior on precision matrix for
    the random effects.}
    
    \item{nu0}{The shape parameter for the Gamma prior on the
    conditional error precision.} 

    \item{delta0}{The scale  parameter for the Gamma prior on
    the conditional error precision.}

    \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCMCpanel} simulates from the posterior density sample using 
  the blocked Gibbs sampler of Chib and Carlin (1999), Algorithm 2. 
  The simulation proper
  is done in compiled C++ code to maximize efficiency.  Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior density sample.
  
  The model takes the following form:
  \deqn{y_i = X_i \beta + W_i b_i + \varepsilon_i}{y_i = X_i * beta + W_i *
   b_i + epsilon_i}
  Where the random effects:
  \deqn{b_i \sim \mathcal{N}_q(0,D)}{b_i ~ N_q(0,D)}
  And the errors:
  \deqn{\varepsilon_i \sim \mathcal{N}(0, \sigma^2 I_k)}{epsilon_i ~ N(0,
   sigma^2 I_k)}
  We assume standard, conjugate priors:
  \deqn{\beta \sim \mathcal{N}(b0,B0^{-1})}{beta ~ N(b0,B0^(-1))}
  And:
  \deqn{\sigma^{-2} \sim \mathcal{G}amma(\nu_0/2, \delta_0/2)}{sigma^(-2) ~ 
  Gamma(nu0/2, delta0/2)}
  And:
  \deqn{D^{-1} \sim \mathcal{W}ishart(\eta_0, R_0^{-1})}{D^-1 ~ Wishart(eta0, 
   R0^-1)}
   See Chib and Carlin (1999) or Martin and Saunders (2002) for more details.
   
   \emph{NOTE: Unlike most models in MCMCpack, we do not provide default 
   parameters for the priors on the precision matrix for the random effects.}
   When fitting one of these models, it is of utmost importance to choose a 
   prior that reflects your prior beliefs about the random effects.  Using the
   \code{dwish} and \code{rwish} functions might be useful in choosing these
   values.  Also, the user is not allowed to specify a starting value
   for the \eqn{\beta}{beta} parameters, as they are simulated in the
   first block of the sampler.
   }
 
\references{
   Siddhartha Chib and Bradley P. Carlin. 1999. ``On MCMC Sampling in 
   Hierarchical Longitudinal Models." \emph{Statistics and Computing.} 9: 
   17-26.

   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  \emph{Scythe 
   Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Andrew D. Martin and Kyle L. Saunders. 2002. ``Bayesian Inference for 
   Political Science Panel Data.'' Paper presented at the 2002 Annual Meeting 
   of the American Political Science Association.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}

