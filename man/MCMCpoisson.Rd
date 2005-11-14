\name{MCMCpoisson}
\alias{MCMCpoisson}
\title{Markov Chain Monte Carlo for Poisson Regression}
\description{
  This function generates a sample from the posterior distribution
  of a Poisson regression model using a random walk Metropolis
  algorithm. The user supplies data and priors,
  and a sample from the posterior distribution is returned as an mcmc
  object, which can be subsequently analyzed with functions 
  provided in the coda package.
  }
  
\usage{
MCMCpoisson(formula, data = parent.frame(), burnin = 1000, mcmc = 10000,
   thin = 1, tune = 1.1, verbose = 0, seed = NA,  beta.start = NA,
   b0 = 0, B0 = 0, ...) }

\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Metropolis iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    mcmc iterations must be divisible by this value.}

    \item{tune}{Metropolis tuning parameter. Can be either a positive
      scalar or a \eqn{k}{k}-vector, where \eqn{k}{k} is the length of
      \eqn{\beta}{beta}.Make sure that the
      acceptance rate is satisfactory (typically between 0.20 and 0.5)
      before using the posterior sample for inference.}

    
    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0 the iteration number,
    the current beta vector, and the Metropolis acceptance rate  are
    printed to the screen every \code{verbose}th iteration.} 

    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details.}
    
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
    as the prior precision of \eqn{\beta}{beta}. Default value of 0 is
    equivalent to an improper uniform prior for beta.}
     
    \item{\ldots}{further arguments to be passed}       

}

\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{\code{MCMCpoisson} simulates from the posterior distribution of
  a Poisson regression model using a random walk Metropolis
  algorithm. The simulation proper is done in compiled C++ code to
  maximize efficiency.  Please consult the coda documentation for a
  comprehensive list of functions that can be used to analyze the
  posterior sample. 
  
  The model takes the following form:
  \deqn{y_i \sim \mathcal{P}oisson(\mu_i)}{y_i ~ Poisson(mu_i)}
  Where the inverse link function:
  \deqn{\mu_i = \exp(x_i'\beta)}{mu_i = exp(x_i'beta)}
  We assume a multivariate Normal prior on \eqn{\beta}{beta}:
    \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}{beta ~ N(b0,B0^(-1))}

  The Metropois proposal distribution is centered at the current value of
  \eqn{\theta}{theta} and has variance-covariance \eqn{V = T
    (B_0 + C^{-1})^{-1} T }{V = T (B0 + C^{-1})^{-1} T}, where
  \eqn{T}{T} is a the diagonal positive definite matrix formed from the
  \code{tune}, \eqn{B_0}{B0} is the prior precision, and \eqn{C}{C} is
  the large sample variance-covariance matrix of the MLEs. This last
  calculation is done via an initial call to \code{glm}.
 }
  
\references{
      
  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  
   \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}


\examples{
   \dontrun{
   counts <- c(18,17,15,20,10,20,25,13,12)
   outcome <- gl(3,1,9)
   treatment <- gl(3,3)
   posterior <- MCMCpoisson(counts ~ outcome + treatment)
   plot(posterior)
   summary(posterior)
   }
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}, \code{\link[base]{glm}}}

