\name{check.parameters}
\alias{check.parameters}
\title{MCMCpack Helper Function: Check Simulation Parameters}
\description{
  Check the simulation parameters to make sure they are
  appropriate. This is called by the MCMC* model functions.
}

\usage{
  check.parameters(burnin, mcmc, thin, fcn, tune)
}

\arguments{
  \item{burnin}{The number of burn-in iterations for the sampler.}
  
  \item{mcmc}{The number of MCMC iterations for the sampler.}
  
  \item{thin}{The thinning interval used in the simulation.  The
    number of MCMC iterations must be divisible by this value.}
  
  \item{fcn}{String containing function name (with no parentheses).
    It is used to throw informative errors.}
  
  \item{tune}{(Optional) tuning parameter for models that have a
    single Metropolis-Hastings tuning parameter.}
}
    

\value{
  Returns zero if exits cleanly.
}

\examples{
  check.parameters(100,1000,10,"MCMCregress")
}

\keyword{programming}
