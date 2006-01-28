\name{MCmultinomdirichlet}
\alias{MCmultinomdirichlet}
\title{Monte Carlo Simulation from a Multinomial Likelihood with a 
   Dirichlet Prior}

\description{
  This function generates a sample from the posterior distribution
  of a multinomial likelihood with a Dirichlet prior.
  }
  
\usage{
MCmultinomdirichlet(y, alpha0, mc=1000, ...)
}

\arguments{
    \item{y}{A vector of data (number of successes for each category).}

    \item{alpha0}{The vector of parameters of the Dirichlet prior.}

    \item{mc}{The number of Monte Carlo draws to make.}
        
    \item{...}{further arguments to be passed}
    
}

\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCmultinomdirichlet} directly simulates from the posterior distribution. 
  This model is designed primarily for instructional use.  \eqn{\pi}{pi}
  is the parameter of interest of the multinomial distribution.  It is of
  dimension \eqn{(d \times 1)}{(d x 1)}. We assume
  a conjugate Dirichlet prior:
  \deqn{\pi \sim \mathcal{D}irichlet(\alpha_0)}{pi ~ Dirichlet(alpha0)}
  \eqn{y} is a \eqn{(d \times 1)}{(d x 1)} vector of observed data.
  }
  
\examples{
\dontrun{
## Example from Gelman, et. al. (1995, p. 78)
posterior <- MCmultinomdirichlet(c(727,583,137), c(1,1,1), mc=10000)
bush.dukakis.diff <- posterior[,1] - posterior[,2]
cat("Pr(Bush > Dukakis): ",
   sum(bush.dukakis.diff > 0) / length(bush.dukakis.diff), "\n")
hist(bush.dukakis.diff)
}
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}}

