\name{MCpoissongamma}
\alias{MCpoissongamma}
\title{Monte Carlo Simulation from a Poisson Likelihood with a Gamma Prior}

\description{
  This function generates a sample from the posterior distribution
  of a Poisson likelihood with a Gamma prior.
  }
  
\usage{
MCpoissongamma(y, alpha, beta, mc=1000, ...)
}

\arguments{
    \item{y}{A vector of counts (must be non-negative).}

    \item{alpha}{Gamma prior distribution shape parameter.}

    \item{beta}{Gamma prior distribution scale parameter.}

    \item{mc}{The number of Monte Carlo draws to make.}
        
    \item{...}{further arguments to be passed}
    
}

\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCpoissongamma} directly simulates from the posterior distribution. 
  This model is designed primarily for instructional use.  \eqn{\lambda}{lambda}
  is the parameter of interest of the Poisson distribution.
  We assume
  a conjugate Gamma prior:
  \deqn{\lambda \sim \mathcal{G}amma(\alpha, \beta)}{lambda ~ Gamma(alpha, beta)}
  \eqn{y} is a vector of counts.
  }
  
\examples{
\dontrun{
data(quine)
posterior <- MCpoissongamma(quine$Days, 15, 1, 5000)
summary(posterior)
plot(posterior)
grid <- seq(14,18,0.01)
plot(grid, dgamma(grid, 15, 1), type="l", col="red", lwd=3, ylim=c(0,1.3),
  xlab="lambda", ylab="density")
lines(density(posterior), col="blue", lwd=3)
legend(17, 1.3, c("prior", "posterior"), lwd=3, col=c("red", "blue"))
}
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}}

