\name{MCMCmetrop1R}
\alias{MCMCmetrop1R}
\title{Metropolis Sampling from User-Written R function}
\description{
  This function allows a user to construct a sample from a user-defined
  R function using a random walk Metropolis algorithm. It assumes the
  parameters to be sampled are in a single block. 
}
\usage{
MCMCmetrop1R(fun, theta.init, burnin = 500, mcmc = 20000, thin = 1,
             tune = 1, verbose = 0, seed=NA, logfun = TRUE,
             optim.trace = 0, optim.REPORT = 10, optim.maxit = 500, ...)
}
\arguments{
  \item{fun}{The (log)density from which to take a sample. This should
    be a function defined in R and it should take a single
    argument. Additional arguments can be passed implicitly by either
    putting them in the global environment or by passing them as
    additional arguments to \code{MCMCmetrop1R()}. See the Details section 
    and the examples below. }
  \item{theta.init}{Starting values for the sampling. Must be of the
    appropriate dimension.}

  \item{burnin}{The number of burn-in iterations for the sampler.}

  \item{mcmc}{The number of MCMC iterations after burnin.}

  \item{thin}{The thinning interval used in the simulation.  The number of
    MCMC iterations must be divisible by this value.}

  \item{tune}{Can be either a positive scalar or a
    \eqn{k}{k}-vector, where \eqn{k}{k} is the length of
    \eqn{\theta}{theta}.}

  \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0 the iteration number, the
    \eqn{\theta}{theta} vector, the function value, and the Metropolis
    acceptance rate are sent to the screen every \code{verbose}th iteration.}
  
  \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details.} 
    
  \item{logfun}{Logical indicating whether \code{fun} returns the natural log
    of a density function (TRUE) or a density (FALSE).}

  \item{optim.trace}{The value of the \code{trace} parameter sent to
    \code{optim} during an initial maximization of \code{fun}. }

  \item{optim.REPORT}{The value of the \code{REPORT} parameter sent to
    \code{optim} during an initial maximization of \code{fun}.}

  \item{optim.maxit}{The value of the \code{maxit} parameter sent to
    \code{optim} during an initial maximization of \code{fun}.}

  \item{\dots}{Additional arguments.}
}
\details{
MCMCmetrop1R produces a sample from a user-defined (log)density using a
random walk Metropolis algorithm with multivariate normal proposal
distribution. See Gelman et al. (2003) and Robert & Casella (2004) for
details of the random walk Metropolis algorithm.

The proposal distribution is centered at the current value of
\eqn{\theta}{theta} and has variance-covariance \eqn{V = T
  (-1\cdot H)^{-1} T }{V = T (-1*H)^{-1} T}, where \eqn{T}{T} is a the
diagonal positive definite matrix formed from the \code{tune} and
\eqn{H}{H} is the approximate Hessian of \code{fun} evaluated at it's
mode. This last calculation is done via an initial call to
\code{optim}.

Passing data into a (log)density function does not conform to standard R 
practice.  The (log)density function should have only one argument -- the 
parameter vector.  The data used by the (log)density need to be either in the 
global environment or defined in the call to \code{MCMCmetrop1R()}, which will 
create an environment that contains the data and then call the (log)density function in that environment.  The names of the data objects 
are not checked, so be careful that these match the names used with the 
(log)density function.  \emph{NOTE:} This interface will likely be changed in 
a future implementation of MCMCpack.
}
\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}
\examples{
  \dontrun{
    
    ## logistic regression with an improper uniform prior
    ## X and y are passed as args to MCMCmetrop1R

    logitfun <- function(beta){
      eta <- X \%*\% beta
      p <- 1.0/(1.0+exp(-eta))
      sum( y * log(p) + (1-y)*log(1-p) )
    }
   
    x1 <- rnorm(1000)
    x2 <- rnorm(1000)
    Xdata <- cbind(1,x1,x2)
    p <- exp(.5 - x1 + x2)/(1+exp(.5 - x1 + x2))
    yvector <- rbinom(1000, 1, p)
    
    post.samp <- MCMCmetrop1R(logitfun, theta.init=c(0,0,0),
                              X=Xdata, y=yvector,
                              thin=1, mcmc=40000, burnin=500,
                              tune=c(1.5, 1.5, 1.5),
                              verbose=500, logfun=TRUE, optim.maxit=100)

    raftery.diag(post.samp)
    plot(post.samp)
    summary(post.samp)
    ## ##################################################
			      
    
    ## logistic regression with an improper uniform prior
    ## X and y are now in the global environment

    logitfun <- function(beta){
      eta <- X \%*\% beta
      p <- 1.0/(1.0+exp(-eta))
      sum( y * log(p) + (1-y)*log(1-p) )
    }
    
    x1 <- rnorm(1000)
    x2 <- rnorm(1000)
    X <- cbind(1,x1,x2)
    p <- exp(.5 - x1 + x2)/(1+exp(.5 - x1 + x2))
    y <- rbinom(1000, 1, p)
    
    post.samp <- MCMCmetrop1R(logitfun, theta.init=c(0,0,0),
                              thin=1, mcmc=40000, burnin=500,
                              tune=c(1.5, 1.5, 1.5),
                              verbose=500, logfun=TRUE, optim.maxit=100)

    raftery.diag(post.samp)
    plot(post.samp)
    summary(post.samp)
    ## ##################################################


    ##  negative binomial regression with an improper unform prior
    ## X and y are passed as args to MCMCmetrop1R
    negbinfun <- function(theta){
      k <- length(theta)
      beta <- theta[1:(k-1)]
      alpha <- exp(theta[k])
      mu <- exp(X \%*\% beta)
      log.like <- sum(
                      lgamma(y+alpha) - lfactorial(y) - lgamma(alpha) +
                      alpha * log(alpha/(alpha+mu)) +
                      y * log(mu/(alpha+mu))
                     )
    }

    n <- 1000
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    XX <- cbind(1,x1,x2)
    mu <- exp(1.5+x1+2*x2)*rgamma(n,1)
    yy <- rpois(n, mu)

    post.samp <- MCMCmetrop1R(negbinfun, theta.init=c(0,0,0,0), y=yy, X=XX,
                              thin=1, mcmc=35000, burnin=1000,
                              tune=1.5, verbose=500, logfun=TRUE,
                              optim.maxit=500, seed=list(NA,1))
    raftery.diag(post.samp)
    plot(post.samp)
    summary(post.samp)
    ## ##################################################
 
  }
}
\references{
  Andrew Gelman, John B. Carlin, Hal S. Stern, and Donald
  B. Rubin. 2003. \emph{Bayesian Data Analysis}. 2nd Edition. Boca
  Raton: Chapman & Hall/CRC.
  
  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  
   \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.

   Christian P. Robert and George Casella. 2004. \emph{Monte Carlo
     Statistical Methods}. 2nd Edition. New York: Springer.   
}
\seealso{\code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}, \code{\link[stats]{optim}}}

\keyword{models}
