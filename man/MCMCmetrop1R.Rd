% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCmetrop1R.R
\name{MCMCmetrop1R}
\alias{MCMCmetrop1R}
\title{Metropolis Sampling from User-Written R function}
\usage{
MCMCmetrop1R(
  fun,
  theta.init,
  burnin = 500,
  mcmc = 20000,
  thin = 1,
  tune = 1,
  verbose = 0,
  seed = NA,
  logfun = TRUE,
  force.samp = FALSE,
  V = NULL,
  optim.method = "BFGS",
  optim.lower = -Inf,
  optim.upper = Inf,
  optim.control = list(fnscale = -1, trace = 0, REPORT = 10, maxit = 500),
  ...
)
}
\arguments{
\item{fun}{The unnormalized (log)density of the distribution from which to
take a sample. This must be a function defined in R whose first argument is
a continuous (possibly vector) variable. This first argument is the point in
the state space at which the (log)density is to be evaluated. Additional
arguments can be passed to \code{fun()} by inserting them in the call to
\code{MCMCmetrop1R()}. See the Details section and the examples below for
more information.}

\item{theta.init}{Starting values for the sampling. Must be of the
appropriate dimension. It must also be the case that \code{fun(theta.init,
...)} is greater than \code{-Inf} if \code{fun()} is a logdensity or greater
than 0 if \code{fun()} is a density.}

\item{burnin}{The number of burn-in iterations for the sampler.}

\item{mcmc}{The number of MCMC iterations after burnin.}

\item{thin}{The thinning interval used in the simulation.  The number of
MCMC iterations must be divisible by this value.}

\item{tune}{The tuning parameter for the Metropolis sampling.  Can be either
a positive scalar or a \eqn{k}-vector, where \eqn{k} is the length of
\eqn{\theta}.}

\item{verbose}{A switch which determines whether or not the progress of the
sampler is printed to the screen.  If \code{verbose} is greater than 0 the
iteration number, the \eqn{\theta} vector, the function value, and
the Metropolis acceptance rate are sent to the screen every \code{verbose}th
iteration.}

\item{seed}{The seed for the random number generator.  If NA, the Mersenne
Twister generator is used with default seed 12345; if an integer is passed
it is used to seed the Mersenne twister.  The user can also pass a list of
length two to use the L'Ecuyer random number generator, which is suitable
for parallel computation.  The first element of the list is the L'Ecuyer
seed, which is a vector of length six or NA (if NA a default seed of
\code{rep(12345,6)} is used).  The second element of list is a positive
substream number. See the MCMCpack specification for more details.}

\item{logfun}{Logical indicating whether \code{fun} returns the natural log
of a density function (TRUE) or a density (FALSE).}

\item{force.samp}{Logical indicating whether the sampling should proceed if
the Hessian calculated from the initial call to \code{optim} routine to
maximize the (log)density is not negative definite. If
\code{force.samp==TRUE} and the Hessian from \code{optim} is non-negative
definite, the Hessian is rescaled by subtracting small values from it's main
diagonal until it is negative definite. Sampling proceeds using this
rescaled Hessian in place of the original Hessian from \code{optim}. By
default, if \code{force.samp==FALSE} and the Hessian from \code{optim} is
non-negative definite, an error message is printed and the call to
\code{MCMCmetrop1R} is terminated.

\emph{Please note that a non-negative Hessian at the mode is often an
indication that the function of interest is not a proper density. Thus,
\code{force.samp} should only be set equal to \code{TRUE} with great
caution.}}

\item{V}{The variance-covariance matrix for the Gaussian proposal
distribution. Must be a square matrix or \code{NULL}. If a square matrix,
\code{V} must have dimension equal to the length of \code{theta.init}. If
\code{NULL}, \code{V} is calculated from \code{tune} and an initial call to
\code{optim}. See the Details section below for more information. Unless the
log-posterior is expensive to compute it will typically be best to use the
default \code{V = NULL}.}

\item{optim.method}{The value of the \code{method} parameter sent to
\code{optim} during an initial maximization of \code{fun}. See \code{optim}
for more details.}

\item{optim.lower}{The value of the \code{lower} parameter sent to
\code{optim} during an initial maximization of \code{fun}. See \code{optim}
for more details.}

\item{optim.upper}{The value of the \code{upper} parameter sent to
\code{optim} during an initial maximization of \code{fun}. See \code{optim}
for more details.}

\item{optim.control}{The value of the \code{control} parameter sent to
\code{optim} during an initial maximization of \code{fun}. See \code{optim}
for more details.}

\item{\dots}{Additional arguments.}
}
\value{
An mcmc object that contains the posterior sample.  This object can
be summarized by functions provided by the coda package.
}
\description{
This function allows a user to construct a sample from a user-defined
continuous distribution using a random walk Metropolis algorithm.
}
\details{
MCMCmetrop1R produces a sample from a user-defined distribution using a
random walk Metropolis algorithm with multivariate normal proposal
distribution. See Gelman et al. (2003) and Robert & Casella (2004) for
details of the random walk Metropolis algorithm.

The proposal distribution is centered at the current value of
\eqn{\theta} and has variance-covariance \eqn{V}. If \eqn{V} is
specified by the user to be \code{NULL} then \eqn{V} is calculated
as: \eqn{V = T (-1\cdot H)^{-1} T }, where \eqn{T} is a the
diagonal positive definite matrix formed from the \code{tune} and
\eqn{H} is the approximate Hessian of \code{fun} evaluated at its
mode. This last calculation is done via an initial call to
\code{optim}.
}
\examples{

  \dontrun{

    ## logistic regression with an improper uniform prior
    ## X and y are passed as args to MCMCmetrop1R

    logitfun <- function(beta, y, X){
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
                              verbose=500, logfun=TRUE)

    raftery.diag(post.samp)
    plot(post.samp)
    summary(post.samp)
    ## ##################################################


    ##  negative binomial regression with an improper unform prior
    ## X and y are passed as args to MCMCmetrop1R
    negbinfun <- function(theta, y, X){
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
                              seed=list(NA,1))
    raftery.diag(post.samp)
    plot(post.samp)
    summary(post.samp)
    ## ##################################################


    ## sample from a univariate normal distribution with
    ## mean 5 and standard deviation 0.1
    ##
    ## (MCMC obviously not necessary here and this should
    ##  really be done with the logdensity for better
    ##  numerical accuracy-- this is just an illustration of how
    ##  MCMCmetrop1R works with a density rather than logdensity)

    post.samp <- MCMCmetrop1R(dnorm, theta.init=5.3, mean=5, sd=0.1,
                          thin=1, mcmc=50000, burnin=500,
                          tune=2.0, verbose=5000, logfun=FALSE)

    summary(post.samp)

  }

}
\references{
Siddhartha Chib; Edward Greenberg. 1995. ``Understanding the
Metropolis-Hastings Algorithm."  \emph{The American Statistician}, 49,
327-335.

Andrew Gelman, John B. Carlin, Hal S. Stern, and Donald B. Rubin. 2003.
\emph{Bayesian Data Analysis}. 2nd Edition. Boca Raton: Chapman & Hall/CRC.

Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
42(9): 1-21.  \doi{10.18637/jss.v042.i09}.

Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
Statistical Library 1.0.} \url{http://scythe.wustl.edu.s3-website-us-east-1.amazonaws.com/}.

Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
\url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.

Christian P. Robert and George Casella. 2004. \emph{Monte Carlo Statistical
Methods}. 2nd Edition. New York: Springer.
}
\seealso{
\code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}},
\code{\link[stats]{optim}}, \code{\link[mcmc]{metrop}}
}
\keyword{models}
