\name{MCMCpoissonChangepoint}
\alias{MCMCpoissonChangepoint}

\title{Markov Chain Monte Carlo for a Poisson Multiple Changepoint Model}
\description{
  This function generates a sample from the posterior distribution
  of a Poisson model with multiple changepoints. The function uses
  the Markov chain Monte Carlo method of Chib (1998).
  The user supplies data and priors, and
  a sample from the posterior distribution is returned as an mcmc
  object, which can be subsequently analyzed with functions
  provided in the coda package.
}

\usage{MCMCpoissonChangepoint(data,  m = 1, c0 = NA, d0 = NA, a = NULL, b = NULL,
            burnin = 10000, mcmc = 10000, thin = 1, verbose = 0, seed = NA,
            lambda.start = NA, P.start = NA,
            marginal.likelihood = c("none", "Chib95"), ...)  }

\arguments{
    \item{data}{The data.}

    \item{m}{The number of changepoints.}

     \item{c0}{\eqn{c_0}{c0} is the shape parameter for Gamma prior on \eqn{\lambda}{lambda}
    (the mean). No default value is provided.}

    \item{d0}{\eqn{d_0}{d0} is the scale parameter for Gamma prior on \eqn{\lambda}{lambda}
    (the mean). No default value is provided.}

    \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burn-in.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0, the iteration number and the posterior density samples are printed to the screen every \code{verbose}th iteration.}

    \item{seed}{The seed for the random number generator.  If NA, current R
    system seed is used.}
	
    \item{lambda.start}{The starting values for the lambda vector. This can either be a scalar or a column vector with dimension equal to the number of lambdas. The default value of NA will use draws from the Uniform distribution with the same boundary with the data as the starting value. If this is a scalar, that value will serve as the starting value mean for all of the lambdas.}	

    \item{P.start}{The starting values for the transition matrix. A user should provide a square matrix with dimension equal to the number of states. By default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper transition matrix for each raw except the last raw.}	

    \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated, and
    \code{Chib95} in which case the method of Chib (1995) is used.}

    \item{...}{further arguments to be passed}
}

\value{
   An mcmc object that contains the posterior sample. This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} for each period, and
   the log-marginal likelihood of the model (\code{logmarglike}).
}

\details{
  \code{MCMCpoissonChangepoint} simulates from the posterior distribution of
  a Poisson model with multiple changepoints.

  The model takes the following form:
  \deqn{Y_t \sim \mathcal{P}oisson(\lambda_i),\;\; i = 1, \ldots, k}{Y_t ~ Poisson(lambda_i), i = 1,...,k.}
  Where \eqn{k}{k} is the number of states.

  We assume Gamma priors for \eqn{\lambda_{i}}{lambda_i} and Beta priors for transition probabilities:
  \deqn{\lambda_i \sim \mathcal{G}amma(c_0, d_0)}{lambda_i ~ Gamma(c0, d0)}
  \deqn{p_{ii} \sim \mathcal{B}eta{a}{b},\;\; i = 1, \ldots, k}{p_ii ~ Beta(a, b), i = 1,...,k.}
  Where \eqn{k}{k} is the number of states.

  The simulation is done in compiled C++ code to maximize efficiency.
  }

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\references{
 Siddhartha Chib. 1995. "Marginal Likelihood from the Gibbs Output."
   \emph{Journal of the American Statistical Association}. 90:
   1313-1321.

 Siddhartha Chib. 1998. "Estimation and comparison of multiple change-point models."
   \emph{Journal of Econometrics}. 86: 221-241.
}

\examples{
    \dontrun{
    set.seed(1973)
    n           <-  100
    true.lambda <-  c(5, 8)
    y1          <-  rpois(n/2, true.lambda[1])
    y2          <-  rpois(n/2, true.lambda[2])
    y           <-  c(y1, y2)

    model1 <-  MCMCpoissonChangepoint(y, m=1, c0=6.85, d0=1, verbose = 10000, marginal.likelihood = "Chib95")
    model2 <-  MCMCpoissonChangepoint(y, m=2, c0=6.85, d0=1, verbose = 10000, marginal.likelihood = "Chib95")
    model3 <-  MCMCpoissonChangepoint(y, m=3, c0=6.85, d0=1, verbose = 10000, marginal.likelihood = "Chib95")
    model4 <-  MCMCpoissonChangepoint(y, m=4, c0=6.85, d0=1, verbose = 10000, marginal.likelihood = "Chib95")
    model5 <-  MCMCpoissonChangepoint(y, m=5, c0=6.85, d0=1, verbose = 10000, marginal.likelihood = "Chib95")

    print(BayesFactor(model1, model2, model3, model4, model5))

    ## Draw plots using the "right" model
    plotState(model1)
    plotChangepoint(model1)
    }
}

\keyword{models}

\seealso{\code{\link{plotState}}, \code{\link{plotChangepoint}}}
