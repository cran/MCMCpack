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

\usage{MCMCpoissonChangepoint(data,  m = 1, burnin = 1000, mcmc = 1000,
        thin = 1, verbose = 0, seed = NA, c0, d0, a = NULL, b = NULL,
        marginal.likelihood = c("none", "Chib95"), ...) }

\arguments{
    \item{data}{The data.}

    \item{m}{The number of changepoints.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burn-in.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0, the iteration number and the posterior density samples are printed to the screen every \code{verbose}th iteration.}

    \item{seed}{The seed for the random number generator.  If NA, current R 
    system seed is used.}

    \item{c0}{\eqn{c_0}{c0} is the shape parameter for Gamma prior on \eqn{\lambda}{lambda} 
    (the mean).}

    \item{d0}{\eqn{d_0}{d0} is the scale parameter for Gamma prior on \eqn{\lambda}{lambda} 
    (the mean).}

    \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default, 
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default, 
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated, and
    \code{Chib95} in which case the method of Chib (1995) is used.}

    \item{...}{further arguments to be passed}
}

\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.  The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} for each period, and the 
   log-likelihood of the model (\code{log.like}).
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
  
  Note that no default value is provided for Gamma priors.
  
  The simulation in this model-fitting is performed in R.
  }

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\references{
 Siddhartha Chib. 1995. "Marginal Likelihood from the Gibbs Output."
   \emph{Journal of the American Statistical Association}. 90:
   1313-1321.

 Siddhartha Chib. 1998. "Estimation and comparison of multiple change-point models."
   \emph{Journal of Econometrics}. 86: 221-241.


 Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}

\examples{
    \dontrun{
    ## generate event count time series data with two breakpoints at 100 and 200
    set.seed(1973)
    n           <-  300
    true.lambda <-  c(3, 2, 4)
    y1          <-  rpois(100, true.lambda[1]) 
    y2          <-  rpois(100, true.lambda[2]) 
    y3          <-  rpois(100, true.lambda[3]) 
    y           <-  c(y1, y2, y3)
   
    ## run the example   
    model1 <-  MCMCpoissonChangepoint(y, m=1, burnin = 100, mcmc = 100, 
                thin=1,  verbose=5, c0=3, d0=1,  marginal.likelihood =
                c("Chib95"))  
    model2 <-  MCMCpoissonChangepoint(y, m=2, burnin = 100, mcmc = 100, 
                thin=1, verbose=5, c0=3, d0=1,  marginal.likelihood =
                c("Chib95"))
    model3 <-  MCMCpoissonChangepoint(y, m=3, burnin = 100, mcmc = 100,
                thin=1, verbose=5, c0=3, d0=1,  marginal.likelihood =
                c("Chib95"))
    model4 <-  MCMCpoissonChangepoint(y, m=4, burnin = 100, mcmc = 100,
                thin=1, verbose=5, c0=3, d0=1,  marginal.likelihood = 
                c("Chib95"))
                
    ## perform model comparison  
    BayesFactor(model1, model2, model3, model4)
       
    ## Draw plots using the "right" model
    plotPostState(model2, y, m=2, legend.control=NULL) 
    plotPostChangepoint(model2, y, m=2)    
    }
}

\keyword{models}

\seealso{\code{\link{plotPostState}}, \code{\link{plotPostChangepoint}}}

