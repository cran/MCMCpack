\name{MCMCintervention}
\alias{MCMCintervention}

\title{Markov Chain Monte Carlo for a linear Gaussian Multiple Changepoint Model}
\description{
  This function generates a sample from the posterior distribution
  of a linear Gaussian model with multiple changepoints. The function uses
  the Markov chain Monte Carlo method of Chib (1998).
  The user supplies data and priors, and
  a sample from the posterior distribution is returned as an mcmc
  object, which can be subsequently analyzed with functions
  provided in the coda package.
}

\usage{MCMCintervention(y, data=parent.frame(),  m = 1, intervention,
        prediction.type=c("trend","ar"), change.type=c("fixed", "random", "all"),		   
        b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, sigma.mu = NA, sigma.var = NA,
        a = NULL, b = NULL, mcmc = 1000, burnin = 1000, thin = 1, verbose = 0,
        seed = NA, beta.start = NA, P.start = NA,
        marginal.likelihood = c("none", "Chib95"), ...)}

\arguments{
    \item{y}{data.}

    \item{data}{Data frame.}

    \item{m}{The number of changepoints.}
    
    \item{intervention}{The timing of intervention measured by its place in the response vector. 
    It should be larger than 0 and smaller than the length of the response vector. No default value is provided.}
    
    \item{prediction.type}{The type of local process. "trend" denotes the linear trend model and "ar" 
    denotes AR(1) process. By default, MCMCintervention uses the linear trend model.}
    
    \item{change.type}{The tyep of parameteric breaks. "all" denotes that all parameters have breaks, 
    "fixed" denotes that only the intercept and the slope have breaks, and "random" denotes that
    only the variance has breaks. By default, MCMCintervetnion assumes that all parameters have breaks.}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a
    scalar or a
    column vector with dimension equal to the number of betas. If this
    takes a scalar  value, then that value will serve as the prior
    mean for all of the betas.}

    \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either be a
    scalar or a square matrix with dimensions equal to the number of betas.
    If this
    takes a scalar value, then that value times an identity matrix serves
    as the prior precision of beta. Default value of 0 is equivalent to
    an improper uniform prior for beta.}

    \item{c0}{\eqn{c_0/2}{c0/2} is the shape parameter for the inverse
    Gamma prior on \eqn{\sigma^2}{sigma^2} (the variance of the
    disturbances). The amount of information in the inverse Gamma prior
    is something like that from \eqn{c_0}{c0} pseudo-observations.}

    \item{d0}{\eqn{d_0/2}{d0/2} is the scale parameter for the
    inverse Gamma prior on \eqn{\sigma^2}{sigma^2} (the variance of the
    disturbances). In constructing the inverse Gamma prior,
    \eqn{d_0}{d0} acts like the sum of squared errors from the
    \eqn{c_0}{c0} pseudo-observations.}

  	\item{sigma.mu}{The mean of the
    inverse Gamma prior on \eqn{\sigma^2}{sigma^2}. 
    \eqn{sigma.mu}{sigma.mu} and \eqn{sigma.var}{sigma.var} allow users to choose the inverse Gamma prior
    by choosing its mean and variance. }
  
  	\item{sigma.var}{The variacne of the
    inverse Gamma prior on \eqn{\sigma^2}{sigma^2}.
    \eqn{sigma.mu}{sigma.mu} and \eqn{sigma.var}{sigma.var} allow users to choose the inverse Gamma prior
    by choosing its mean and variance. }

    \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burnin.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0 the iteration number, the
    \eqn{\beta}{beta} vector, and the error variance are printed to
    the screen every \code{verbose}th iteration.}

    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA
    a default seed of \code{rep(12345,6)} is used).  The second element of
    list is a positive substream number. See the MCMCpack
    specification for more details.}

    \item{beta.start}{The starting values for the \eqn{\beta}{beta} vector.
    This can either be a scalar or a
    column vector with dimension equal to the number of betas.
    The default value of of NA will use the MLE
    estimate of \eqn{\beta}{beta} as the starting value.  If this is a
    scalar,  that value will serve as the starting value
    mean for all of the betas.}
    
    \item{P.start}{The starting values for the transition matrix. 
    A user should provide a square matrix with dimension equal to the number of states. 
    By default, draws from the \code{Beta(0.9, 0.1)} are used to construct 
    a proper transition matrix for each raw except the last raw.}

    \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated, and
    \code{Chib95} in which case the method of Chib (1995) is used.}

    \item{...}{further arguments to be passed}
}

\value{
   An mcmc object that contains the posterior sample.  This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} for each period, the
   log-likelihood of the model (\code{loglike}), and
   the log-marginal likelihood of the model (\code{logmarglike}).
}

\details{
  \code{MCMCintervention} simulates from the posterior distribution of
  a binary model with multiple changepoints.

  The model takes the following form:
  \deqn{y_t = x_t ' \beta_i + I(s_t = i)\varepsilon_{t},\;\; i = 1, \ldots, k}{y_t = x_t'beta_i + I(s_t = i)epsilon_t, i = 1,...,k.}
  Where \eqn{k}{k} is the number of states and \eqn{I(s_t = i)}{I(s_t = i)} is an indicator function
  that becomes 1 when a state at \eqn{t}{t} is \eqn{i}{i} and otherwise 0.

  The errors are assumed to be Gaussian in each regime:
  \deqn{I(s_t = i)\varepsilon_{t} \sim \mathcal{N}(0, \sigma^2_i)}{I(s_t = i)epsilon_t ~ N(0,
    sigma^2_i)}

  We assume standard, semi-conjugate priors:
  \deqn{\beta_i \sim \mathcal{N}(b_0,B_0^{-1}),\;\; i = 1, \ldots, k}{beta_i ~ N(b0,B0^(-1)), i = 1,...,k.}
  And:
  \deqn{\sigma^{-2}_i \sim \mathcal{G}amma(c_0/2, d_0/2),\;\; i = 1, \ldots, k}{sigma^(-2)_i ~
    Gamma(c0/2, d0/2), i = 1,...,k.}
  Where \eqn{\beta_i}{beta_i} and \eqn{\sigma^{-2}_i}{sigma^(-2)_i} are assumed
  \emph{a priori} independent.

  The simulation proper is done in compiled C++ code to maximize efficiency.   }

\references{
 Jong Hee Park. 2012. "A Change-point Approach to Intervention Analysis Using Bayesian Inference"
 Presented at the 2012 Annual Meeting of Korean Statistical Society.    

 Siddhartha Chib. 1998. "Estimation and comparison of multiple change-point models."
   \emph{Journal of Econometrics}. 86: 221-241.
}

\examples{
    \dontrun{
    Nile.std <- (Nile - mean(Nile))/sd(Nile)
    set.seed(1973)
    b0 <- matrix(c(0, 0), 2, 1);
    B0 <- diag(rep(0.25, 2))
    c0 = 2; d0 = 1
    
    ## Model Comparison
    ar0 <- MCMCintervention(Nile.std, m=0, prediction.type="ar", change.type = "all",
                        b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                        intervention = 29, marginal.likelihood = "Chib95")
    ar1all <- MCMCintervention(Nile.std, m=1, prediction.type="ar", change.type = "all",
                         b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                         intervention = 29, marginal.likelihood = "Chib95")
    ar1fixed <- MCMCintervention(Nile.std, m=1, prediction.type="ar", change.type = "fixed",
                         b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                         intervention = 29, marginal.likelihood = "Chib95")
    ar1random <- MCMCintervention(Nile.std, m=1, prediction.type="ar", change.type = "random",
                         b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                         intervention = 29, marginal.likelihood = "Chib95")
    tr0 <- MCMCintervention(Nile.std, m=0, prediction.type="trend", change.type = "all",
                        b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                        intervention = 29, marginal.likelihood = "Chib95")
    tr1all <- MCMCintervention(Nile.std, m=1, prediction.type="trend", change.type = "all",
                         b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                         intervention = 29, marginal.likelihood = "Chib95")
    tr1fixed <- MCMCintervention(Nile.std, m=1, prediction.type="trend", change.type = "fixed",
                         b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                         intervention = 29, marginal.likelihood = "Chib95")
    tr1random <- MCMCintervention(Nile.std, m=1, prediction.type="trend", change.type = "random",
                         b0=b0, B0=B0, c0=c0, d0=d0, mcmc = 1000, burnin = 1000, verbose = 500,
                         intervention = 29, marginal.likelihood = "Chib95")
    
    BayesFactor(ar0, ar1all, ar1fixed, ar1random, tr0, tr1all, tr1fixed, tr1random)

    par(mfrow=c(1,3))
    plotState(ar1fixed, start=1871, main="Hidden Regime Change")
    plotIntervention(ar1fixed, start=1871, main="Forward Analysis", alpha= 0.5, ylab="Nile River flow", xlab="Year")
    plotIntervention(ar1fixed, forward=FALSE, start=1871, main="Backward Analysis", alpha= 0.5, ylab="Nile River flow", xlab="Year")
    }
}

\keyword{models}

\seealso{\code{\link{plotIntervention}}}
