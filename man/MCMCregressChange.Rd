\name{MCMCregressChange}
\alias{MCMCregressChange}

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

\usage{MCMCregressChange(formula, data=parent.frame(),  m=1,
        b0=0, B0=0, c0=0.001, d0=0.001, sigma.mu=NA, sigma.var=NA, 
        a=NULL, b=NULL, mcmc=1000, burnin=1000, thin=1, verbose=0,
        seed=NA, beta.start=NA, P.start=NA,
        marginal.likelihood=c("none", "Chib95"), ...)}

\arguments{
   \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{m}{The number of changepoints.}

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
    By default, draws from the \code{Beta(0.9, 0.1)} 
    are used to construct a proper transition matrix for each raw except the last raw.}	



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
  \code{MCMCregressChange} simulates from the posterior distribution of
  a binary model with multiple changepoints.

  The model takes the following form:
  \deqn{y_t=x_t ' \beta_i + I(s_t=i)\varepsilon_{t},\;\; i=1, \ldots, k}{y_t=x_t'beta_i + I(s_t=i)epsilon_t, i=1,...,k.}
  Where \eqn{k}{k} is the number of states and \eqn{I(s_t=i)}{I(s_t=i)} is an indicator function
  that becomes 1 when a state at \eqn{t}{t} is \eqn{i}{i} and otherwise 0.

  The errors are assumed to be Gaussian in each regime:
  \deqn{I(s_t=i)\varepsilon_{t} \sim \mathcal{N}(0, \sigma^2_i)}{I(s_t=i)epsilon_t ~ N(0,
    sigma^2_i)}

  We assume standard, semi-conjugate priors:
  \deqn{\beta_i \sim \mathcal{N}(b_0,B_0^{-1}),\;\; i=1, \ldots, k}{beta_i ~ N(b0,B0^(-1)), i=1,...,k.}
  And:
  \deqn{\sigma^{-2}_i \sim \mathcal{G}amma(c_0/2, d_0/2),\;\; i=1, \ldots, k}{sigma^(-2)_i ~
    Gamma(c0/2, d0/2), i=1,...,k.}
  Where \eqn{\beta_i}{beta_i} and \eqn{\sigma^{-2}_i}{sigma^(-2)_i} are assumed
  \emph{a priori} independent.

  The simulation proper is done in compiled C++ code to maximize efficiency.   }

\references{
 Jong Hee Park. 2012. "A Change-point Approach to Intervention Analysis Using Bayesian Inference"
 Presented at the 2012 Annual Meeting of Korean Statistical Society.    

 Siddhartha Chib. 1995. "Marginal Likelihood from the Gibbs Output."
   \emph{Journal of the American Statistical Association}. 90:
   1313-1321.

 Siddhartha Chib. 1998. "Estimation and comparison of multiple change-point models."
   \emph{Journal of Econometrics}. 86: 221-241.


  Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.
  ``Output Analysis and Diagnostics for MCMC (CODA)'',
   \emph{R News}. 6(1): 7-11.
  \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.

}

\examples{
\dontrun{
set.seed(1119)
n <- 100
x1 <- runif(n)
true.beta1 <- c(2, -2)
true.beta2 <- c(0,  2)
true.Sigma <- c(1, 2)
true.s <- rep(1:2, each=n/2)

mu1 <- cbind(1, x1[true.s==1])%*%true.beta1
mu2 <- cbind(1, x1[true.s==2])%*%true.beta2

y <- as.ts(c(rnorm(n/2, mu1, sd=sqrt(true.Sigma[1])), rnorm(n/2, mu2, sd=sqrt(true.Sigma[2]))))
formula=y ~ x1

ols1 <- lm(y[true.s==1] ~x1[true.s==1])
ols2 <- lm(y[true.s==2] ~x1[true.s==2])
    
## prior
b0 <- 0
B0 <- 1
sigma.mu=sd(y)
sigma.var=var(y)
		
## models
model0 <-  MCMCregressChange(formula, m=0, b0=b0, B0=B0, mcmc=100, burnin=100, 
           sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
model1 <-  MCMCregressChange(formula, m=1, b0=b0, B0=B0, mcmc=100, burnin=100, 
           sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
model2 <-  MCMCregressChange(formula, m=2, b0=b0, B0=B0, mcmc=100, burnin=100,
           sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
model3 <-  MCMCregressChange(formula, m=3, b0=b0, B0=B0, mcmc=100, burnin=100,
           sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
model4 <-  MCMCregressChange(formula, m=4, b0=b0, B0=B0, mcmc=100, burnin=100,
           sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
model5 <-  MCMCregressChange(formula, m=5, b0=b0, B0=B0, mcmc=100, burnin=100,
           sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")

print(BayesFactor(model0, model1, model2, model3, model4, model5))
plotState(model1)
plotChangepoint(model1)

}
}

\keyword{models}

\seealso{\code{\link{plotState}}, \code{\link{plotChangepoint}}}
