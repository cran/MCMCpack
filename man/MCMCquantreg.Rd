\name{MCMCquantreg}

\alias{MCMCquantreg}

\title{Bayesian quantile regression using Gibbs sampling}

\description{
  This function fits quantile regression models under Bayesian inference.  
  The function samples from the posterior distribution using Gibbs sampling with data augmentation. 
 A multivariate normal prior is assumed for \eqn{\beta}{beta} and an inverse Gamma prior 
 is assumed for \eqn{\sigma}{sigma}. The user supplies the prior parameters.
 A sample of the posterior distribution is returned as an mcmc object, 
 which can then be analysed by functions in the coda package.
}

\usage{
MCMCquantreg(formula, p=0.5, data = NULL, burnin = 1000,
   mcmc = 10000, thin = 1, verbose = 0, seed = NA,
   beta.start = NA, b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, ...)
}

\arguments{

  \item{formula}{ Model formula. }

  \item{p}{The quantile of interest. Must be between 0 and 1. The default value of 0.5 corresponds to median regression.}


  \item{data}{ Data frame. }

  \item{burnin}{ The number of burn-in iterations for the sampler. }

  \item{mcmc}{ The number of MCMC iterations after burnin. }

  \item{thin}{ The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value. }
  
  \item{verbose}{ A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0 the iteration number and the most recently sampled values of  
    \eqn{\beta}{beta} and \eqn{\sigma}{sigma} are printed to 
    the screen every \code{verbose}th iteration. }

  \item{seed}{ The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details. }

  \item{beta.start}{ The starting values for \eqn{\beta}{beta}.
    This can either be a scalar or a
    column vector with dimension equal to the dimension of \eqn{\beta}{beta}.
    The default value of NA will use the OLS
    estimate \eqn{\hat{\beta}}{beta^hat} with 
    \eqn{\hat{\sigma}\Phi^{-1}(p)}{sigma^hat*Phi^(-1)(p)} added on to the first element of \eqn{\hat{\beta}}{beta^hat} as the starting value. 
    (\eqn{\hat{\sigma}^2}{(sigma^hat)^2} denotes the usual unbiased estimator of
    \eqn{\sigma^2}{sigma^2} under ordinary mean regression and 
    \eqn{\Phi^{-1}(p)}{Phi^(-1)(p)} denotes the inverse of the
    cumulative density function of the standard normal distribution.)  
     Note that the default value assume that an intercept is included in the model.
    If a scalar is given,  that value will serve as the starting value
    for all \eqn{\beta}{beta}. }

  \item{b0}{ The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a
    column vector with dimension equal to the dimension of \eqn{\beta}{beta}. If this
    takes a scalar value, then that value will serve as the prior
    mean for all \eqn{\beta}{beta}. }

  \item{B0}{ The prior precision of \eqn{\beta}{beta}.  This can either be a 
    scalar or a square matrix with dimensions equal to the number of betas.  
    If this
    takes a scalar value, then that value times an identity matrix serves
    as the prior precision of \eqn{\beta}{beta}. Default value of 0 is equivalent to
    an improper uniform prior for \eqn{\beta}{beta}. }

  \item{c0}{ \eqn{c_0/2}{c0/2} is the shape parameter for the inverse
    Gamma prior on \eqn{\sigma}{sigma}. }

  \item{d0}{ \eqn{d_0/2}{d0/2} is the scale parameter for the
    inverse Gamma prior on \eqn{\sigma}{sigma}.}

    \item{\dots}{ further arguments to be passed }
}

\details{
  \code{MCMCquantreg} simulates from the posterior distribution using 
  a partially collapsed Gibbs sampler with data augmentation (see \url{http://people.brunel.ac.uk/~mastkky/}). 
  This involves a multivariate Normal draw for \eqn{\beta}{beta}, and an
  inverse Gamma draw for \eqn{\sigma}{sigma}. 
  The augmented data are drawn conditionally from the inverse Gaussian distribution. The simulation
  is carried out in compiled C++ code to maximise efficiency.  Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyse the posterior sample. 
  
  The model takes the following form:
  \deqn{y_i = x_i ' \beta + \varepsilon_i}{y_i = x_i'beta + epsilon_i}
  The errors are assumed to have an Asymmetric Laplace distribution with pth quantile equal to zero 
  and shape parameter equal to \eqn{\sigma}{sigma}:
  \deqn{\varepsilon_i(p) \sim \mathcal{AL}(0, \sigma, p)}{epsilon_i(p) ~ AL(0,
    sigma, p),}
  where \eqn{\beta}{beta} and \eqn{\sigma}{sigma} depend on p. 
  We assume standard, semi-conjugate priors:
  \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}{beta ~ N(b0,B0^(-1)),}
  and
  \deqn{\sigma^{-1} \sim \mathcal{G}amma(c_0/2, d_0/2)}{sigma^(-1) ~
    Gamma(c0/2, d0/2),}
  where \eqn{\beta}{beta} and \eqn{\sigma^{-1}}{sigma^(-1)} are assumed 
  \emph{a priori} independent of each other and the parameters associated with the other quantiles. Only starting values for
  \eqn{\beta}{beta} are allowed for this sampler.
}

\value{
  An mcmc object that contains the posterior sample.  This 
   object can be summarised by functions provided by the coda package.
}

\references{ Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  
   \emph{Scythe Statistical Library 1.2.} \url{http://scythe.wustl.edu}.
   
   Craig Reed and Keming Yu. 2009. ``A Partially Collapsed Gibbs Sampler for Bayesian Quantile Regression." Technical Report.
   
   Keming Yu and Jin Zhang. 2005. ``A Three Parameter Asymmetric Laplace Distribution and it's extensions."
\emph{Communications in Statistics - Theory and Methods}, 34, 1867-1879.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.}
     
\author{Craig Reed}

\seealso{ \code{\link[MCMCpack]{MCMCregress}}, \code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}, \code{\link[stats]{lm}}, \code{\link[quantreg]{rq}} }

\examples{
\dontrun{

x<-rep(1:10,5)
y<-rnorm(50,mean=x)
posterior_50 <- MCMCquantreg(y~x)
posterior_95 <- MCMCquantreg(y~x, p=0.95, verbose=10000,
    mcmc=50000, thin=10, seed=2)
plot(posterior_50)
plot(posterior_95)
raftery.diag(posterior_50)
autocorr.plot(posterior_95)
summary(posterior_50)
summary(posterior_95)
}
}

\keyword{models}
