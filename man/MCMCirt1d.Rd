\name{MCMCirt1d}
\alias{MCMCirt1d}
\title{Markov chain Monte Carlo for One Dimensional Item Response Theory
   Model}
\description{
  This function generates a posterior density sample from a one
  dimentional item response theory (IRT) model, with Normal
  priors on the subject abilities (ideal points), and
  multivariate Normal priors on the item parameters.  The user
  supplies data and priors, and a sample from the posterior
  density is returned as an mcmc object, which can be
  subsequently analyzed with functions provided in the coda
  package.

  \emph{NOTE:} This implementation of this model
  may be deprecated in future releases.  We have implemented
  a general K-dimensional item response theory model that
  allows the user to place arbitrary constraints on item and
  subject parameters.  
  }
  
\usage{
MCMCirt1d(datamatrix, theta.fixed = 1, burnin = 500, mcmc = 1000,
   thin=5, verbose = FALSE, seed = 0, theta.start = NA, 
   alpha.start = NA, beta.start = NA, t0 = 0, T0 = 1, b0.alpha = 0,
   b0.beta = 0, B0.alpha = 1, B0.beta = 1, B0.corr = 0,
   store.item = FALSE, ... )  }

\arguments{
    \item{datamatrix}{The matrix of data.  Must be 0, 1, or missing values.  
    It is of dimensionality items by subjects.}

    \item{theta.fixed}{Identifying restriction.  This is the
    subject whose subject ability (ideal point) is constrained to
    be negative.  It makes most sense to choose someone who is
    extreme on the latent scale.  Make sure to check the
    posterior density sample to ensure the sampler is sampling
    from only one posterior mode.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Gibbs iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number and
    the subject abilities (ideal points) are printed to the screen.}
    
    \item{seed}{The seed for the random number generator.  The code uses the
    Mersenne Twister, which requires an integer as an input.  If nothing
    is provided, the Scythe default seed is used.}

    \item{theta.start}{The starting values for the subject
    abilities (ideal points).   This can either be a scalar or a
    column vector with dimension equal to the number of voters.  
    If this takes a scalar value, then that value will serve as
    the starting value for all of the thetas.  The default value
    of NA will choose the starting values randomly.}

    \item{alpha.start}{The starting values for the
    \eqn{\alpha}{alpha} difficulty parameters. This can either be
    a scalar or a column vector with dimension equal to the
    number of items.   If this takes a scalar value, then that
    value will serve as the starting value for all of the alphas.
     The default value of NA will choose the starting values
    randomly.}

    \item{beta.start}{The starting values for the
    \eqn{\beta}{beta} discrimination parameters. This can either
    be a scalar or a column vector with dimension equal to the
    number of items.   If this takes a scalar value, then that
    value will serve as the starting value for all of the betas. 
    The default value of NA will choose the starting values
    randomly.}
    
    \item{t0}{The prior means of the subject abilities (ideal
    points), stacked for all subjects.  This can either be a
    scalar or a column vector with dimension equal to the number
    of thetas. If this takes a scalar value, then that value will
    serve as the prior mean for all of the betas.}
    
    \item{T0}{The prior variances of the subject abilities (ideal
    points), stacked for all subjects.   This can either be a
    scalar or a column vector with dimension equal to the number
    of thetas. If this takes a scalar  value, then that value
    will serve as the prior variance for all of the thetas.}
    
    \item{b0.alpha}{The prior means of the difficulty parameters,
    stacked for all items. This can either be a scalar or a
    column vector with dimension equal to the number of alphas.
    If this takes a scalar value, then that value will serve as
    the prior mean for all of the alphas.}

    \item{b0.beta}{The prior means of the discrimination
    parameters, stacked for all items. This can either be a
    scalar or a column vector with dimension equal to the number
    of betas. If this takes a scalar value, then that value will
    serve as the prior mean for all of the betas.}
    
    \item{B0.alpha}{The prior variances of the difficulty
    parameters, stacked for all items. This can either be a
    scalar or a column vector with dimension equal to the number
    of alphas. If this takes a scalar  value, then that value
    will serve as the prior variance for all of the alphas.}
    
    \item{B0.beta}{The prior variances of the discrimination
    parameters, stacked for all items. This can either be a
    scalar or a column vector with dimension equal to the number
    of betas. If this takes a scalar  value, then that value will
    serve as the prior variance for all of the betas.}
    
    \item{B0.corr}{The prior correlations of the difficulty and
    discrimination parameters, stacked for all items.  These are
    converted into covariances to complete the multivariate
    Normal prior.  This can either be a scalar or a column vector
    with dimension equal to the number of items. If this takes a
    scalar  value, then that value will serve as the prior
    correlation for all of the items.}
    
    \item{store.item}{A switch that determines whether or not to
    store the item parameters for posterior analysis. 
    \emph{NOTE: This takes an enormous amount of memory, so
    should only be used if the chain is thinned heavily, or for
    applications with a small number of items}.  By default, the
    item parameters are not stored.}
    
    \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCMCirt1d} simulates from the posterior density using
  standard Gibbs sampling using data augmentation (a Normal draw
  for the subject abilities, a multivariate Normal
  draw for the item parameters, and a truncated Normal draw for
  the latent utilities). The simulation proper is done in
  compiled C++ code to maximize efficiency.  Please consult the
  coda documentation for a comprehensive list of functions that
  can be used to analyze the posterior density sample.
  
  The default number of burnin and mcmc iterations is much
  smaller than the typical default values in MCMCpack.  This is
  because fitting this model is extremely computationally
  expensive.  It does not mean that this small of a number of
  scans will yield good estimates. If the verbose option is
  chosen, output will be printed to the screen every fifty
  iterations.  The priors of this model need to be proper for
  identification purposes.  The user is asked to provide prior
  means and variances \emph{(not precisions)} for the item
  parameters and the subject parameters.
  
  The model takes the following form.  We assume that each
  subject has an subject ability (ideal point) denoted
  \eqn{\theta_j}{theta_j} and that each item has a difficulty
  parameter \eqn{\alpha_i}{alpha_i} and discrimination parameter
  \eqn{\beta_i}{beta_i}.  The observed choice by subject
  \eqn{j}{j} on item \eqn{i}{i} is the observed data matrix which
  is \eqn{(I \times J)}{(I * J)}.  We assume that the choice is
  dictated by an unobserved utility: \deqn{z_{i,j} = \alpha_i +
  \beta_i \theta_j + \varepsilon_{i,j}}{z_ij = alpha_i +
  beta_i*theta_j + epsilon_ij} Where the errors are assumed to be
  distributed standard Normal.  The parameters of interest are
  the subject abilities (ideal points) and the item parameters.

  We assume the following priors.  For the subject abilities (ideal points):
  \deqn{\theta_j \sim \mathcal{N}(t_{0,j},T_{0,j})}{theta_j ~ N(t_0j, T0j)}
  Note that this implies a separate prior mean and variance for each
  subject.  For the item parameters, the prior is:
  \deqn{\left[\alpha_i, \beta_i \right]' \sim \mathcal{N}_2 
  (b_{0,i},B_{0,i})}{[alpha_i beta_i]' ~ N_2 (b_0,i, B_0,i)}
  Again, there is a separate prior for each item parameter.  When 
  supplying priors of the item parameters to the function, the user provides 
  each element of the mean vector and the covariance matrix (taking 
  correlations instead of covariances for convenience).
  
  The model is identified by the proper priors on the subject
  abilities (ideal points) and item parameters.  The
  \code{theta.fixed} solves the rotational invariance problem by
  constraining the sampler to one of the two identical posterior
  modes.  This function differs from the \code{MCMCirtKd} function
  in the manner in which the model is identified and in the assumed
  priors.
  }
  
\references{
   James H. Albert. 1992. ``Bayesian Estimation of Normal Ogive Item Response 
   Curves Using Gibbs Sampling." \emph{Journal of Educational Statistics}.  
   17: 251-269.
   
   Joshua Clinton, Simon Jackman, and Douglas Rivers. 2000. ``The Statistical 
   Analysis of Legislative Behavior: A Unified Approach." Paper presented at 
   the Annual Meeting of the Political Methodology Society.
   
   Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling." 
   Springer: New York.

   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2003.  
   \emph{Scythe Statistical Library 0.4.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}


\examples{
   \dontrun{
   data(SupremeCourt)
   posterior1 <- MCMCirt1d(SupremeCourt, verbose=TRUE,
      burnin=10000, mcmc=50000)
   plot(posterior1)
   summary(posterior1)
   
   data(Senate)
   posterior2 <- MCMCirt1d(t(Senate[,6:677]), theta.fixed = 9, verbose=TRUE,
      burnin=1000, mcmc=5000)
   plot(posterior2)
   summary(posterior2)
   }
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
\code{\link[MCMCpack]{MCMCirtKd}}}

