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

  If you are interested in fitting  K-dimensional item response theory
  models, or would rather identify the model by placing constraints
  on the item parameters, please see \code{\link[MCMCpack]{MCMCirtKd}}.
  }

\usage{
MCMCirt1d(datamatrix, theta.constraints=list(), burnin = 1000,
   mcmc = 20000, thin=1, verbose = FALSE, seed = NA, theta.start = NA, 
   alpha.start = NA, beta.start = NA, t0 = 0, T0 = 1, ab0=0, AB0=.25,
   store.item = FALSE, ... )  }

\arguments{
    \item{datamatrix}{The matrix of data.  Must be 0, 1, or missing values.  
     The rows of \code{datamatrix} correspond to subjects and the
     columns correspond to items.}

    \item{theta.constraints}{A list specifying possible simple equality
    or inequality constraints on the ability parameters. A typical
    entry in the list has one of three forms: \code{varname=c} which
    will constrain the ability parameter for the subject named
    \code{varname} to be equal to c, \code{varname="+"} which will
    constrain the ability parameter for the subject named \code{varname}
    to be positive, and  \code{varname="-"} which will constrain the
    ability parameter for the subject named \code{varname} to be
    negative. If x is a matrix without row names defaults names of
    ``V1",``V2", ... , etc will be used. See Rivers (2003) for a
    thorough discussion of identification of IRT models.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Gibbs iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
      the sampler is printed to the screen.  If TRUE, the iteration number
      is printed to the screen every 100 iterations.}
    
    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details.}

    \item{theta.start}{The starting values for the subject
    abilities (ideal points). This can either be a scalar or a
    column vector with dimension equal to the number of voters.  
    If this takes a scalar value, then that value will serve as
    the starting value for all of the thetas.  The default value
    of NA will choose the starting values based on an
    eigenvalue-eigenvector decomposition of the aggreement score matrix
    formed from the \code{datamatrix}.}

    \item{alpha.start}{The starting values for the
    \eqn{\alpha}{alpha} difficulty parameters. This can either be
    a scalar or a column vector with dimension equal to the
    number of items.   If this takes a scalar value, then that
    value will serve as the starting value for all of the alphas.
     The default value of NA will set the starting values based on
     a series of probit regressions that condition on the starting
     values of theta.}

    \item{beta.start}{The starting values for the
    \eqn{\beta}{beta} discrimination parameters. This can either
    be a scalar or a column vector with dimension equal to the
    number of items.   If this takes a scalar value, then that
    value will serve as the starting value for all of the betas. 
    The default value of NA will set the starting values based on a
    series of probit regressions that condition on the starting values
    of theta.}
    
    \item{t0}{A scalar parameter giving the prior mean of the subject
      abilities (ideal points).} 
    
    \item{T0}{A scalar parameter giving the prior precision (inverse
      variance) of the subject abilities (ideal points).}
    
    \item{ab0}{The prior mean of \code{(alpha, beta)}. Can be either a
      scalar or a 2-vector. If a scalar both means will be set to the
      passed value. The prior mean is assumed to be the same across all
      items.} 

    \item{AB0}{The prior precision of \code{(alpha, beta)}.This can
      either be ascalar or a 2 by 2 matrix. If this takes a scalar
      value, then that value times an identity matrix serves as the
      prior precision. The prior precision is assumed to be the same
      across all items.}
    
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
  
  The model takes the following form.  We assume that each
  subject has an subject ability (ideal point) denoted
  \eqn{\theta_j}{theta_j} and that each item has a difficulty
  parameter \eqn{\alpha_i}{alpha_i} and discrimination parameter
  \eqn{\beta_i}{beta_i}.  The observed choice by subject
  \eqn{j}{j} on item \eqn{i}{i} is the observed data matrix which
  is \eqn{(I \times J)}{(I * J)}.  We assume that the choice is
  dictated by an unobserved utility: \deqn{z_{i,j} = -\alpha_i +
  \beta_i \theta_j + \varepsilon_{i,j}}{z_ij = -alpha_i +
  beta_i*theta_j + epsilon_ij} Where the errors are assumed to be
  distributed standard Normal.  The parameters of interest are
  the subject abilities (ideal points) and the item parameters.

  We assume the following priors.  For the subject abilities (ideal points):
  \deqn{\theta_j \sim \mathcal{N}(t_{0},T_{0}^{-1})}{theta_j ~ N(t0, T0^{-1})}
  For the item parameters, the prior is:
  \deqn{\left[\alpha_i, \beta_i \right]' \sim \mathcal{N}_2 
  (ab_{0},AB_{0}^{-1})}{[alpha_i beta_i]' ~ N_2 (ab0, AB0^{-1})}
   
  The model is identified by the proper priors on the item parameters
  and constraints placed on the ability parameters.
}
  
\references{
   James H. Albert. 1992. ``Bayesian Estimation of Normal Ogive Item Response 
   Curves Using Gibbs Sampling." \emph{Journal of Educational Statistics}.  
   17: 251-269.
   
   Joshua Clinton, Simon Jackman, and Douglas Rivers. 2004. ``The Statistical 
   Analysis of Roll Call Data."  \emph{American Political Science Review}.
   98: 355-370.
   
   Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling." 
   Springer: New York.

   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  
   \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
   
   Douglas Rivers.  2004.  ``Identification of Multidimensional Item-Response
   Models."  Stanford University, typescript.
}


\examples{
   \dontrun{
   ## US Supreme Court Example with inequality constraints
   data(SupremeCourt)
   posterior1 <- MCMCirt1d(t(SupremeCourt),
                   theta.constraints=list(Scalia="+", Ginsburg="-"),
                   B0.alpha=.2, B0.beta=.2,
                   burnin=500, mcmc=100000, thin=20, verbose=TRUE,
                   store.item=TRUE)
   geweke.diag(posterior1)
   plot(posterior1)
   summary(posterior1)

   ## US Senate Example with equality constraints
   data(Senate)
   Sen.rollcalls <- Senate[,6:677]
   rownames(Sen.rollcalls) <- as.character(Senate$member)
   posterior2 <- MCMCirt1d(Sen.rollcalls,
                    theta.constraints=list(KENNEDY=-2, HELMS=2),
                    burnin=2000, mcmc=100000, thin=20, verbose=TRUE)
   geweke.diag(posterior2)
   plot(posterior2)
   summary(posterior2)
   }
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
\code{\link[MCMCpack]{MCMCirtKd}}}

