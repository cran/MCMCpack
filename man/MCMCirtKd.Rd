\name{MCMCirtKd}
\alias{MCMCirtKd}
\title{Markov chain Monte Carlo for K-Dimensional Item Response Theory
   Model}
\description{
  This function generates a posterior density sample from a
  K-dimensional item response theory (IRT) model, with standard
  Normal priors on the subject abilities (ideal points), and
  Normal priors on the item parameters.  The user
  supplies data and priors, and a sample from the posterior
  density is returned as an mcmc object, which can be
  subsequently analyzed with functions provided in the coda
  package.
}
  
\usage{
MCMCirtKd(datamatrix, dimensions, item.constraints=list(),
   burnin = 1000, mcmc = 10000, thin=1, verbose = FALSE, seed = NA,
   alphabeta.start = NA, b0 = 0, B0=0, store.item = FALSE,
   store.ability=TRUE, drop.constantvars=TRUE, ... )  }

\arguments{
    \item{datamatrix}{The matrix of data.  Must be 0, 1, or missing values.  
    It is of dimensionality items by subjects.}

    \item{dimensions}{The number of dimensions in the latent space.}
    
    \item{item.constraints}{List of lists specifying possible equality
    or simple inequality constraints on the item parameters. A typical
    entry in the list has one of three forms: \code{rowname=list(d,c)}
    which will constrain the dth item parameter for the item named
    rowname to be equal to c, \code{rowname=list(d,"+")} which will
    constrain the dth item parameter for the item named rowname to be
    positive, and\code{rowname=list(d, "-")} which will constrain the dth
    item parameter for the item named varname to be negative. If x is a
    matrix without row names defaults names of ``V1", ``V2", ... , etc
    will be used. In a d dimensional model, the first item parameter for
    item \eqn{i}{i} is the difficulty parameter (\eqn{\alpha_i}{alpha_i}),
    the second item parameter is the discrimation parameter on dimension
    1 (\eqn{\beta_{i,1}{beta_{i,1}}}), the third item parameter is the
    discrimation parameter on dimension 2
    (\eqn{\beta_{i,2}{beta_{i,2}}}), ...,  and the (d+1)th item parameter
    is the discrimation parameter on dimension d
    (\eqn{\beta_{i,1}{beta_{i,1}}}). 
    The item difficulty parameters (\eqn{\alpha}{alpha}) should
    generally not be constrained. 
    }

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number and
    the subject abilities (ideal points) are printed to the screen.}
    
    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
      Twister generator is used with default seed 12345; if an integer is 
      passed it is used to seed the Mersenne twister.  The user can also
      pass a list of length two to use the L'Ecuyer random number generator,
      which is suitable for parallel computation.  The first element of the
      list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
      a default seed of \code{rep(12345,6)} is used).  The second element of 
      list is a positive substream number. See the MCMCpack 
      specification for more details.}

    \item{alphabeta.start}{The starting values for the
    \eqn{\alpha}{alpha} and \eqn{\beta}{beta} difficulty and
    discrimination parameters. If \code{alphabeta.start} is set to a
    scalar the starting value for all unconstrained item parameters will
    be set to that scalar. If \code{alphabeta.start} is a matrix of
    dimension \eqn{(K+1) \times items}{(K+1) x items} then the
    \code{alphabeta.start} matrix is used as the starting values (except
    for equality-constrained elements). If \code{alphabeta.start} is set
    to \code{NA} (the default) then starting values for unconstrained
    elements are set to values generated from a series of proportional
    odds logistic regression fits, and starting values for inequality
    constrained elements are set to either 1.0 or -1.0 depending on the
    nature of the constraints. }

    \item{b0}{The prior means of the
    \eqn{\alpha}{alpha} and \eqn{\beta}{beta} difficulty and
    discrimination parameters, stacked for all items.
    If a scalar is passed, it
    is used as the prior mean for all items.}
    
    \item{B0}{The prior precisions (inverse variances) of the
    independent Normal prior on the item parameters.
    Can be either a scalar or a matrix of dimension
    \eqn{(K+1) \times items}{(K+1) x items}.}

    \item{store.item}{A switch that determines whether or not to
    store the item parameters for posterior analysis. 
    \emph{NOTE: This takes an enormous amount of memory, so
    should only be used if the chain is thinned heavily, or for
    applications with a small number of items}.  By default, the
    item parameters are not stored.}
 
   \item{store.ability}{A switch that determines whether or not to store
    the subject abilities for posterior analysis.  By default, the
    item parameters are all stored.}

  \item{drop.constantvars}{A switch that determines whether or not
        items and subjects that have no variation
	should be deleted before fitting the model. Default = TRUE.}
  
    \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCMCirtKd} simulates from the posterior density using
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
  means and precisions \emph{(not variances)} for the item
  parameters and the subject parameters.
  
  The model takes the following form.  We assume that each
  subject has an subject ability (ideal point) denoted
  \eqn{\theta_j}{theta_j} \eqn{(K \times 1)}{(K x 1)},
  and that each item has a difficulty
  parameter \eqn{\alpha_i}{alpha_i} and discrimination parameter
  \eqn{\beta_i}{beta_i} \eqn{(K \times 1)}{(K x 1)}.
  The observed choice by subject
  \eqn{j}{j} on item \eqn{i}{i} is the observed data matrix which
  is \eqn{(I \times J)}{(I * J)}.  We assume that the choice is
  dictated by an unobserved utility: \deqn{z_{i,j} = \alpha_i +
  \beta_i' \theta_j + \varepsilon_{i,j}}{z_ij = alpha_i +
  beta_i'*theta_j + epsilon_ij} Where the errors are assumed to be
  distributed standard Normal.  The parameters of interest are
  the subject abilities (ideal points) and the item parameters.

  We assume the following priors.  For the subject abilities (ideal points)
  we assume independent standard Normal priors:
  \deqn{\theta_{j,k} \sim \mathcal{N}(0,1)}{theta_j,k ~ N(0, 1)}
  These cannot be changed by the user.
  For each item parameter, we assume independent Normal priors:
  \deqn{\left[\alpha_i, \beta_i \right]' \sim \mathcal{N}_{(K+1)} 
  (b_{0,i},B_{0,i})}{[alpha_i beta_i]' ~ N_(K+1) (b_0,i, B_0,i)}
  Where \eqn{B_{0,i}}{B_0,i} is a diagonal matrix.
  One can specify a separate prior mean and precision
  for each item parameter.
  
  The model is identified by the constraints on the item parameters
  (see Jackman 2001).  The user cannot place constraints on the subect
  abilities.  This identification scheme differs from that in
  \code{MCMCirt1d}, which uses a single directional constraint on
  one subject ability.  However, in our experience, using subject 
  ability constraints for models in greater than one dimension does not work 
  particularly well.
  }
  
\references{
   James H. Albert. 1992. ``Bayesian Estimation of Normal Ogive Item Response 
   Curves Using Gibbs Sampling." \emph{Journal of Educational Statistics}.  
   17: 251-269.
   
   Joshua Clinton, Simon Jackman, and Douglas Rivers. 2000. ``The Statistical 
   Analysis of Legislative Behavior: A Unified Approach." Paper presented at 
   the Annual Meeting of the Political Methodology Society.
   
   Simon Jackman. 2001. ``Multidimensional Analysis of Roll Call Data
   via Bayesian Simulation.'' \emph{Political Analysis.} 9: 227-241.
   
   Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling." 
   Springer: New York.

   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  
   \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}


\examples{
   \dontrun{
   data(SupremeCourt)
   # note that the rownames (the item names) are "1", "2", etc
   posterior1 <- MCMCirtKd(SupremeCourt, dimensions=1,
                   burnin=5000, mcmc=50000, thin=10,
                   B0=.25, store.item=TRUE,
                   item.constraints=list("1"=list(2,"-")))
   plot(posterior1)
   summary(posterior1)


   data(Senate)
   rownames(Senate) <- Senate$member
   # note that we need to transpose the data to get
   # the bills on the rows
   posterior2 <- MCMCirtKd(t(Senate[,6:677]), dimensions=2,
                   burnin=5000, mcmc=50000, thin=10,
                   item.constraints=list(rc2=list(2,"-"), rc2=c(3,0),
                                         rc3=list(3,"-")),
                   B0=.25)
   plot(posterior2)
   summary(posterior2)
   }
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
\code{\link[MCMCpack]{MCMCirt1d}}, \code{\link[MCMCpack]{MCMCordfactanal}}
}

