\name{MCMCfactanal}
\alias{MCMCfactanal}
\title{Markov chain Monte Carlo for Normal Theory Factor Analysis Model}
\description{
  This function generates a posterior density sample from Normal theory
  factor analysis model. Normal priors are assumed on the factor
  loadings and factor scores while inverse Gamma priors are assumed for
  the uniquenesses. The user supplies data and parameters for the prior
  distributions, and a sample from the posterior density is returned as
  an mcmc object, which can be subsequently analyzed with
  functions provided in the coda package.
}
  
\usage{
MCMCfactanal(x, factors, lambda.constraints=list(),
             data=list(), burnin = 1000, mcmc = 10000,
             thin=5, verbose = FALSE, seed = 0,
             lambda.start = NA, psi.start = NA,
             l0=0, L0=0, a0=0.001, b0=0.001,
             store.scores = FALSE, std.var=TRUE, ... )
 }

\arguments{
  \item{x}{Either a formula or a numeric matrix containing the
    manifest variables.}

  \item{factors}{The number of factors to be fitted.}
  
  \item{lambda.constraints}{List of lists specifying possible simple equality
    or inequality constraints on the factor loadings. A typical
    entry in the list has one of three forms: \code{varname=list(d,c)} which
    will constrain the dth loading for the variable named varname to
    be equal to c, \code{varname=list(d,"+")} which will constrain the dth
    loading for the variable named varname to be positive, and
    \code{varname=list(d, "-")} which will constrain the dth loading for the
    variable named varname to be negative. If x is a matrix without
    column names defaults names of ``V1",``V2", ... , etc will be
    used.}
  
  \item{data}{A data frame.}
  
  \item{burnin}{The number of burn-in iterations for the sampler.}
  
  \item{mcmc}{The number of iterations for the sampler.}
  
  \item{thin}{The thinning interval used in the simulation.  The number of
    iterations must be divisible by this value.}
  
  \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If TRUE, the iteration number and
    the factor loadings and uniquenesses are printed to the screen.}
  
  \item{seed}{The seed for the random number generator.  The code uses the
    Mersenne Twister, which requires an integer as an input.  If nothing
    is provided, the Scythe default seed is used.}
  
  \item{lambda.start}{Starting values for the factor loading matrix
    Lambda. If \code{lambda.start} is set to a scalar the starting value for
    all unconstrained loadings will be set to that scalar. If
    \code{lambda.start} is a matrix of the same dimensions as Lambda then the
    \code{lambda.start} matrix is used as the starting values (except
    for equality-constrained elements). If \code{lambda.start} is set to
    \code{NA} (the default) then starting values for unconstrained
    elements are set to 0, and starting values for inequality
    constrained elements are set to either 0.5 or -0.5 depending on the
    nature of the constraints.}
  
  \item{psi.start}{Starting values for the uniquenesses. If
    \code{psi.start} is set to a scalar then the starting value for all
    diagonal elements of \code{Psi} are set to this value. If
    \code{psi.start} is a \eqn{k}{k}-vector (where \eqn{k}{k} is the
    number of manifest variables) then the staring value of \code{Psi}
    has \code{psi.start} on the main diagonal. If \code{psi.start} is
    set to \code{NA} (the default) the starting values of all the
    uniquenesses are set to 0.5.}

  \item{l0}{The means of the independent Normal prior on the factor
    loadings. Can be either a scalar or a matrix with the same
    dimensions as \code{Lambda}.}

  \item{L0}{The precisions (inverse variances) of the independent Normal
  prior on the factor loadings. Can be either a scalar or a matrix with
  the same dimensions as \code{Lambda}.}

  \item{a0}{Controls the shape of the inverse Gamma prior on the
    uniqueness. Can be either a scalar or a \eqn{k}{k}-vector.}

  \item{b0}{Controls the scale of the inverse Gamma prior on the
    uniquenesses. Can be either a scalar or a \eqn{k}{k}-vector.}
  
  \item{store.scores}{A switch that determines whether or not to
    store the factor scores for posterior analysis. 
    \emph{NOTE: This takes an enormous amount of memory, so
      should only be used if the chain is thinned heavily, or for
      applications with a small number of observations}.  By default, the
    factor scores are not stored.}

  \item{std.var}{If \code{TRUE} (the default) the manifest variables are
  rescaled to have zero mean and unit variance. Otherwise, the manifest
  variables are rescaled to have zero mean but retain their observed
  variances.}
  
  \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{The model takes the following form:

  \deqn{x_i = \Lambda \phi_i + \epsilon_i}{x_i = Lambda phi_i + epsilon_i}
  \deqn{\epsilon_i \sim \mathcal{N}(0,\Psi)}{epsilon_i ~ N(0, Psi)}

  where \eqn{x_i}{x_i} is the \eqn{k}{k}-vector of observed variables
  specific to observation \eqn{i}{i}, \eqn{\Lambda}{Lambda} is the
  \eqn{k \times d}{k by d} matrix of factor loadings, \eqn{\phi_i}{phi_i} is
  the \eqn{d}{d}-vector of latent factor scores, and \eqn{\Psi}{Psi} is
  a diagonal, positive definite matrix. Traditional factor analysis
  texts refer to the diagonal elements of \eqn{\Psi}{Psi} as
  uniquenesses.  

  The implementation used here assumes independent conjugate priors for
  each element of \eqn{\Lambda}{Lambda}, each \eqn{\phi_i}{phi_i}, and
  each diagonal element of \eqn{\Psi}{Psi}. More specifically we assume:

  \deqn{\Lambda_{ij} \sim \mathcal{N}(l_{0_{ij}}, L_{0_{ij}}^{-1}),
    i=1,\ldots,k,  j=1,\ldots,d}{Lambda_ij ~ N(l0_ij,  L0_ij^-1),
    i=1,...,k, j=1,...,d} 

    \deqn{\phi_i \sim \mathcal{N}(0, I),  i=1,\dots,n}{phi_i ~ N(0, I),
      i=1,...,n} 

    \deqn{\Psi_{ii} \sim \mathcal{IG}(a_{0_i}/2, b_{0_i}/2), 
      i=1,\ldots,k}{Psi_ii ~ IG(a0_i/2, b0_i/2), i=1,...,k}
    
  \code{MCMCfactanal} simulates from the posterior density using
  standard Gibbs sampling. The simulation proper is done in
  compiled C++ code to maximize efficiency.  Please consult the
  coda documentation for a comprehensive list of functions that
  can be used to analyze the posterior density sample.     
  }
}

\references{
   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2003.  
   \emph{Scythe Statistical Library 0.4.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}


\examples{
   \dontrun{
   ### An example using the formula interface
   data(swiss)
   posterior <- MCMCfactanal(~Agriculture+Examination+Education+Catholic
                             +Infant.Mortality, factors=2,
                             lambda.constraints=list(Examination=list(1,"+"),
                             Examination=list(2,"-"), Education=c(2,0),
                             Infant.Mortality=c(1,0)),
                             verbose=FALSE, store.scores=FALSE, a0=1, b0=0.15,
                             data=swiss, burnin=5000, mcmc=50000, thin=20)
   plot(posterior)
   summary(posterior)

   ### An example using the matrix interface
   Lambda <- matrix(runif(45,-.5,.5), 15, 3)
   Psi <- diag(1 - apply(Lambda ^2, 1, sum))
   Sigma <- Lambda \%*\% t(Lambda) + Psi 
   Y <- t(t(chol(Sigma)) \%*\% matrix(rnorm(500*15), 15, 500))

   posterior <- MCMCfactanal(Y, factors=3,
                             lambda.constraints=list(V1=c(2,0),
                             V1=c(3,0), V2=c(3,0), V3=list(1,"+"),
                             V3=list(2,"+"), V3=list(3,"+")),
                             verbose=FALSE)
   plot(posterior)
   summary(posterior)        
   }
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},\code{\link[mva]{factanal}}}

