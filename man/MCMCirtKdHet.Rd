\name{MCMCirtKdHet}
\Rdversion{1.1}
\alias{MCMCirtKdHet}
\title{
Markov Chain Monte Carlo for Heteroskedastic K-Dimensional Item Response Theory Model}
\description{
This function generates a sample from the posterior distribution of a
  heteroskedastic K-dimensional item response theory (IRT) model, with standard
  normal priors on the subject abilities (ideal points),
  normal priors on the item parameters, and inverse-gamma priors on 
  subject error variances.  To maintain identification and comparability
  with results of the homoskedastic estimator, the mean root subject
  error precision is constrained to one.  The user supplies data and 
  priors, and a sample from the posterior distribution is returned as an mcmc 
  object, which can be subsequently analyzed with functions provided in the coda
  package.}
\usage{
MCMCirtKdHet(datamatrix, dimensions, item.constraints = list(), burnin = 1000,
mcmc = 1000, thin = 1, verbose = 0, seed = NA, alphabeta.start = NA, b0 = 0,
B0 = 0.04, c0 = 0, d0 = 0, store.item = FALSE, store.ability = TRUE,
store.sigma = TRUE, drop.constant.items = TRUE)
}
\arguments{
     \item{datamatrix}{The matrix of data.  Must be 0, 1, or NA.  
    It is of dimensionality subjects by items.}

    \item{dimensions}{The number of dimensions in the latent space.}
    
    \item{item.constraints}{List of lists specifying possible equality
    or simple inequality constraints on the item parameters. A typical
    entry in the list has one of three forms: \code{rowname=list(d,c)}
    which will constrain the dth item parameter for the item named
    rowname to be equal to c, \code{rowname=list(d,"+")} which will
    constrain the dth item parameter for the item named rowname to be
    positive, and\code{rowname=list(d, "-")} which will constrain the dth
    item parameter for the item named rowname to be negative. If x is a
    matrix without row names defaults names of ``V1", ``V2", ... , etc
    will be used. In a K dimensional model, the first item parameter for
    item \eqn{i}{i} is the difficulty parameter (\eqn{\alpha_i}{alpha_i}),
    the second item parameter is the discrimation parameter on dimension
    1 (\eqn{\beta_{i,1}}{beta_{i,1}}), the third item parameter is the
    discrimation parameter on dimension 2
    (\eqn{\beta_{i,2}}{beta_{i,2}}), ...,  and the (K+1)th item parameter
    is the discrimation parameter on dimension K
    (\eqn{\beta_{i,1}}{beta_{i,1}}). 
    The item difficulty parameters (\eqn{\alpha}{alpha}) should
    generally not be constrained. 
    }

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.   If \code{verbose} is greater
    than 0 then every \code{verbose}th iteration will be printed to the
    screen.}
    
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
    independent normal prior on the item parameters.
    Can be either a scalar or a matrix of dimension
    \eqn{(K+1) \times items}{(K+1) x items}.}
  
    \item{c0}{The first parameter of the inverse gamma prior on the 
    subject-specific variance parameters.  This can be thought of as
    the number of bills that the prior information is equivalent to.
    This scalar value is common across all subjects (legislators) and
    defaults to an uninformative prior.  NOTE: regardless of the value 
    provided, identification is provided by a constraint on the mean 
    root subject specific variance.}

    \item{d0}{The second parameter of the inverse gamma prior on the
    subject-specific variance parameters.  This can be thought of as
    the sum of square error that the prior information is equivalent to.
    This scalar value is common across all subjects (legislators) and
    defaults to an uninformative prior.  NOTE: regardless of the value 
    provided, identification is provided by a constraint on the mean 
    root subject specific variance.}

    \item{store.item}{A switch that determines whether or not to
    store the item parameters for posterior analysis. 
    \emph{NOTE: In applications with many items
      this takes an enormous amount of memory. If you have many items
    and want to want to store the item parameters you may want to thin
    the chain heavily}.  By default, the
    item parameters are not stored.}
 
   \item{store.ability}{A switch that determines whether or not to store
     the subject abilities for posterior analysis. \emph{NOTE: In
     applications with many subjects this takes an enormous amount of
     memory. If you have many subjects and want to want to store the ability
     parameters you may want to thin the chain heavily}. By default, the
     ability parameters are all stored.}

  \item{store.sigma}{A switch that determines whether or not to store
     the subject-specific variances for posterior analysis. \emph{NOTE: In
     applications with many subjects this takes an enormous amount of
     memory. If you have many subjects and want to want to store the ability
     parameters you may want to thin the chain heavily}. By default, the
     subject-specific variance parameters are all stored.}

  \item{drop.constant.items}{A switch that determines whether or not
        items that have no variation
    should be deleted before fitting the model. Default = TRUE.}
}
\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}
\references{
 Benjamin E. Lauderdale. ``Unpredictable Voters in Ideal Point Estimation''
\emph{Political Analysis.} forthcoming.
}
\author{
Benjamin E. Lauderdale, \email{blauderd@princeton.edu},
\url{http://www.princeton.edu/~blauderd/}.

Modified from \code{\link[MCMCpack]{MCMCirtKd}} and \code{\link[MCMCpack]{MCMCordfactanal}}.  Suggestions for additional options are welcome.
}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
\code{\link[MCMCpack]{MCMCirtKd}}
}

\examples{

data(Senate)
Y <- as.matrix(Senate[,6:677])

Hompost <- MCMCirtKd(Y,1,b0=0,B0=0.04,burn=1000,mcmc=1000,thin=1,verbose=250)
Hetpost <- MCMCirtKdHet(Y,1,b0=0,B0=0.04,burn=1000,mcmc=1000,thin=1,verbose=250)

SenatorNames <- Senate[,5]
HomoskedasticIdealPointEstimates <- colMeans(Hompost)[1:102]
HeteroskedasticIdealPointEstimates <- colMeans(Hetpost)[1:102]
HeteroskedasticSigmaEstimates <- colMeans(Hetpost)[103:204]

plot(HomoskedasticIdealPointEstimates, HeteroskedasticIdealPointEstimates,
cex= HeteroskedasticSigmaEstimates,xlab="Ideal Points (Homoskedastic)",
ylab="Ideal Points (Heteroskedastic)",
main="Comparison of Ideal Point Estimates for the 106th Senate",
xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))
for (i in 1:102){
	if (rank(-HeteroskedasticSigmaEstimates)[i] <= 10){
		text(HomoskedasticIdealPointEstimates[i], 
		HeteroskedasticIdealPointEstimates[i],SenatorNames[i],
		pos=3-sign(HomoskedasticIdealPointEstimates[i]),cex=0.75)
	}
}
legend(x="topleft",legend=c("Point sizes proportional to estimated legislator",
"variance under heteroskedastic model.","Some legislators with large variance have",
"more extreme estimated ideal points under the","heteroskedastic model because their",
"deviations from the party line are attributable","to idiosyncrasy rather than moderation."),cex=0.5)

}

\keyword{models}
