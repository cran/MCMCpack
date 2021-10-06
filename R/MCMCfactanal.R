##########################################################################
## sample from the posterior distribution of a factor analysis model
## model in R using linked C++ code in Scythe.
##
## The model is:
##
## x_i = \Lambda \phi_i + \epsilon_i,   \epsilon_i \sim N(0, \Psi)
##
## where \Psi is diagonal and the priors are:
##
## \lambda_{ij} \sim N(l_{ij}, L^{-1}_{ij})
## \phi_i \sim N(0,I)
## \psi^2_{jj} \sim IG(a0_j/2, b0_j/2)
##
##
## Andrew D. Martin
## Washington University
##
## Kevin M. Quinn
## Harvard University
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## May 7, 2003
## revised to accomodate new spec 7/5/2004 KQ
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for Normal Theory Factor Analysis Model
#'
#' This function generates a sample from the posterior distribution of
#' a normal theory factor analysis model. Normal priors are assumed on
#' the factor loadings and factor scores while inverse Gamma priors
#' are assumed for the uniquenesses. The user supplies data and
#' parameters for the prior distributions, and a sample from the
#' posterior distribution is returned as an mcmc object, which can be
#' subsequently analyzed with functions provided in the coda package.
#'
#' The model takes the following form:
#'
#' \deqn{x_i = \Lambda \phi_i + \epsilon_i}
#'
#' \deqn{\epsilon_i \sim \mathcal{N}(0,\Psi)}
#'
#' where \eqn{x_i} is the \eqn{k}-vector of observed variables
#' specific to observation \eqn{i}, \eqn{\Lambda} is the \eqn{k \times
#' d} matrix of factor loadings, \eqn{\phi_i} is the \eqn{d}-vector of
#' latent factor scores, and \eqn{\Psi} is a diagonal, positive
#' definite matrix. Traditional factor analysis texts refer to the
#' diagonal elements of \eqn{\Psi} as uniquenesses.
#'
#' The implementation used here assumes independent conjugate priors
#' for each element of \eqn{\Lambda} each \eqn{\phi_i}, and each
#' diagonal element of \eqn{\Psi}. More specifically we assume:
#'
#' \deqn{\Lambda_{ij} \sim \mathcal{N}(l_{0_{ij}}, L_{0_{ij}}^{-1}),
#' i=1,\ldots,k, j=1,\ldots,d}
#'
#' \deqn{\phi_i \sim \mathcal{N}(0, I), i=1,\dots,n}
#'
#' \deqn{\Psi_{ii} \sim \mathcal{IG}(a_{0_i}/2, b_{0_i}/2),
#' i=1,\ldots,k}
#'
#' \code{MCMCfactanal} simulates from the posterior distribution using
#' standard Gibbs sampling. The simulation proper is done in compiled
#' C++ code to maximize efficiency.  Please consult the coda
#' documentation for a comprehensive list of functions that can be
#' used to analyze the posterior sample.
#'
#' As is the case with all measurement models, make sure that you have
#' plenty of free memory, especially when storing the scores.
#'
#' @param x Either a formula or a numeric matrix containing the
#'   manifest variables.
#'
#' @param factors The number of factors to be fitted.
#'
#' @param lambda.constraints List of lists specifying possible simple
#'   equality or inequality constraints on the factor loadings. A
#'   typical entry in the list has one of three forms:
#'   \code{varname=list(d,c)} which will constrain the dth loading for
#'   the variable named \code{varname} to be equal to c,
#'   \code{varname=list(d,"+")} which will constrain the dth loading
#'   for the variable named \code{varname} to be positive, and
#'   \code{varname=list(d, "-")} which will constrain the dth loading
#'   for the variable named \code{varname} to be negative. If x is a
#'   matrix without column names defaults names of ``V1",``V2", ... ,
#'   etc will be used.
#'
#' @param data A data frame.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The
#'   number of iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the
#'   progress of the sampler is printed to the screen.  If
#'   \code{verbose} is greater than 0 the iteration number and the
#'   factor loadings and uniquenesses are printed to the screen every
#'   \code{verbose}th iteration.
#'
#' @param seed The seed for the random number generator.  If NA, the
#'   Mersenne Twister generator is used with default seed 12345; if an
#'   integer is passed it is used to seed the Mersenne twister.  The
#'   user can also pass a list of length two to use the L'Ecuyer
#'   random number generator, which is suitable for parallel
#'   computation.  The first element of the list is the L'Ecuyer seed,
#'   which is a vector of length six or NA (if NA a default seed of
#'   \code{rep(12345,6)} is used).  The second element of list is a
#'   positive substream number. See the MCMCpack specification for
#'   more details.
#'
#' @param lambda.start Starting values for the factor loading matrix
#'   Lambda. If \code{lambda.start} is set to a scalar the starting
#'   value for all unconstrained loadings will be set to that
#'   scalar. If \code{lambda.start} is a matrix of the same dimensions
#'   as Lambda then the \code{lambda.start} matrix is used as the
#'   starting values (except for equality-constrained elements). If
#'   \code{lambda.start} is set to \code{NA} (the default) then
#'   starting values for unconstrained elements are set to 0, and
#'   starting values for inequality constrained elements are set to
#'   either 0.5 or -0.5 depending on the nature of the constraints.
#'
#' @param psi.start Starting values for the uniquenesses. If
#'   \code{psi.start} is set to a scalar then the starting value for
#'   all diagonal elements of \code{Psi} are set to this value. If
#'   \code{psi.start} is a \eqn{k}-vector (where \eqn{k} is the number
#'   of manifest variables) then the staring value of \code{Psi} has
#'   \code{psi.start} on the main diagonal. If \code{psi.start} is set
#'   to \code{NA} (the default) the starting values of all the
#'   uniquenesses are set to 0.5.
#'
#' @param l0 The means of the independent Normal prior on the factor loadings.
#' Can be either a scalar or a matrix with the same dimensions as
#' \code{Lambda}.
#'
#' @param L0 The precisions (inverse variances) of the independent Normal prior
#' on the factor loadings. Can be either a scalar or a matrix with the same
#' dimensions as \code{Lambda}.
#'
#' @param a0 Controls the shape of the inverse Gamma prior on the uniqueness.
#' The actual shape parameter is set to \code{a0/2}. Can be either a scalar or
#' a \eqn{k}-vector.
#'
#' @param b0 Controls the scale of the inverse Gamma prior on the uniquenesses.
#' The actual scale parameter is set to \code{b0/2}. Can be either a scalar or
#' a \eqn{k}-vector.
#'
#' @param store.scores A switch that determines whether or not to store the
#' factor scores for posterior analysis.  \emph{NOTE: This takes an enormous
#' amount of memory, so should only be used if the chain is thinned heavily, or
#' for applications with a small number of observations}.  By default, the
#' factor scores are not stored.
#'
#' @param std.var If \code{TRUE} (the default) the manifest variables are
#' rescaled to have zero mean and unit variance. Otherwise, the manifest
#' variables are rescaled to have zero mean but retain their observed
#' variances.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the sample from the posterior
#' distribution. This object can be summarized by functions provided by the
#' coda package.
#'
#' @export
#'
#' @seealso
#' \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},\code{\link[stats]{factanal}}
#'
#' @references Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#' ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical
#' Software}. 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
#'    ### An example using the formula interface
#'    data(swiss)
#'    posterior <- MCMCfactanal(~Agriculture+Examination+Education+Catholic
#'                     +Infant.Mortality, factors=2,
#'                     lambda.constraints=list(Examination=list(1,"+"),
#'                        Examination=list(2,"-"), Education=c(2,0),
#'                        Infant.Mortality=c(1,0)),
#'                     verbose=0, store.scores=FALSE, a0=1, b0=0.15,
#'                     data=swiss, burnin=5000, mcmc=50000, thin=20)
#'    plot(posterior)
#'    summary(posterior)
#'
#'    ### An example using the matrix interface
#'    Y <- cbind(swiss$Agriculture, swiss$Examination,
#'               swiss$Education, swiss$Catholic,
#'               swiss$Infant.Mortality)
#'    colnames(Y) <- c("Agriculture", "Examination", "Education", "Catholic",
#'                     "Infant.Mortality")
#'    post <- MCMCfactanal(Y, factors=2,
#'                         lambda.constraints=list(Examination=list(1,"+"),
#'                           Examination=list(2,"-"), Education=c(2,0),
#'                           Infant.Mortality=c(1,0)),
#'                         verbose=0, store.scores=FALSE, a0=1, b0=0.15,
#'                         burnin=5000, mcmc=50000, thin=20)
#'    }
#'
"MCMCfactanal" <-
  function(x, factors, lambda.constraints=list(),
           data=NULL, burnin = 1000, mcmc = 20000,
           thin=1, verbose = 0, seed = NA,
           lambda.start = NA, psi.start = NA,
           l0=0, L0=0, a0=0.001, b0=0.001,
           store.scores = FALSE, std.var=TRUE, ... ) {

    ## check for an offset
    check.offset(list(...))

    ## get data matrix and associated constants
    if (is.matrix(x)){
      X <- x
      xvars <- dimnames(X)[[2]]
      xobs <- dimnames(X)[[1]]
      N <- nrow(X)
      K <- ncol(X)
    }
    else {
      holder <- parse.formula(formula=x, data=data,
                              intercept=FALSE, justX=TRUE)
      X <- holder[[2]]
      xvars <- holder[[3]]
      xobs <- holder[[4]]
      N <- nrow(X)
      K <- ncol(X)
    }

    ## standardize X
    if (std.var){
      for (i in 1:K){
        X[,i] <- (X[,i]-mean(X[,i]))/sd(X[,i])
      }
    }
    else{
      for (i in 1:K){
        X[,i] <- X[,i]-mean(X[,i])
      }
    }

    ## take care of the case where X has no row names
    if (is.null(xobs)){
      xobs <- 1:N
    }

    ## check mcmc parameters
    check.mcmc.parameters(burnin, mcmc, thin)

    # seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## setup constraints on Lambda
    holder <- build.factor.constraints(lambda.constraints, X, K, factors)
    Lambda.eq.constraints <- holder[[1]]
    Lambda.ineq.constraints <- holder[[2]]
    X.names <- holder[[3]]

    ## setup prior on Lambda
    holder <- form.factload.norm.prior(l0, L0, K, factors, X.names)
    Lambda.prior.mean <- holder[[1]]
    Lambda.prior.prec <- holder[[2]]

    ## setup and check prior on Psi
    holder <- form.ig.diagmat.prior(a0, b0, K)
    a0 <- holder[[1]]
    b0 <- holder[[2]]

    ## starting values for Lambda
    Lambda <- factload.start(lambda.start, K, factors,
                             Lambda.eq.constraints, Lambda.ineq.constraints)

    ## starting values for Psi
    Psi <- factuniqueness.start(psi.start, X)


    ## define holder for posterior sample
    if(store.scores == FALSE) {
      sample <- matrix(data=0, mcmc/thin, K*factors+K)
    }
    else {
      sample <- matrix(data=0, mcmc/thin, K*factors+K+factors*N)
    }
    posterior <- NULL

    ## call C++ code to do the sampling
    auto.Scythe.call(output.object="posterior", cc.fun.name="cMCMCfactanal",
                     sample.nonconst=sample, X=X, burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose),
                     Lambda=Lambda, Psi=Psi, Lameq=Lambda.eq.constraints,
                     Lamineq=Lambda.ineq.constraints,
                     Lampmean=Lambda.prior.mean, Lampprec=Lambda.prior.prec,
                     a0=a0, b0=b0, storescores=as.integer(store.scores))

    Lambda.names <- paste(paste("Lambda",
                                rep(X.names,
                                    each=factors), sep=""),
                          rep(1:factors,K), sep="_")
    Psi.names <- paste("Psi", X.names, sep="")
    par.names <- c(Lambda.names, Psi.names)
    if (store.scores==TRUE){
      phi.names <- paste(paste("phi",
                               rep(xobs, each=factors), sep="_"),
                         rep(1:factors,factors), sep="_")
      par.names <- c(par.names, phi.names)
    }
    title <- "MCMCpack Factor Analysis Posterior Sample"
    output <- form.mcmc.object(posterior, par.names, title)
    ## get rid of columns for constrained parameters
    output.df <- as.data.frame(as.matrix(output))
    output.sd <- apply(output.df, 2, sd)
    output.df <- output.df[,output.sd != 0]

    output <- mcmc(as.matrix(output.df), start=burnin+1, end=mcmc+burnin,
                   thin=thin)

    ## add constraint info so this isn't lost
    attr(output, "constraints") <- lambda.constraints
    attr(output, "n.manifest") <- K
    attr(output, "n.factors") <- factors
    return(output)
  }
