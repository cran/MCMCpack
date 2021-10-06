##########################################################################
## sample from the posterior distribution of a one-dimensional item
## response theory model in R using linked C++ code in Scythe.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## ADM and KQ 1/23/2003
## updated extensively ADM & KQ 7/28/2004
## store.ability arg added KQ 1/27/2006
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for One Dimensional Item Response Theory Model
#'
#' This function generates a sample from the posterior distribution of a one
#' dimensional item response theory (IRT) model, with Normal priors on the
#' subject abilities (ideal points), and multivariate Normal priors on the item
#' parameters.  The user supplies data and priors, and a sample from the
#' posterior distribution is returned as an mcmc object, which can be
#' subsequently analyzed with functions provided in the coda package.
#'
#' If you are interested in fitting K-dimensional item response theory models,
#' or would rather identify the model by placing constraints on the item
#' parameters, please see \code{\link[MCMCpack]{MCMCirtKd}}.
#'
#' \code{MCMCirt1d} simulates from the posterior distribution using standard
#' Gibbs sampling using data augmentation (a Normal draw for the subject
#' abilities, a multivariate Normal draw for the item parameters, and a
#' truncated Normal draw for the latent utilities). The simulation proper is
#' done in compiled C++ code to maximize efficiency.  Please consult the coda
#' documentation for a comprehensive list of functions that can be used to
#' analyze the posterior sample.
#'
#' The model takes the following form.  We assume that each subject has an
#' subject ability (ideal point) denoted \eqn{\theta_j} and that each
#' item has a difficulty parameter \eqn{\alpha_i} and discrimination
#' parameter \eqn{\beta_i}.  The observed choice by subject \eqn{j}
#' on item \eqn{i} is the observed data matrix which is \eqn{(I \times
#' J)}.  We assume that the choice is dictated by an unobserved
#' utility:
#'
#' \deqn{z_{i,j} = -\alpha_i + \beta_i \theta_j + \varepsilon_{i,j}}
#'
#'
#' Where the errors are assumed to be distributed standard Normal.
#' The parameters of interest are the subject abilities (ideal points)
#' and the item parameters.
#'
#' We assume the following priors.  For the subject abilities (ideal points):
#'
#' \deqn{\theta_j \sim \mathcal{N}(t_{0},T_{0}^{-1})}
#'
#' For the item parameters, the prior is:
#'
#' \deqn{\left[\alpha_i, \beta_i \right]' \sim \mathcal{N}_2 (ab_{0},AB_{0}^{-1})}
#'
#' The model is identified by the proper priors on the item parameters and
#' constraints placed on the ability parameters.
#'
#' As is the case with all measurement models, make sure that you have plenty
#' of free memory, especially when storing the item parameters.
#'
#' @param datamatrix The matrix of data.  Must be 0, 1, or missing values.  The
#' rows of \code{datamatrix} correspond to subjects and the columns correspond
#' to items.
#'
#' @param theta.constraints A list specifying possible simple equality or
#' inequality constraints on the ability parameters. A typical entry in the
#' list has one of three forms: \code{varname=c} which will constrain the
#' ability parameter for the subject named \code{varname} to be equal to c,
#' \code{varname="+"} which will constrain the ability parameter for the
#' subject named \code{varname} to be positive, and \code{varname="-"} which
#' will constrain the ability parameter for the subject named \code{varname} to
#' be negative. If x is a matrix without row names defaults names of
#' ``V1",``V2", ... , etc will be used. See Rivers (2003) for a thorough
#' discussion of identification of IRT models.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of Gibbs iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' Gibbs iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 then
#' every \code{verbose}th iteration will be printed to the screen.
#'
#' @param seed The seed for the random number generator.  If NA, the Mersenne
#' Twister generator is used with default seed 12345; if an integer is passed
#' it is used to seed the Mersenne twister.  The user can also pass a list of
#' length two to use the L'Ecuyer random number generator, which is suitable
#' for parallel computation.  The first element of the list is the L'Ecuyer
#' seed, which is a vector of length six or NA (if NA a default seed of
#' \code{rep(12345,6)} is used).  The second element of list is a positive
#' substream number. See the MCMCpack specification for more details.
#'
#' @param theta.start The starting values for the subject abilities (ideal
#' points). This can either be a scalar or a column vector with dimension equal
#' to the number of voters.  If this takes a scalar value, then that value will
#' serve as the starting value for all of the thetas.  The default value of NA
#' will choose the starting values based on an eigenvalue-eigenvector
#' decomposition of the aggreement score matrix formed from the
#' \code{datamatrix}.
#'
#' @param alpha.start The starting values for the \eqn{\alpha}
#' difficulty parameters. This can either be a scalar or a column vector with
#' dimension equal to the number of items.  If this takes a scalar value, then
#' that value will serve as the starting value for all of the alphas.  The
#' default value of NA will set the starting values based on a series of probit
#' regressions that condition on the starting values of theta.
#'
#' @param beta.start The starting values for the \eqn{\beta}
#' discrimination parameters. This can either be a scalar or a column vector
#' with dimension equal to the number of items.  If this takes a scalar value,
#' then that value will serve as the starting value for all of the betas.  The
#' default value of NA will set the starting values based on a series of probit
#' regressions that condition on the starting values of theta.
#'
#' @param t0 A scalar parameter giving the prior mean of the subject abilities
#' (ideal points).
#'
#' @param T0 A scalar parameter giving the prior precision (inverse variance)
#' of the subject abilities (ideal points).
#'
#' @param ab0 The prior mean of \code{(alpha, beta)}. Can be either a scalar or
#' a 2-vector. If a scalar both means will be set to the passed value. The
#' prior mean is assumed to be the same across all items.
#'
#' @param AB0 The prior precision of \code{(alpha, beta)}.This can either be
#' ascalar or a 2 by 2 matrix. If this takes a scalar value, then that value
#' times an identity matrix serves as the prior precision. The prior precision
#' is assumed to be the same across all items.
#'
#' @param store.item A switch that determines whether or not to store the item
#' parameters for posterior analysis.  \emph{NOTE: In situations with many
#' items storing the item parameters takes an enormous amount of memory, so
#' \code{store.item} should only be \code{FALSE} if the chain is thinned
#' heavily, or for applications with a small number of items}.  By default, the
#' item parameters are not stored.
#'
#' @param store.ability A switch that determines whether or not to store the
#' ability parameters for posterior analysis.  \emph{NOTE: In situations with
#' many individuals storing the ability parameters takes an enormous amount of
#' memory, so \code{store.ability} should only be \code{TRUE} if the chain is
#' thinned heavily, or for applications with a small number of individuals}.
#' By default, the item parameters are stored.
#'
#' @param drop.constant.items A switch that determines whether or not items
#' that have no variation should be deleted before fitting the model. Default =
#' TRUE.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the sample from the posterior
#' distribution. This object can be summarized by functions provided by the
#' coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MCMCpack]{MCMCirtKd}}
#'
#' @references James H. Albert. 1992. ``Bayesian Estimation of Normal Ogive
#' Item Response Curves Using Gibbs Sampling." \emph{Journal of Educational
#' Statistics}.  17: 251-269.
#'
#' Joshua Clinton, Simon Jackman, and Douglas Rivers. 2004. ``The Statistical
#' Analysis of Roll Call Data."  \emph{American Political Science Review}.  98:
#' 355-370.
#'
#' Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling."
#' Springer: New York.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' Douglas Rivers.  2004.  ``Identification of Multidimensional Item-Response
#' Models."  Stanford University, typescript.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
#'    ## US Supreme Court Example with inequality constraints
#'    data(SupremeCourt)
#'    posterior1 <- MCMCirt1d(t(SupremeCourt),
#'                    theta.constraints=list(Scalia="+", Ginsburg="-"),
#'                    B0.alpha=.2, B0.beta=.2,
#'                    burnin=500, mcmc=100000, thin=20, verbose=500,
#'                    store.item=TRUE)
#'    geweke.diag(posterior1)
#'    plot(posterior1)
#'    summary(posterior1)
#'
#'    ## US Senate Example with equality constraints
#'    data(Senate)
#'    Sen.rollcalls <- Senate[,6:677]
#'    posterior2 <- MCMCirt1d(Sen.rollcalls,
#'                     theta.constraints=list(KENNEDY=-2, HELMS=2),
#'                     burnin=2000, mcmc=100000, thin=20, verbose=500)
#'    geweke.diag(posterior2)
#'    plot(posterior2)
#'    summary(posterior2)
#'    }
#'
"MCMCirt1d" <-
  function(datamatrix, theta.constraints=list(), burnin = 1000,
           mcmc = 20000, thin=1, verbose = 0, seed = NA,
           theta.start = NA, alpha.start = NA, beta.start = NA, t0 = 0,
           T0 = 1, ab0=0, AB0=.25, store.item = FALSE, store.ability=TRUE,
           drop.constant.items=TRUE, ... ) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    ## check vote matrix and convert to work with C++ code
    if (drop.constant.items==TRUE){
      x.col.var <- apply(datamatrix, 2, var, na.rm=TRUE)
      keep.inds <- x.col.var>0
      keep.inds[is.na(keep.inds)] <- FALSE
      datamatrix <- datamatrix[,keep.inds]
    }
    datamatrix <- as.matrix(datamatrix)
    K <- ncol(datamatrix)   # cases, bills, items, etc
    J <- nrow(datamatrix)   # justices, legislators, subjects, etc
    if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) != (J*K)) {
      cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
      stop("Please check data and try MCMCirt1d() again.\n",
              call.=FALSE)
    }
    datamatrix[is.na(datamatrix)] <- 9
    item.names <- colnames(datamatrix)
    subject.names <- rownames(datamatrix)

    ## setup constraints on theta
    if(length(theta.constraints) != 0) {
       for (i in 1:length(theta.constraints)){
          theta.constraints[[i]] <-
             list(as.integer(1), theta.constraints[[i]][1])
       }
    }
    holder <- build.factor.constraints(theta.constraints, t(datamatrix), J, 1)
    theta.eq.constraints <- holder[[1]]
    theta.ineq.constraints <- holder[[2]]
    subject.names <- holder[[3]]
    ## names
    item.names <- colnames(datamatrix)
    if (is.null(item.names)){
      item.names <- paste("item", 1:K, sep="")
    }

    ## prior for theta
    holder <- form.mvn.prior(t0, T0, 1)
    t0 <- holder[[1]]
    T0 <- holder[[2]]

    ## prior for (alpha, beta)
    holder <- form.mvn.prior(ab0, AB0, 2)
    ab0 <- holder[[1]]
    AB0 <- holder[[2]]

    ## starting values for theta error checking
    theta.start <- factor.score.start.check(theta.start, datamatrix,
                                            t0, T0,
                                            theta.eq.constraints,
                                            theta.ineq.constraints, 1)

    ## starting values for (alpha, beta)
    ab.starts <- matrix(NA, K, 2)
    for (i in 1:K){
      local.y <-  datamatrix[,i]
      local.y[local.y==9] <- NA
      if (var(na.omit(local.y))==0){
        ab.starts[i,] <- c(0,10)
      }
      else {
        ab.starts[i,] <- coef(suppressWarnings(glm(local.y~theta.start,
                                                   family=binomial(link="probit"),
                                                   control=glm.control(
                                                     maxit=8, epsilon=1e-3)
                                                   )))
      }
    }
    ab.starts[,1] <- -1 * ab.starts[,1] # make this into a difficulty param

    ## starting values for alpha and beta error checking
    if (is.na(alpha.start)) {
      alpha.start <- ab.starts[,1]
    }
    else if(is.null(dim(alpha.start))) {
      alpha.start <- alpha.start * matrix(1,K,1)
    }
    else if((dim(alpha.start)[1] != K) || (dim(alpha.start)[2] != 1)) {
      cat("Error: Starting value for alpha not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n",
           call.=FALSE)
    }
    if (is.na(beta.start)) {
      beta.start <- ab.starts[,2]
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,K,1)
    }
    else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
      cat("Error: Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n",
              call.=FALSE)
    }

    ## define holder for posterior sample
    if(store.item == FALSE & store.ability == TRUE) {
      sample <- matrix(data=0, mcmc/thin, J)
    }
    else if (store.item == TRUE & store.ability == FALSE){
      sample <- matrix(data=0, mcmc/thin, 2*K)
    }
    else if (store.item == TRUE & store.ability == TRUE){
      sample <- matrix(data=0, mcmc/thin, J + 2 * K)
    }
    else{
      cat("Error: store.item == FALSE & store.ability == FALSE.\n")
      stop("Please respecify and call MCMCirt1d() again.\n",
              call.=FALSE)
    }

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    # call C++ code to draw sample
    posterior <- .C("cMCMCirt1d",
                    sampledata = as.double(sample),
                    samplerow = as.integer(nrow(sample)),
                    samplecol = as.integer(ncol(sample)),
                    Xdata = as.integer(datamatrix),
                    Xrow = as.integer(nrow(datamatrix)),
                    Xcol = as.integer(ncol(datamatrix)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    thetastartdata = as.double(theta.start),
                    thetastartrow = as.integer(nrow(theta.start)),
                    thetastartcol = as.integer(ncol(theta.start)),
                    astartdata = as.double(alpha.start),
                    astartrow = as.integer(length(alpha.start)),
                    astartcol = as.integer(1),
                    bstartdata = as.double(beta.start),
                    bstartrow = as.integer(length(beta.start)),
                    bstartcol = as.integer(1),
                    t0 = as.double(t0),
                    T0 = as.double(T0),
                    ab0data = as.double(ab0),
                    ab0row = as.integer(nrow(ab0)),
                    ab0col = as.integer(ncol(ab0)),
                    AB0data = as.double(AB0),
                    AB0row = as.integer(nrow(AB0)),
                    AB0col = as.integer(ncol(AB0)),
                    thetaeqdata = as.double(theta.eq.constraints),
                    thetaeqrow = as.integer(nrow(theta.eq.constraints)),
                    thetaeqcol = as.integer(ncol(theta.eq.constraints)),
                    thetaineqdata = as.double(theta.ineq.constraints),
                    thetaineqrow = as.integer(nrow(theta.ineq.constraints)),
                    thetaineqcol = as.integer(ncol(theta.ineq.constraints)),
                    storei = as.integer(store.item),
                    storea = as.integer(store.ability),
                    PACKAGE="MCMCpack"
                  )

    theta.names <- paste("theta.", subject.names, sep = "")
    alpha.beta.names <- paste(rep(c("alpha.","beta."), K),
                              rep(item.names, each = 2),
                              sep = "")

    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol,
                     byrow=FALSE)
    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)

    names <- NULL
    if(store.ability == TRUE) {
      names <- c(names, theta.names)
    }
    if (store.item == TRUE){
      names <- c(names, alpha.beta.names)
    }
    varnames(output) <- names
    attr(output,"title") <-
      "MCMCirt1d Posterior Sample"
    return(output)

  }
