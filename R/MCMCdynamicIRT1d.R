##########################################################################
## samples from the posterior distribution of a dynamic 1d IRT model
## a la Martin and Quinn (2002)
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Assumes a local level model for the evolution of theta
##
##  y_{jkt}^* = -alpha_k + beta_k * theta_{jt} + epsilon_{jkt}
##  theta_{jt} ~ N(theta_{j(t-1)}, tau^2)
##
## Kevin Quinn
## 1/28/2008
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for Dynamic One Dimensional Item Response
#' Theory Model
#'
#' This function generates a sample from the posterior distribution of
#' a dynamic one dimensional item response theory (IRT) model, with
#' Normal random walk priors on the subject abilities (ideal points),
#' and multivariate Normal priors on the item parameters. The user
#' supplies data and priors, and a sample from the posterior
#' distribution is returned as an mcmc object, which can be
#' subsequently analyzed with functions provided in the coda package.
#'
#' \code{MCMCdynamicIRT1d} simulates from the posterior distribution
#' using the algorithm of Martin and Quinn (2002). The simulation
#' proper is done in compiled C++ code to maximize efficiency.  Please
#' consult the coda documentation for a comprehensive list of
#' functions that can be used to analyze the posterior sample.
#'
#' The model takes the following form. We assume that each subject has
#' an subject ability (ideal point) denoted \eqn{\theta_{j,t}} (where
#' \eqn{j} indexes subjects and \eqn{t} indexes time periods) and that
#' each item has a difficulty parameter \eqn{\alpha_i} and
#' discrimination parameter \eqn{\beta_i}. The observed choice by
#' subject \eqn{j} on item \eqn{i} is the observed data matrix which
#' is \eqn{(I \times J)}. We assume that the choice is dictated by an
#' unobserved utility:
#'
#' \deqn{z_{i,j,t} = -\alpha_i + \beta_i \theta_{j,t} +
#' \varepsilon_{i,j,t}}
#'
#' Where the disturbances are assumed to be distributed standard
#' Normal. The parameters of interest are the subject abilities (ideal
#' points) and the item parameters.
#'
#' We assume the following priors.  For the subject abilities (ideal
#' points):
#'
#' \deqn{\theta_{j,t} \sim \mathcal{N}(\theta_{j,t-1}, \tau^2_j)}
#'
#' with
#'
#' \deqn{\theta_{j,0} \sim \mathcal{N}(e0, E0)}.
#'
#' The evolution variance has the following prior:
#'
#' \deqn{\tau^2_j \sim \mathcal{IG}(c0/2, d0/2)}.
#'
#' For the item parameters in the standard model, the prior is:
#'
#' \deqn{\alpha_i \sim \mathcal{N}(a0, A0^{-1})}
#'
#' and
#'
#' \deqn{\beta_i \sim \mathcal{N}(b0, B0^{-1})}.
#'
#' The model is identified by the proper priors on the item parameters
#' and constraints placed on the ability parameters.
#'
#' As is the case with all measurement models, make sure that you have
#' plenty of free memory, especially when storing the item parameters.
#'
#' @aliases MCMCdynamicIRT1d MCMCdynamicIRT1d_b
#'
#' @param datamatrix The matrix of data.  Must be 0, 1, or missing
#'   values.  The rows of \code{datamatrix} correspond to subjects and
#'   the columns correspond to items.
#'
#' @param item.time.map A vector that relates each item to a time
#'   period.  Each element of \code{item.time.map} gives the time
#'   period of the corresponding column of \code{datamatrix}. It is
#'   assumed that the minimum value of \code{item.time.map} is 1.
#'
#' @param theta.constraints A list specifying possible simple equality
#'   or inequality constraints on the ability parameters. A typical
#'   entry in the list has one of three forms: \code{varname=c} which
#'   will constrain the ability parameter for the subject named
#'   \code{varname} to be equal to c, \code{varname="+"} which will
#'   constrain the ability parameter for the subject named
#'   \code{varname} to be positive, and \code{varname="-"} which will
#'   constrain the ability parameter for the subject named
#'   \code{varname} to be negative. If x is a matrix without row names
#'   defaults names of ``V1",``V2", ... , etc will be used. See Rivers
#'   (2003) for a thorough discussion of identification of IRT models.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of Gibbs iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The
#'   number of Gibbs iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the
#'   progress of the sampler is printed to the screen.  If
#'   \code{verbose} is greater than 0 then every \code{verbose}th
#'   iteration will be printed to the screen.
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
#' @param theta.start The starting values for the subject abilities
#'   (ideal points). This can either be a scalar or a column vector
#'   with dimension equal to the number of voters.  If this takes a
#'   scalar value, then that value will serve as the starting value
#'   for all of the thetas.  The default value of NA will choose the
#'   starting values based on an eigenvalue-eigenvector decomposition
#'   of the aggreement score matrix formed from the \code{datamatrix}.
#'
#' @param alpha.start The starting values for the \eqn{\alpha}
#'   difficulty parameters. This can either be a scalar or a column
#'   vector with dimension equal to the number of items.  If this
#'   takes a scalar value, then that value will serve as the starting
#'   value for all of the alphas.  The default value of NA will set
#'   the starting values based on a series of probit regressions that
#'   condition on the starting values of theta.
#'
#' @param beta.start The starting values for the \eqn{\beta}
#'   discrimination parameters. This can either be a scalar or a
#'   column vector with dimension equal to the number of items.  If
#'   this takes a scalar value, then that value will serve as the
#'   starting value for all of the betas.  The default value of NA
#'   will set the starting values based on a series of probit
#'   regressions that condition on the starting values of theta.
#'
#' @param tau2.start The starting values for the evolution variances
#'   (the variance of the random walk increments for the ability
#'   parameters / ideal points. Order corresponds to the rows of
#'   \code{datamatrix}.
#'
#' @param a0 A vector containing the prior mean of each of the
#'   difficulty parameters \eqn{\alpha}. Should have as many
#'   elements as items / roll calls. Order corresponds to the columns
#'   of \code{datamatrix}. If a scalar is passed it is assumed that
#'   all elements of \code{a0} are equal to the scalar.
#'
#' @param A0 A vector containing the prior precision (inverse
#'   variance) of each of the difficulty parameters \eqn{\alpha}.
#'   Should have as many elements as items / roll calls. Order
#'   corresponds to the columns of \code{datamatrix}. If a scalar is
#'   passed it is assumed that all elements of \code{A0} are equal to
#'   the scalar.
#'
#' @param b0 A vector containing the prior mean of each of the
#'   discrimination parameters \eqn{\beta}. Should have as many
#'   elements as items / roll calls. Order corresponds to the columns
#'   of \code{datamatrix}. If a scalar is passed it is assumed that
#'   all elements of \code{b0} are equal to the scalar.
#'
#' @param B0 A vector containing the prior precision (inverse
#'   variance) of each of the discrimination parameters
#'   \eqn{\beta}. Should have as many elements as items / roll
#'   calls. Order corresponds to the columns of \code{datamatrix}. If
#'   a scalar is passed it is assumed that all elements of \code{B0}
#'   are equal to the scalar.
#'
#' @param c0 \eqn{c_{0/2}} is the shape parameter for the inverse
#'   Gamma prior on \eqn{\tau^2} (the variance of the random walk
#'   increments).  The amount of information in the inverse Gamma
#'   prior is something like that from \eqn{c_0}
#'   pseudo-observations. \code{c0} can be either a vector with an
#'   element for each subject or a scalar. If \code{c0} is negative
#'   then \eqn{\tau^2} is not estimated-- the values in
#'   \code{tau2.start} are used throughout the sampling.
#'
#' @param d0 \eqn{d_{0/2}} is the scale parameter for the inverse
#'   Gamma prior on \eqn{\tau^2} (the variance of the random walk
#'   increments).  In constructing the inverse Gamma prior, \eqn{d_0}
#'   acts like the sum of squared errors from the \eqn{c_0}
#'   pseudo-observations. \code{d0} can be either a vector with an
#'   element for each subject or a scalar. If \code{d0} is negative
#'   then \eqn{\tau^2} is not estimated-- the values in
#'   \code{tau2.start} are used throughout the sampling.
#'
#' @param e0 A vector containing the prior mean of the initial ability
#'   parameter / ideal point for each subject. Should have as many
#'   elements as subjects. Order corresponds to the rows of
#'   \code{datamatrix}. If a scalar is passed it is assumed that all
#'   elements of \code{e0} are equal to the scalar.
#'
#' @param E0 A vector containing the prior variance of the initial
#'   ability parameter / ideal point for each subject. Should have as
#'   many elements as subjects. Order corresponds to the rows of
#'   \code{datamatrix}. If a scalar is passed it is assumed that all
#'   elements of \code{E0} are equal to the scalar.
#'
#' @param store.ability A switch that determines whether or not to
#'   store the ability parameters for posterior analysis.
#'   \emph{NOTE}: In situations with many individuals storing the
#'   ability parameters takes an enormous amount of memory, so
#'   \code{store.ability} should only be \code{TRUE} if the chain is
#'   thinned heavily, or for applications with a small number of
#'   individuals.  By default, the item parameters are stored.
#'
#' @param store.item A switch that determines whether or not to store
#'   the item parameters for posterior analysis.  \emph{NOTE}: In
#'   situations with many items storing the item parameters takes an
#'   enormous amount of memory, so \code{store.item} should only be
#'   \code{FALSE} if the chain is thinned heavily, or for applications
#'   with a small number of items.  By default, the item parameters
#'   are not stored.
#'
#' @param \dots further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample. This
#'   object can be summarized by functions provided by the coda
#'   package.
#'
#' @export
#'
#' @author Kevin M. Quinn
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MCMCpack]{MCMCirt1d}}
#'
#' @references Andrew D. Martin and Kevin M. Quinn. 2002. "Dynamic Ideal Point
#' Estimation via Markov Chain Monte Carlo for the U.S. Supreme Court,
#' 1953-1999." \emph{Political Analysis.} 10: 134-153. <doi:10.1093/pan/10.2.134>
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' @keywords models
#'
#' @examples
#'
#'   \dontrun{
#' 	data(Rehnquist)
#'
#' 	## assign starting values
#' 	theta.start <- rep(0, 9)
#' 	theta.start[2] <- -3 ## Stevens
#' 	theta.start[7] <- 2  ## Thomas
#'
#' 	out <- MCMCdynamicIRT1d(t(Rehnquist[,1:9]),
#' 	                        item.time.map=Rehnquist$time,
#' 	                        theta.start=theta.start,
#' 	                        mcmc=50000, burnin=20000, thin=5,
#' 	                        verbose=500, tau2.start=rep(0.1, 9),
#' 	                        e0=0, E0=1,
#' 	                        a0=0, A0=1,
#' 	                        b0=0, B0=1, c0=-1, d0=-1,
#' 	                        store.item=FALSE,
#' 	                        theta.constraints=list(Stevens="-", Thomas="+"))
#'
#' 	summary(out)
#'   }
#'
"MCMCdynamicIRT1d" <- function(datamatrix, item.time.map,
                               theta.constraints=list(),
                               burnin=1000, mcmc=20000, thin=1,
                               verbose=0, seed=NA, theta.start=NA,
                               alpha.start=NA, beta.start=NA,
                               tau2.start=1,
                               a0=0, A0=.1, b0=0, B0=.1,
                               c0=-1, d0=-1, e0=0, E0=1,
                               store.ability=TRUE,
                               store.item=TRUE, ...
                               ){


  datamatrix <- as.matrix(datamatrix)

  nitems <- ncol(datamatrix)
  nsubj  <- nrow(datamatrix)
  ntime <- max(item.time.map)

  ## checks
  check.offset(list(...))
  check.mcmc.parameters(burnin, mcmc, thin)


  if (nitems != length(item.time.map)){
    cat("Number of rows of datamatrix not equal to length of item.time.map\n")
    stop("Please check data and try MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }
  if (min(item.time.map) != 1){
    cat("Minimum value in item.time.map not equal to 1\n")
    stop("Please check data and try MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }
  if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) !=
     (nitems * nsubj)) {
    cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
    stop("Please check data and try MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }


  if (A0 < 0){
    cat("Error: A0 (prior precision for alpha) is less than 0).\n")
    stop("Please respecify and try MCMCdynamicIRT1d() again.\n")
  }
  if (B0 < 0){
    cat("Error: B0 (prior precision for beta) is less than 0).\n")
    stop("Please respecify and try MCMCdynamicIRT1d() again.\n")
  }






  ## setup constraints on theta
  if(length(theta.constraints) != 0) {
    for (i in 1:length(theta.constraints)){
      theta.constraints[[i]] <-
        list(as.integer(1), theta.constraints[[i]][1])
    }
  }
  holder <- build.factor.constraints(theta.constraints, t(datamatrix), nsubj, 1)
  theta.eq.constraints <- holder[[1]]
  theta.ineq.constraints <- holder[[2]]
  ##subject.names <- holder[[3]]
  ## names
  item.names <- colnames(datamatrix)
  if (is.null(item.names)){
    item.names <- paste("item", 1:nitems, sep="")
  }






  ## starting values for theta error checking
  theta.start <- factor.score.start.check(theta.start, datamatrix,
                                          matrix(0,1,1), matrix(1,1,1),
                                          theta.eq.constraints,
                                          theta.ineq.constraints, 1)

  ## starting values for (alpha, beta)
  ab.starts <- matrix(NA, nitems, 2)
  for (i in 1:nitems){
    local.y <-  datamatrix[,i]
    ##local.y[local.y==9] <- NA
    if (length(na.omit(local.y)) <= 1){
      ab.starts[i,] <- c(0, 10)
    }
    else if (var(na.omit(local.y))==0){
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
    alpha.start <- alpha.start * matrix(1,nitems,1)
  }
  else if((dim(alpha.start)[1] != nitems) || (dim(alpha.start)[2] != 1)) {
    cat("Error: Starting value for alpha not conformable.\n")
    stop("Please respecify and call MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }
  if (is.na(beta.start)) {
    beta.start <- ab.starts[,2]
  }
  else if(is.null(dim(beta.start))) {
    beta.start <- beta.start * matrix(1,nitems,1)
  }
  else if((dim(beta.start)[1] != nitems) || (dim(beta.start)[2] != 1)) {
    cat("Error: Starting value for beta not conformable.\n")
    stop("Please respecify and call MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }




  ## generate time-specific theta information and create extended theta.start
  subject.names <- NULL
  theta.start.new <- NULL
  ## theta.info has:
  ## col1: subj ID, col2: #time periods, col3: offset (first term C indexing)
  theta.info <- matrix(NA, nsubj, 3)
  for (s in 1:nsubj){
    namestub <- rownames(datamatrix)[s]
    theta.info[s,1] <- s
    count <- 0
    holder <- NULL
    for (i in 1:nitems){
      if (!is.na(datamatrix[s,i])){
        holder <- c(holder, item.time.map[i])
      }
    }
    minholder <- min(holder)
    maxholder <- max(holder)
    theta.info[s,2] <- maxholder - minholder + 1
    theta.info[s,3] <- minholder - 1
    theta.start.new <- c(theta.start.new, rep(theta.start[s], theta.info[s,2]))
    subject.names <- c(subject.names,
                       paste(namestub, ".t", minholder:maxholder, sep=""))
  }
  nthetas <- length(subject.names)
  theta.start <- theta.start.new



  if (length(c0) < nsubj){
    c0 <- array(c0, nsubj)
  }
  if (length(d0) < nsubj){
    d0 <- array(d0, nsubj)
  }
  if (length(tau2.start) < nsubj){
    tau2.start <- array(tau2.start, nsubj)
  }
  if (length(e0) < nsubj){
    e0 <- array(e0, nsubj)
  }
  if (length(E0) < nsubj){
    E0 <- array(E0, nsubj)
  }
  E0inv <- 1/E0


  subject.names.short <- rownames(datamatrix)


  ## convert datamatrix into a sparse format where the missing values are
  ## not explicitly represented
  ##
  ## 1st column: subject index (C indexing)
  ## 2nd column: item index (C indexing)
  ## 3rd column: vote
  data.sparse.si <- NULL
  for (i in 1:nsubj){
    for (j in 1:nitems){
      if (!is.na(datamatrix[i,j])){
        data.sparse.si <- rbind(data.sparse.si, c(i-1, j-1, datamatrix[i,j]))
      }
    }
  }
  ## 1st column: item index (C indexing)
  ## 2nd column: subject index (C indexing)
  ## 3rd column: vote
##  data.sparse.is <- NULL
##  for (i in 1:nitems){
##    for (j in 1:nsubj){
##      if (!is.na(datamatrix[j,i])){
##        data.sparse.is <- rbind(data.sparse.is, c(i-1, j-1, datamatrix[j,i]))
##      }
##    }
##  }

  rm(datamatrix)






  ## define storage matrix for posterior theta draws
  thetadraws <- matrix(0, mcmc/thin, nthetas)
  if (store.ability != TRUE){
    thetadraws <- matrix(0, 1, 1)
  }
  alphadraws <- matrix(1, mcmc/thin, nitems)
  betadraws  <- matrix(2, mcmc/thin, nitems)
  if (store.item != TRUE){
    alphadraws <- matrix(1, 1, 1)
    betadraws  <- matrix(2, 1, 1)
  }
  tau2draws <- matrix(0, mcmc/thin, nsubj)


  ## seeds
  seeds <- form.seeds(seed)
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]

 ## print(theta.eq.constraints)
 ## print(theta.ineq.constraints)


#  return(list(theta=theta.start, alpha=alpha.start, beta=beta.start))


  ## call C++ code to draw sample
  posterior <- .C("cMCMCdynamicIRT1d",

                  thetadata = as.double(thetadraws),
                  thetarow = as.integer(nrow(thetadraws)),
                  thetacol = as.integer(ncol(thetadraws)),

                  alphadata = as.double(alphadraws),
                  alpharow = as.integer(nrow(alphadraws)),
                  alphacol = as.integer(ncol(alphadraws)),

                  betadata = as.double(betadraws),
                  betarow = as.integer(nrow(betadraws)),
                  betacol = as.integer(ncol(betadraws)),

                  tau2data = as.double(tau2draws),
                  tau2row = as.integer(nrow(tau2draws)),
                  tau2col = as.integer(ncol(tau2draws)),

                  nsubj = as.integer(nsubj),
                  nitems = as.integer(nitems),
                  ntime = as.integer(ntime),

                  Ydata = as.integer(data.sparse.si),
                  Yrow = as.integer(nrow(data.sparse.si)),
                  Ycol = as.integer(ncol(data.sparse.si)),

                  IT = as.integer(item.time.map-1),
                  ITlength = as.integer(length(item.time.map)),

                  burnin = as.integer(burnin),
                  mcmc = as.integer(mcmc),
                  thin = as.integer(thin),

                  lecuyer = as.integer(lecuyer),
                  seedarray = as.integer(seed.array),
                  lecuyerstream = as.integer(lecuyer.stream),

                  verbose = as.integer(verbose),

                  thetastartdata = as.double(theta.start),
                  thetastartlength = as.integer(length(theta.start)),

                  thetainfo = as.integer(theta.info),
                  thetainforow = as.integer(nrow(theta.info)),
                  thetainfocol = as.integer(ncol(theta.info)),

                  astartdata = as.double(alpha.start),
                  astartlength = as.integer(length(alpha.start)),

                  bstartdata = as.double(beta.start),
                  bstartlength = as.integer(length(beta.start)),

                  tau2data = as.double(tau2.start),
                  tau2length = as.integer(length(tau2.start)),

                  c0data = as.double(c0),
                  c0length = as.integer(length(c0)),

                  d0data = as.double(d0),
                  d0length = as.integer(length(d0)),

                  a0data = as.double(a0),
                  A0data = as.double(A0),

                  b0data = as.double(b0),
                  B0data = as.double(B0),

                  e0data = as.double(e0),
                  E0invdata = as.double(E0inv),

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





  tau2samp <- matrix(posterior$tau2data,
                      posterior$tau2row,
                      posterior$tau2col)
  colnames(tau2samp) <- paste("tau2.", subject.names.short, sep="")



  if (store.item & store.ability){
    thetasamp <- matrix(posterior$thetadata,
                        posterior$thetarow,
                        posterior$thetacol)
    colnames(thetasamp) <- paste("theta.", subject.names, sep="")

    alphasamp <- matrix(posterior$alphadata,
                        posterior$alpharow,
                        posterior$alphacol)
    colnames(alphasamp) <- paste("alpha.", item.names, sep="")


    betasamp <- matrix(posterior$betadata,
                       posterior$betarow,
                       posterior$betacol)
    colnames(betasamp) <- paste("beta.", item.names, sep="")

    outmat <- mcmc(cbind(thetasamp, alphasamp, betasamp, tau2samp),
                   start=1, end=mcmc, thin=thin)
  }
  else if (store.item){
    alphasamp <- matrix(posterior$alphadata,
                        posterior$alpharow,
                        posterior$alphacol)
    colnames(alphasamp) <- paste("alpha.", item.names, sep="")


    betasamp <- matrix(posterior$betadata,
                       posterior$betarow,
                       posterior$betacol)
    colnames(betasamp) <- paste("beta.", item.names, sep="")

    outmat <- mcmc(cbind(alphasamp, betasamp, tau2samp),
                   start=1, end=mcmc, thin=thin)
  }
  else if (store.ability){
    thetasamp <- matrix(posterior$thetadata,
                        posterior$thetarow,
                        posterior$thetacol)
    colnames(thetasamp) <- paste("theta.", subject.names, sep="")

    outmat <- mcmc(cbind(thetasamp, tau2samp),
                   start=1, end=mcmc, thin=thin)
  }
  else {
    outmat <- mcmc(tau2samp,
                   start=1, end=mcmc, thin=thin)
  }

  return(outmat)


}
