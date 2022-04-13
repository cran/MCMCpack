##########################################################################
## MCMCmnl.R samples from the posterior distribution of a multinomial
## logit model using Metropolis-Hastings.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 12/22/2004
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

## parse formula and return a list that contains the model response
## matrix as element one, the model matrix as element two,
## the column names of X as element three, the rownames of
## X and y as element four, and the number of choices in the largest
## choice set in element five
"parse.formula.mnl" <- function(formula, data, baseline=NULL,
                                intercept=TRUE){

  ## extract Y, X, and variable names for model formula and frame
  mt <- terms(formula, data=data)
  if(missing(data)) data <- sys.frame(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- mf$baseline <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")

  mf <- eval(mf, sys.frame(sys.parent()))
  if (!intercept){
    attributes(mt)$intercept <- 0
  }


  ## deal with Y
  Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
  if (ncol(Y)==1){
    Y <- factor(Y)
    number.choices <- length(unique(Y))
    choice.names <- sort(unique(Y))
    Ymat <- matrix(NA, length(Y), number.choices)
    colnames(Ymat) <- choice.names
    for (i in 1:(number.choices)){
      Ymat[,i] <- as.numeric(Y==choice.names[i])
    }
  }
  else{
    ## this block will allow for nonconstant choice sets
    number.choices <- ncol(Y)
    Ytemp <- Y
    Ytemp[Y== -999] <- NA
    if ( min(unique(array(Y)) %in% c(-999,0,1))==0 ||
        min(apply(Ytemp, 1, sum, na.rm=TRUE) == rep(1, nrow(Y)))==0){
      stop("Y is a matrix but it is not composed of 0/1/-999 values\n and/or rows do not sum to 1\n")
    }
    Ymat <- Y
    choice.names <- colnames(Y)
  }
  colnames(Ymat) <- choice.names
  rownames(Ymat) <- 1:nrow(Ymat)

  #Y.long <- matrix(t(Ymat), length(Ymat), 1)
  #colnames(Y.long) <- "Y"
  #rownames(Y.long) <- rep(1:nrow(Ymat), rep(length(choice.names), nrow(Ymat)))
  #rownames(Y.long) <- paste(rownames(Y.long), choice.names, sep=".")
  #group.id <- rep(1:nrow(Ymat), rep(ncol(Ymat), nrow(Ymat)))

  ## deal with X
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)
  X <- as.matrix(X)         # X matrix
  xvars <- dimnames(X)[[2]] # X variable names
  xobs  <- dimnames(X)[[1]] # X observation names

  if (is.null(baseline)){
    baseline <- choice.names[1]
  }
  if (! baseline %in% choice.names){
    stop("'baseline' not consistent with observed choice levels in y\n")
  }

  ## deal with choice specific covariates
  choicevar.indic <- rep(FALSE, length(xvars)) # indicators for choice
                                               # specific variables
  choicevar.indic[grep("^choicevar\\(", xvars)] <- TRUE
  if (sum(choicevar.indic) > 0){
    cvarname1.vec <- rep(NA, sum(choicevar.indic))
    cvarname2.vec <- rep(NA, sum(choicevar.indic))
    counter <- 0
    for (i in 1:length(xvars)){
      if (choicevar.indic[i]){
        counter <- counter + 1
        vn2 <- strsplit(xvars[i], ",")
        vn3 <- strsplit(vn2[[1]], "\\(")
        vn4 <- strsplit(vn3[[3]], "=")
        cvarname1 <- vn3[[2]][1]
        cvarname1 <- strsplit(cvarname1, "\"")[[1]]
        cvarname1 <- cvarname1[length(cvarname1)]
        cvarname2 <- vn4[[1]][length(vn4[[1]])]
        cvarname2 <- strsplit(cvarname2, "\"")[[1]][2]
        if (! cvarname2 %in% choice.names){
          stop("choicelevel that was set in choicevar() not consistent with\n observed choice levels in y")
        }
        cvarname1.vec[counter] <- cvarname1
        cvarname2.vec[counter] <- cvarname2
        xvars[i] <- paste(cvarname1, cvarname2, sep=".")
      }
    }

    X.cho <- X[, choicevar.indic]
    X.cho.mat <- matrix(NA, length(choice.names)*nrow(X),
                        length(unique(cvarname1.vec)))
    rownames(X.cho.mat) <- rep(rownames(X), rep(length(choice.names), nrow(X)))
    rownames(X.cho.mat) <- paste(rownames(X.cho.mat), choice.names, sep=".")
    colnames(X.cho.mat) <- unique(cvarname1.vec)
    choice.names.n <- rep(choice.names, nrow(X))
    for (j in 1:length(unique(cvarname1.vec))){
      for (i in 1:length(cvarname2.vec)){
        if (colnames(X.cho.mat)[j] == cvarname1.vec[i]){
          X.cho.mat[choice.names.n==cvarname2.vec[i], j] <-
            X.cho[,i]
        }
      }
    }
  }


  ## deal with individual specific covariates
  xvars.ind.mat <- rep(xvars[!choicevar.indic],
                       rep(length(choice.names),
                           sum(!choicevar.indic)))
  xvars.ind.mat <- paste(xvars.ind.mat, choice.names, sep=".")

  if (sum(!choicevar.indic) > 0){
    X.ind.mat <- X[,!choicevar.indic] %x% diag(length(choice.names))
    colnames(X.ind.mat) <- xvars.ind.mat
    rownames(X.ind.mat) <- rep(rownames(X), rep(length(choice.names), nrow(X)))
    rownames(X.ind.mat) <- paste(rownames(X.ind.mat), choice.names, sep=".")

    ## delete columns correpsonding to the baseline choice
    ivarname1 <- strsplit(xvars.ind.mat, "\\.")
    ivar.keep.indic <- rep(NA, ncol(X.ind.mat))
    for (i in 1:ncol(X.ind.mat)){
      ivar.keep.indic[i] <- ivarname1[[i]][length(ivarname1[[i]])] != baseline
    }
    X.ind.mat <- X.ind.mat[,ivar.keep.indic]
  }

  if (sum(choicevar.indic) > 0 & sum(!choicevar.indic) > 0){
    X <- cbind(X.cho.mat, X.ind.mat)
  }
  else if (sum(!choicevar.indic) > 0){
    X <- X.ind.mat
  }
  else if (sum(choicevar.indic) > 0){
    X <- X.cho.mat
  }
  else {
    stop("X matrix appears to have neither choice-specific nor individual-specific variables.\n")
  }
  #Y <- Y.long
  xvars <- colnames(X)
  xobs <- rownames(X)

  return(list(Ymat, X, xvars, xobs, number.choices))

}

## dummy function used to handle choice-specific covariates
#' Handle Choice-Specific Covariates in Multinomial Choice Models
#'
#' This function handles choice-specific covariates in multinomial choice
#' models.  See the example for an example of useage.
#'
#' @param var The is the name of the variable in the dataframe.
#'
#' @param varname The name of the new variable to be created.
#'
#' @param choicelevel The level of \code{y} that the variable corresponds to.
#'
#' @export
#'
#' @return The new variable used by the \code{MCMCmnl()} function.
#'
#' @seealso \code{\link{MCMCmnl}}
#'
#' @keywords manip
"choicevar" <- function(var, varname, choicelevel){
  junk1 <- varname
  junk2 <- choicelevel
  return(var)
}


## MNL log-posterior function (used to get starting values)
## vector Y without NAs
"mnl.logpost.noNA" <- function(beta, new.Y, X, b0, B0){
  nobs <- length(new.Y)
  ncat <- nrow(X) / nobs
  Xb <- X %*% beta
  Xb <- matrix(Xb, byrow=TRUE, ncol=ncat)
  indices <- cbind(1:nobs, new.Y)
  Xb.reform <- Xb[indices]
  eXb <- exp(Xb)
  #denom <- log(apply(eXb, 1, sum))
  z <- rep(1, ncat)
  denom <- log(eXb %*% z)

  log.prior <- 0.5 * t(beta - b0) %*% B0 %*% (beta - b0)

  return(sum(Xb.reform - denom) + log.prior)
}

## MNL log-posterior function (used to get starting values)
## matrix Y with NAs
"mnl.logpost.NA" <- function(beta, Y, X, b0, B0){
    k <- ncol(X)
    numer <- exp(X %*% beta)
    numer[Y== -999] <- NA
    numer.mat <- matrix(numer, nrow(Y), ncol(Y), byrow=TRUE)
    denom <- apply(numer.mat, 1, sum, na.rm=TRUE)
    choice.probs <- numer.mat / denom
    Yna <- Y
    Yna[Y == -999] <- NA
    log.like.mat <- log(choice.probs) * Yna
    log.like <- sum(apply(log.like.mat, 1, sum, na.rm=TRUE))

    log.prior <- 0.5 * t(beta - b0) %*% B0 %*% (beta - b0)

    return(log.like + log.prior)
}



#' Markov Chain Monte Carlo for Multinomial Logistic Regression
#'
#' This function generates a sample from the posterior distribution of a
#' multinomial logistic regression model using either a random walk Metropolis
#' algorithm or a slice sampler. The user supplies data and priors, and a
#' sample from the posterior distribution is returned as an mcmc object, which
#' can be subsequently analyzed with functions provided in the coda package.
#'
#' \code{MCMCmnl} simulates from the posterior distribution of a multinomial
#' logistic regression model using either a random walk Metropolis algorithm or
#' a univariate slice sampler. The simulation proper is done in compiled C++
#' code to maximize efficiency.  Please consult the coda documentation for a
#' comprehensive list of functions that can be used to analyze the posterior
#' sample.
#'
#' The model takes the following form:
#'
#' \deqn{y_i \sim \mathcal{M}ultinomial(\pi_i)}
#'
#' where:
#'
#' \deqn{\pi_{ij} = \frac{\exp(x_{ij}'\beta)}{\sum_{k=1}^p\exp(x_{ik}'\beta)}}
#'
#' We assume a multivariate Normal prior on \eqn{\beta}:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}
#'
#' The Metropolis proposal distribution is centered at the current
#' value of \eqn{\beta} and has variance-covariance
#' \eqn{V = T(B_0 + C^{-1})^{-1} T}, where \eqn{T} is a the
#' diagonal positive definite matrix formed from the \code{tune},
#' \eqn{B_0} is the prior precision, and \eqn{C} is the large sample
#' variance-covariance matrix of the MLEs. This last calculation is
#' done via an initial call to \code{optim}.
#'
#' @param formula Model formula.
#'
#' If the choicesets do not vary across individuals, the \code{y} variable
#' should be a factor or numeric variable that gives the observed choice of
#' each individual. If the choicesets do vary across individuals, \code{y}
#' should be a \eqn{n \times p} matrix where \eqn{n} is the number of
#' individuals and \eqn{p} is the maximum number of choices in any
#' choiceset.  Here each column of \code{y} corresponds to a particular
#' observed choice and the elements of \code{y} should be either \code{0} (not
#' chosen but available), \code{1} (chosen), or \code{-999} (not available).
#'
#' Choice-specific covariates have to be entered using the syntax:
#' \code{choicevar(cvar, "var", "choice")} where \code{cvar} is the name of a
#' variable in \code{data}, \code{"var"} is the name of the new variable to be
#' created, and \code{"choice"} is the level of \code{y} that \code{cvar}
#' corresponds to. Specifying each choice-specific covariate will typically
#' require \eqn{p} calls to the \code{choicevar} function in the formula.
#'
#' Individual specific covariates can be entered into the formula normally.
#'
#' See the examples section below to see the syntax used to fit various models.
#' @param baseline The baseline category of the response variable.
#'
#' \code{baseline} should be set equal to one of the observed levels of the
#' response variable. If left equal to \code{NULL} then the baseline level is
#' set to the alpha-numerically first element of the response variable. If the
#' choicesets vary across individuals, the baseline choice must be in the
#' choiceset of each individual.
#'
#' @param data The data frame used for the analysis. Each row of the dataframe
#' should correspond to an individual who is making a choice.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of iterations to run the sampler past burn-in.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' mcmc iterations must be divisible by this value.
#'
#' @param mcmc.method Can be set to either "IndMH" (default), "RWM", or "slice"
#' to perform independent Metropolis-Hastings sampling, random walk Metropolis
#' sampling or slice sampling respectively.
#'
#' @param tdf Degrees of freedom for the multivariate-t proposal distribution
#' when \code{mcmc.method} is set to "IndMH". Must be positive.
#'
#' @param tune Metropolis tuning parameter. Can be either a positive scalar or
#' a \eqn{k}-vector, where \eqn{k} is the length of \eqn{\beta}.
#' Make sure that the acceptance rate is satisfactory (typically between 0.20
#' and 0.5) before using the posterior sample for inference.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the current beta vector, and the Metropolis acceptance
#' rate are printed to the screen every \code{verbose}th iteration.
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
#' @param beta.start The starting value for the \eqn{\beta} vector.  This
#' can either be a scalar or a column vector with dimension equal to the number
#' of betas. If this takes a scalar value, then that value will serve as the
#' starting value for all of the betas.  The default value of NA will use the
#' maximum likelihood estimate of \eqn{\beta} as the starting value.
#'
#' @param b0 The prior mean of \eqn{\beta}.  This can either be a scalar
#' or a column vector with dimension equal to the number of betas. If this
#' takes a scalar value, then that value will serve as the prior mean for all
#' of the betas.
#'
#' @param B0 The prior precision of \eqn{\beta}. This can either be a
#' scalar or a square matrix with dimensions equal to the number of betas.  If
#' this takes a scalar value, then that value times an identity matrix serves
#' as the prior precision of \eqn{\beta}. Default value of 0 is
#' equivalent to an improper uniform prior for beta.
#'
#' @param ... Further arguments to be passed.
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[nnet]{multinom}}
#'
#' @references
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Radford Neal. 2003. ``Slice Sampling'' (with discussion). \emph{Annals of
#' Statistics}, 31: 705-767.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' Siddhartha Chib, Edward Greenberg, and Yuxin Chen. 1998.  ``MCMC Methods for
#' Fitting and Comparing Multinomial Response Models."
#'
#' @keywords models
#'
#' @examples
#'
#'   \dontrun{
#'   data(Nethvote)
#'
#'   ## just a choice-specific X var
#'   post1 <- MCMCmnl(vote ~
#'                 choicevar(distD66, "sqdist", "D66") +
#'                 choicevar(distPvdA, "sqdist", "PvdA") +
#'                 choicevar(distVVD, "sqdist", "VVD") +
#'                 choicevar(distCDA, "sqdist", "CDA"),
#'                 baseline="D66", mcmc.method="IndMH", B0=0,
#'                 verbose=500, mcmc=100000, thin=10, tune=1.0,
#'                 data=Nethvote)
#'
#'   plot(post1)
#'   summary(post1)
#'
#'
#'
#'   ## just individual-specific X vars
#'   post2<- MCMCmnl(vote ~
#'                 relig + class + income + educ + age + urban,
#'                 baseline="D66", mcmc.method="IndMH", B0=0,
#'                 verbose=500, mcmc=100000, thin=10, tune=0.5,
#'                 data=Nethvote)
#'
#'   plot(post2)
#'   summary(post2)
#'
#'
#'
#'   ## both choice-specific and individual-specific X vars
#'   post3 <- MCMCmnl(vote ~
#'                 choicevar(distD66, "sqdist", "D66") +
#'                 choicevar(distPvdA, "sqdist", "PvdA") +
#'                 choicevar(distVVD, "sqdist", "VVD") +
#'                 choicevar(distCDA, "sqdist", "CDA") +
#'                 relig + class + income + educ + age + urban,
#'                 baseline="D66", mcmc.method="IndMH", B0=0,
#'                 verbose=500, mcmc=100000, thin=10, tune=0.5,
#'                 data=Nethvote)
#'
#'   plot(post3)
#'   summary(post3)
#'
#'
#'   ## numeric y variable
#'   nethvote$vote <- as.numeric(nethvote$vote)
#'   post4 <- MCMCmnl(vote ~
#'                 choicevar(distD66, "sqdist", "2") +
#'                 choicevar(distPvdA, "sqdist", "3") +
#'                 choicevar(distVVD, "sqdist", "4") +
#'                 choicevar(distCDA, "sqdist", "1") +
#'                 relig + class + income + educ + age + urban,
#'                 baseline="2", mcmc.method="IndMH", B0=0,
#'                 verbose=500, mcmc=100000, thin=10, tune=0.5,
#'                 data=Nethvote)
#'
#'
#'   plot(post4)
#'   summary(post4)
#'
#'
#'
#'   ## Simulated data example with nonconstant choiceset
#'   n <- 1000
#'   y <- matrix(0, n, 4)
#'   colnames(y) <- c("a", "b", "c", "d")
#'   xa <- rnorm(n)
#'   xb <- rnorm(n)
#'   xc <- rnorm(n)
#'   xd <- rnorm(n)
#'   xchoice <- cbind(xa, xb, xc, xd)
#'   z <- rnorm(n)
#'   for (i in 1:n){
#'     ## randomly determine choiceset (c is always in choiceset)
#'     choiceset <- c(3, sample(c(1,2,4), 2, replace=FALSE))
#'     numer <- matrix(0, 4, 1)
#'     for (j in choiceset){
#'       if (j == 3){
#'         numer[j] <- exp(xchoice[i, j] )
#'       }
#'       else {
#'         numer[j] <- exp(xchoice[i, j] - z[i] )
#'       }
#'     }
#'     p <- numer / sum(numer)
#'     y[i,] <- rmultinom(1, 1, p)
#'     y[i,-choiceset] <- -999
#'   }
#'
#'   post5 <- MCMCmnl(y~choicevar(xa, "x", "a") +
#'                   choicevar(xb, "x", "b") +
#'                   choicevar(xc, "x", "c") +
#'                   choicevar(xd, "x", "d") + z,
#'                   baseline="c", verbose=500,
#'                   mcmc=100000, thin=10, tune=.85)
#'
#'   plot(post5)
#'   summary(post5)
#'
#'   }
#'
"MCMCmnl" <-
  function(formula, baseline=NULL, data=NULL,
           burnin = 1000, mcmc = 10000, thin=1,
           mcmc.method = "IndMH", 
           tune = 1.0, tdf=6, verbose = 0, seed = NA,
           beta.start = NA, b0 = 0, B0 = 0, ...) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    if (tdf <= 0){
      stop("degrees of freedom for multivariate-t proposal must be positive.\n Respecify tdf and try again.\n")
    }

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrix
    holder <- parse.formula.mnl(formula=formula, baseline=baseline,
                                data=data)
    Y <- holder[[1]]
    ## check to make sure baseline category is always available in choiceset
    if (is.null(baseline)){
      if (max(Y[,1] == -999) == 1){
        stop("Baseline choice not available in all choicesets.\n Respecify baseline category and try again.\n")
      }
    }
    else{
      if (max(Y[,baseline] == -999) == 1){
        stop("Baseline choice not available in all choicesets.\n Respecify baseline category and try again.\n")
      }
    }
    X <- holder[[2]]
    xnames <- holder[[3]]
    xobs <- holder[[4]]
    number.choices <- holder[[5]]
    K <- ncol(X)  # number of covariates

    ## form the tuning parameter
    tune <- vector.tune(tune, K)

    ## priors and starting values
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    beta.init <- rep(0, K)
    cat("Calculating MLEs and large sample var-cov matrix.\n")
    cat("This may take a moment...\n")
    if (max(is.na(Y))){
      optim.out <- optim(beta.init, mnl.logpost.NA, method="BFGS",
                         control=list(fnscale=-1),
                         hessian=TRUE, Y=Y, X=X, b0=b0, B0=B0)
    }
    else{
      new.Y <- apply(Y==1, 1, which)
      optim.out <- optim(beta.init, mnl.logpost.noNA, method="BFGS",
                         control=list(fnscale=-1),
                         hessian=TRUE, new.Y=new.Y, X=X, b0=b0, B0=B0)
    }
    cat("Inverting Hessian to get large sample var-cov matrix.\n")
    ##V <- solve(-1*optim.out$hessian)
    V <- chol2inv(chol(-1*optim.out$hessian))
    beta.mode <- matrix(optim.out$par, K, 1)


    if (is.na(beta.start) || is.null(beta.start)){
      beta.start <- matrix(optim.out$par, K, 1)
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- matrix(beta.start, K, 1)
    }
    else if (length(beta.start != K)){
      stop("beta.start not of appropriate dimension\n")
    }

    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, dim(X)[2] )
    posterior <- NULL

    if (mcmc.method=="RWM"){
      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCmnlMH",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       tune=tune, lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose),
                       betastart=beta.start, betamode=beta.mode,
                       b0=b0, B0=B0,
                       V=V, RW=as.integer(1), tdf=as.double(tdf))

      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCmnl Posterior Sample")
    }
    else if (mcmc.method=="IndMH"){
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCmnlMH",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       tune=tune, lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose),
                       betastart=beta.start, betamode=beta.mode,
                       b0=b0, B0=B0,
                       V=V, RW=as.integer(0), tdf=as.double(tdf))

      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCmnl Posterior Sample")

    }
    else if (mcmc.method=="slice"){
      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCmnlslice",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose), betastart=beta.start,
                       b0=b0, B0=B0, V=V)

      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCmnl Posterior Sample")

    }
    else{
        cat("\n\nmcmc.method not equal to one of 'RWM', 'IndMH', or 'slice'.\n")
        stop("Please respecifify and call MCMCmnl() again.\n")
    }

    return(output)

  }
