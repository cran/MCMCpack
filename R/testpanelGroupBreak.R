####################################################################
## test group-level breaks from panel residuals
##
## written by Jong Hee Park 03/2009
## modified and integrated with other codes by JHP 07/2011
######################################################################

#' A Test for the Group-level Break using a Multivariate Linear Regression
#' Model with Breaks
#'
#' testpanelGroupBreak fits a multivariate linear regression model with
#' parametric breaks using panel residuals to test the existence of group-level
#' breaks in panel residuals. The details are discussed in Park (2011).
#'
#'
#' \code{testpanelGroupBreak} fits a multivariate linear regression model with
#' parametric breaks using panel residuals to detect the existence of
#' system-level breaks in unobserved factors as discussed in Park (2011).
#'
#' The model takes the following form:
#'
#' \deqn{e_{i} \sim \mathcal{N}(\beta_{m}, \sigma^2_m I)\;\; m = 1, \ldots, M}
#'
#' We assume standard, semi-conjugate priors:
#'
#' \deqn{\beta \sim \mathcal{N}(b0, B0)}
#'
#' And:
#'
#' \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}
#'
#' Where \eqn{\beta} and \eqn{\sigma^{-2}} are
#' assumed \emph{a priori} independent.
#'
#' And:
#'
#' \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}
#'
#' Where \eqn{M} is the number of states.
#'
#' @param subject.id A numeric vector indicating the group number. It should
#' start from 1.
#'
#' @param time.id A numeric vector indicating the time unit. It should start
#' from 1.
#'
#' @param resid A vector of panel residuals
#'
#' @param m The number of changepoints.
#'
#' @param mcmc The number of MCMC iterations after burn-in.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0, the
#' iteration number and the posterior density samples are printed to the screen
#' every \code{verbose}th iteration.
#'
#' @param b0 The prior mean of the residual mean.
#'
#' @param B0 The prior precision of the residual variance
#'
#' @param c0 \eqn{c_0/2} is the shape parameter for the inverse Gamma
#' prior on \eqn{\sigma^2}. The amount of information in the inverse
#' Gamma prior is something like that from \eqn{c_0} pseudo-observations.
#'
#' @param d0 \eqn{d_0/2} is the scale parameter for the inverse Gamma
#' prior on \eqn{\sigma^2}.
#'
#' @param a \eqn{a} is the shape1 beta prior for transition probabilities.
#' By default, the expected duration is computed and corresponding a and b
#' values are assigned. The expected duration is the sample period divided by
#' the number of states.
#'
#' @param b \eqn{b} is the shape2 beta prior for transition probabilities.
#' By default, the expected duration is computed and corresponding a and b
#' values are assigned. The expected duration is the sample period divided by
#' the number of states.
#'
#' @param seed The seed for the random number generator.  If NA, current R
#' system seed is used.
#'
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated and \code{Chib95} in which case the method of Chib (1995) is
#' used.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample. This object can
#' be summarized by functions provided by the coda package.  The object
#' contains an attribute \code{prob.state} storage matrix that contains the
#' probability of \eqn{state_i} for each period, and the log-marginal
#' likelihood of the model (\code{logmarglike}).
#'
#' @export
#'
#' @references Jong Hee Park, 2012. ``Unified Method for Dynamic and
#' Cross-Sectional Heterogeneity: Introducing Hidden Markov Panel Models.''
#' \emph{American Journal of Political Science}.56: 1040-1054.
#' <doi: 10.1111/j.1540-5907.2012.00590.x>
#'
#' Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point
#' models.'' \emph{Journal of Econometrics}. 86: 221-241.
#' <doi: 10.1080/01621459.1995.10476635>
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#'    ## data generating
#'   set.seed(1977)
#'   Q <- 3
#'   true.beta1   <-  c(1, 1, 1) ; true.beta2   <-  c(1, -1, -1)
#'   true.sigma2 <-  c(1, 3); true.D1 <- diag(.5, Q); true.D2 <- diag(2.5, Q)
#'   N=20; T=100;
#'   NT <- N*T
#'   x1 <- rnorm(NT)
#'   x2 <- runif(NT, 5, 10)
#'   X <- cbind(1, x1, x2);   W <- X;   y <- rep(NA, NT)
#'
#'   ## true break numbers are one and at the center
#'   break.point = rep(T/2, N); break.sigma=c(rep(1, N));
#'   break.list <- rep(1, N)
#'   id  <-  rep(1:N, each=NT/N)
#'   K <-  ncol(X);
#'   ruler <- c(1:T)
#'
#'   ## compute the weight for the break
#'   W.mat <- matrix(NA, T, N)
#'   for (i in 1:N){
#'     W.mat[, i] <- pnorm((ruler-break.point[i])/break.sigma[i])
#'   }
#'   Weight <- as.vector(W.mat)
#'
#'   ## data generating by weighting two means and variances
#'   j = 1
#'   for (i in 1:N){
#'     Xi <- X[j:(j+T-1), ]
#'     Wi <- W[j:(j+T-1), ]
#'     true.V1 <- true.sigma2[1]*diag(T) + Wi%*%true.D1%*%t(Wi)
#'     true.V2 <- true.sigma2[2]*diag(T) + Wi%*%true.D2%*%t(Wi)
#'     true.mean1 <- Xi%*%true.beta1
#'     true.mean2 <- Xi%*%true.beta2
#'     weight <- Weight[j:(j+T-1)]
#'     y[j:(j+T-1)] <- (1-weight)*true.mean1 + (1-weight)*chol(true.V1)%*%rnorm(T) +
#'       weight*true.mean2 + weight*chol(true.V2)%*%rnorm(T)
#'     j <- j + T
#'   }
#'   ## model fitting
#'   subject.id <- c(rep(1:N, each=T))
#'   time.id <- c(rep(1:T, N))
#'
#'
#'   resid <- rstandard(lm(y ~X-1 + as.factor(subject.id)))
#'   G <- 100
#'   out0 <- testpanelGroupBreak(subject.id, time.id, resid, m=0,
#'            mcmc=G, burnin=G, thin=1, verbose=G,
#'            b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
#'   out1 <- testpanelGroupBreak(subject.id, time.id, resid, m=1,
#'            mcmc=G, burnin=G, thin=1, verbose=G,
#'            b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
#'   out2 <- testpanelGroupBreak(subject.id, time.id, resid, m=2,
#'            mcmc=G, burnin=G, thin=1, verbose=G,
#'            b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
#'   out3 <- testpanelGroupBreak(subject.id, time.id, resid, m=3,
#'            mcmc=G, burnin=G, thin=1, verbose=G,
#'            b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
#'
#'   ## Note that the code is for a hypothesis test of no break in panel residuals.
#'   ## When breaks exist, the estimated number of break in the mean and variance of panel residuals
#'   ## tends to be larger than the number of break in the data generating process.
#'   ## This is due to the difference in parameter space, not an error of the code.
#'   BayesFactor(out0, out1, out2, out3)
#'
#'   ## In order to identify the number of breaks in panel parameters,
#'   ## use HMMpanelRE() instead.
#'
#' }
#'
"testpanelGroupBreak" <-
  function(subject.id, time.id, resid, m=1,
           mcmc=1000, burnin=1000, thin=1, verbose=0,
           b0, B0, c0, d0, a = NULL, b = NULL,
           seed = NA, marginal.likelihood = c("none", "Chib95"), ...){
    ## beta.start and sigma2.start are not arguments. OLS estimates will be used!
    ## subject.id is a numeric list indicating the group number. It should start from 1.
    ## time.id is a numeric list indicating the time, starting from 1.
    ## subject.offset is the obs number from which a new subject unit starts
    ## time.offset is the obs number from which a new time unit starts when we stack data by time.id

    cl <- match.call()
    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    ns <- m + 1
    nobs <- length(subject.id)
    newY <- matrix(resid, nobs, 1)
    newX <- matrix(1, nobs, 1)
    K <-  1

    ## Sort Data based on time.id
    oldTSCS <- cbind(time.id, subject.id, newY, newX)
    newTSCS <- oldTSCS[order(oldTSCS[,1]),]
    newYT <- as.matrix(newTSCS[,3])
    newXT <- as.matrix(newTSCS[,4])
    b0 <- as.matrix(b0)
    B0 <- as.matrix(B0)

    nstore <- mcmc/thin
    nsubj <- length(unique(subject.id))
    ## subject.groupinfo matrix
    if (unique(subject.id)[1] != 1){
      stop("subject.id should start 1!")
    }
    subject.offset <- c(0, which(diff(sort(subject.id))==1)[-nsubj])
    nsubject.vec <- rep(NA, nsubj)
    for (i in 1:nsubj){
      nsubject.vec[i] <- sum(subject.id==unique(subject.id)[i])
    }
    subject.groupinfo <- cbind(unique(subject.id), subject.offset, nsubject.vec)

    ## time.groupinfo
    if(unique(time.id)[1] != 1){
      time.id <- time.id - unique(time.id)[1] + 1
      cat("time.id does not start from 1. So it is modified by subtracting the first unit of time.")
    }
    ntime <- max(nsubject.vec)## maximum time length
    ntime.vec <- rep(NA, ntime)
    for (i in 1:ntime){
      ntime.vec[i] <- sum(time.id==unique(time.id)[i])
    }
    time.offset <- c(0, which(diff(sort(time.id))==1)[-ntime])
    time.groupinfo <- cbind(unique(time.id), time.offset, ntime.vec)

    ## prior inputs
    if (m > 0){
      P0 <- trans.mat.prior(m=m, n=ntime, a=a, b=b)
    }
    else {
      P0 <- matrix(1, 1, 1)
    }
    betadraws <- matrix(data=0, nstore, ns*K)
    sigmadraws <- matrix(data=0, nstore, ns)
    psdraws <- matrix(data=0, ntime, ns)
    ols <- lm(newY ~ newX - 1)
    beta.start  <- rep(coef(ols)[1], ns)
    sigma2.start <- summary(ols)$sigma^2

    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)

    ## following MCMCregress, set chib as binary
    logmarglike <- loglik <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    ## call C++ code to draw sample
    posterior <- .C("HMMmultivariateGaussian",
                    betadata = as.double(betadraws),
                    betarow = as.integer(nrow(betadraws)),
                    betacol = as.integer(ncol(betadraws)),

                    sigmadata = as.double(sigmadraws),
                    psout = as.double(psdraws),
                    nsubj = as.integer(nsubj), ntime = as.integer(ntime),
                    m = as.integer(m),
                    nobs = as.integer(nobs),
                    subjectid = as.integer(subject.id),
                    timeid = as.integer(time.id),

                    Ydata = as.double(newY), Yrow = as.integer(nrow(newY)), Ycol = as.integer(ncol(newY)),
                    Xdata = as.double(newX), Xrow = as.integer(nrow(newX)), Xcol = as.integer(ncol(newX)),
                    YTdata = as.double(newYT), XTdata = as.double(newXT),
                    burnin = as.integer(burnin), mcmc = as.integer(mcmc),
                    thin = as.integer(thin), verbose = as.integer(verbose),

                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),

                    betastartdata = as.double(beta.start),
                    sigma2start = as.double(sigma2.start),

                    b0data = as.double(b0), B0data = as.double(B0),
                    c0 = as.double(c0), d0 = as.double(d0),
                    P0data = as.double(P0),
                    P0row = as.integer(nrow(P0)),
                    P0col = as.integer(ncol(P0)),

                    subject_groupinfodata = as.double(subject.groupinfo),
                    time_groupinfodata = as.double(time.groupinfo),

                    logmarglikeholder = as.double(0),
                    loglikeholder = as.double(0),
                    chib = as.integer(chib),
                    PACKAGE="MCMCpack"
                    )

    ## pull together matrix and build MCMC object to return
    beta.samp <- matrix(posterior$betadata,
                        posterior$betarow,
                        posterior$betacol)
    ## stored by the order of (11, 12, 13, 21, 22, 23)
    sigma.samp <- matrix(posterior$sigmadata,
                         posterior$betarow,
                         ns)
    xnames <-  sapply(c(1:K), function(i){paste("beta", i, sep = "")})
    output1 <- mcmc(data=beta.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
    output2 <- mcmc(data=sigma.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
    if(m>1){
      varnames(output1)  <- sapply(c(1:ns),
                                   function(i){
                                     paste(xnames, "_regime", i, sep = "")
                                   })
      varnames(output2)  <- sapply(c(1:ns),
                                   function(i){
                                     paste("sigma2_regime", i, sep = "")
                                   })
    }
    output <- as.mcmc(cbind(output1, output2))

    attr(output, "title") <- "testpanelGroupBreak Posterior Sample"
    attr(output, "call")   <- cl
    attr(output, "y")       <- resid[1:ntime]
    attr(output, "m")       <- m
    attr(output, "nsubj")   <- nsubj
    attr(output, "ntime")   <- ntime
    if(m>0){
      ps.holder   <- matrix(posterior$psout, ntime, ns)
      attr(output, "prob.state") <- ps.holder/nstore
    }
    attr(output, "logmarglike") <- posterior$logmarglikeholder
    attr(output, "loglike") <- posterior$loglikeholder

    ## report the results in a simple manner


    return(output)
  }
