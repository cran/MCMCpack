####################################################################
## HMM Gaussian Panel with Time Varying Intercepts
##
## first written by Jong Hee Park on 02/2009
## revised for MCMCpack inclusion on 09/2011
######################################################################
#' Markov Chain Monte Carlo for the Hidden Markov Fixed-effects Model
#'
#' HMMpanelFE generates a sample from the posterior distribution of
#' the fixed-effects model with varying individual effects model
#' discussed in Park (2011).  The code works for both balanced and
#' unbalanced panel data as long as there is no missing data in the
#' middle of each group.  This model uses a multivariate Normal prior
#' for the fixed effects parameters and varying individual effects, an
#' Inverse-Gamma prior on the residual error variance, and Beta prior
#' for transition probabilities. The user supplies data and priors,
#' and a sample from the posterior distribution is returned as an mcmc
#' object, which can be subsequently analyzed with functions provided
#' in the coda package.
#'
#' \code{HMMpanelFE} simulates from the fixed-effect hidden Markov
#' pbject level: \deqn{\varepsilon_{it} \sim \mathcal{N}(\alpha_{im},
#' \sigma^2_{im})}
#'
#' We assume standard, semi-conjugate priors: \deqn{\beta \sim
#' \mathcal{N}(b_0,B_0^{-1})} And: \deqn{\sigma^{-2} \sim
#' \mathcal{G}amma(c_0/2, d_0/2)} And: \deqn{\alpha \sim
#' \mathcal{N}(delta_0,Delta_0^{-1})} \eqn{\beta}, \eqn{\alpha} and
#' \eqn{\sigma^{-2}} are assumed \emph{a priori} independent.
#'
#' And: \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}
#' Where \eqn{M} is the number of states.
#'
#' OLS estimates are used for starting values.
#'
#' @param subject.id A numeric vector indicating the group number. It
#'   should start from 1.
#'
#' @param y The response variable.
#'
#' @param X The model matrix excluding the constant.
#'
#' @param m A vector of break numbers for each subject in the panel.
#'
#' @param mcmc The number of MCMC iterations after burn-in.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The
#'   number of MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the
#'   progress of the sampler is printed to the screen.  If
#'   \code{verbose} is greater than 0, the iteration number and the
#'   posterior density samples are printed to the screen every
#'   \code{verbose}th iteration.
#'
#' @param b0 The prior mean of \eqn{\beta}.  This can either be a
#'   scalar or a column vector with dimension equal to the number of
#'   betas. If this takes a scalar value, then that value will serve
#'   as the prior mean for all of the betas.
#'
#' @param B0 The prior precision of \eqn{\beta}.  This can either be a
#'   scalar or a square matrix with dimensions equal to the number of
#'   betas.  If this takes a scalar value, then that value times an
#'   identity matrix serves as the prior precision of beta. Default
#'   value of 0 is equivalent to an improper uniform prior for beta.
#'
#' @param c0 \eqn{c_0/2} is the shape parameter for the inverse Gamma
#'   prior on \eqn{\sigma^2} (the variance of the disturbances). The
#'   amount of information in the inverse Gamma prior is something
#'   like that from \eqn{c_0} pseudo-observations.
#'
#' @param d0 \eqn{d_0/2} is the scale parameter for the inverse Gamma
#'   prior on \eqn{\sigma^2} (the variance of the disturbances). In
#'   constructing the inverse Gamma prior, \eqn{d_0} acts like the sum
#'   of squared errors from the \eqn{c_0} pseudo-observations.
#'
#' @param delta0 The prior mean of \eqn{\alpha}.
#'
#' @param Delta0 The prior precision of \eqn{\alpha}.
#'
#' @param a \eqn{a} is the shape1 beta prior for transition
#'   probabilities.  By default, the expected duration is computed and
#'   corresponding a and b values are assigned. The expected duration
#'   is the sample period divided by the number of states.
#'
#' @param b \eqn{b} is the shape2 beta prior for transition
#'   probabilities.  By default, the expected duration is computed and
#'   corresponding a and b values are assigned. The expected duration
#'   is the sample period divided by the number of states.
#'
#' @param seed The seed for the random number generator.  If NA,
#'   current R system seed is used.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample. This
#'   object can be summarized by functions provided by the coda
#'   package.  The object contains an attribute \code{sigma} storage
#'   matrix that contains time-varying residual variance, an attribute
#'   \code{state} storage matrix that contains posterior samples of
#'   hidden states, and an attribute \code{delta} storage matrix
#'   containing time-varying intercepts.
#'
#' @export
#'
#' @references Jong Hee Park, 2012. ``Unified Method for Dynamic and
#'   Cross-Sectional Heterogeneity: Introducing Hidden Markov Panel
#'   Models.''  \emph{American Journal of Political Science}.56:
#'   1040-1054. <doi: 10.1111/j.1540-5907.2012.00590.x>
#'
#'   Siddhartha Chib. 1998. ``Estimation and comparison of multiple
#'   change-point models.'' \emph{Journal of Econometrics}. 86: 221-241.
#'   <doi: 10.1016/S0304-4076(97)00115-2>
#'
#'   Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#'   ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical
#'   Software}. 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#'   ## data generating
#'   set.seed(1974)
#'   N <- 30
#'   T <- 80
#'   NT <- N*T
#'
#'   ## true parameter values
#'   true.beta <- c(1, 1)
#'   true.sigma <- 3
#'   x1 <- rnorm(NT)
#'   x2 <- runif(NT, 2, 4)
#'
#'   ## group-specific breaks
#'   break.point = rep(T/2, N); break.sigma=c(rep(1, N));
#'   break.list <- rep(1, N)
#'
#'   X <- as.matrix(cbind(x1, x2), NT, );
#'   y <- rep(NA, NT)
#'   id  <-  rep(1:N, each=NT/N)
#'   K <-  ncol(X);
#'   true.beta <- as.matrix(true.beta, K, 1)
#'
#'   ## compute the break probability
#'   ruler <- c(1:T)
#'   W.mat <- matrix(NA, T, N)
#'   for (i in 1:N){
#'     W.mat[, i] <- pnorm((ruler-break.point[i])/break.sigma[i])
#'   }
#'   Weight <- as.vector(W.mat)
#'
#'   ## draw time-varying individual effects and sample y
#'   j = 1
#'   true.sigma.alpha <- 30
#'   true.alpha1 <- true.alpha2 <- rep(NA, N)
#'   for (i in 1:N){
#'     Xi <- X[j:(j+T-1), ]
#'     true.mean <- Xi  %*% true.beta
#'     weight <- Weight[j:(j+T-1)]
#'     true.alpha1[i] <- rnorm(1, 0, true.sigma.alpha)
#'     true.alpha2[i] <- -1*true.alpha1[i]
#'     y[j:(j+T-1)] <- ((1-weight)*true.mean + (1-weight)*rnorm(T, 0, true.sigma) +
#'     		    (1-weight)*true.alpha1[i]) +
#'     		    (weight*true.mean + weight*rnorm(T, 0, true.sigma) + weight*true.alpha2[i])
#'     j <- j + T
#'   }
#'
#'   ## extract the standardized residuals from the OLS with fixed-effects
#'   FEols <- lm(y ~ X + as.factor(id) -1 )
#'   resid.all <- rstandard(FEols)
#'   time.id <- rep(1:80, N)
#'
#'   ## model fitting
#'   G <- 100
#'   BF <- testpanelSubjectBreak(subject.id=id, time.id=time.id,
#'          resid= resid.all, max.break=3, minimum = 10,
#'          mcmc=G, burnin = G, thin=1, verbose=G,
#'          b0=0, B0=1/100, c0=2, d0=2, Time = time.id)
#'
#'   ## get the estimated break numbers
#'   estimated.breaks <- make.breaklist(BF, threshold=3)
#'
#'   ## model fitting
#'   out <- HMMpanelFE(subject.id = id, y, X=X, m =  estimated.breaks,
#'              mcmc=G, burnin=G, thin=1, verbose=G,
#'              b0=0, B0=1/100, c0=2, d0=2, delta0=0, Delta0=1/100)
#'
#'   ## print out the slope estimate
#'   ## true values are 1 and 1
#'   summary(out)
#'
#'   ## compare them with the result from the constant fixed-effects
#'   summary(FEols)
#' }
#'
"HMMpanelFE" <-
  function(subject.id, y, X, m,
           mcmc=1000, burnin=1000, thin=1, verbose=0,
           b0=0, B0=0.001, c0 = 0.001, d0 = 0.001, delta0=0, Delta0=0.001,
           a = NULL, b = NULL, seed = NA, ...){
    ## m should be a vector with a number of breaks for each group
    ## id is a numeric list
    ## p is a lag order
    ## offset is the first time period from which each group starts

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    Y <- as.matrix(y);
    X <- as.matrix(cbind(1, X)) ## the intercept is not reported
    N <-  nrow(Y);
    K <-  ncol(X)

    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    mvn.prior <- form.mvn.prior(delta0, Delta0, 1)
    delta0 <- mvn.prior[[1]]
    Delta0 <- mvn.prior[[2]]

    nstore <- mcmc/thin
    nsubj <- length(unique(subject.id))

    ## groupinfo matrix
    ## col1: subj ID, col2: offset (first time C indexing), col3: #time periods
    if (unique(subject.id)[1] != 1){
      stop("subject.id should start 1!")
    }
    subject.offset <- c(0, which(diff(sort(subject.id))==1)[-nsubj])
    ## col1: subj ID, col2: offset (C indexing), col3: #time periods in each subject
    nsubject.vec <- rep(NA, nsubj)
    for (i in 1:nsubj){
      nsubject.vec[i] <- sum(subject.id==unique(subject.id)[i])
    }
    subject.groupinfo <- cbind(unique(subject.id), subject.offset, nsubject.vec)

    ## maximum time length
    ntime <- max(nsubject.vec)
    m.max <- max(m)
    m.min <- min(m)

    ## prior inputs
    P0data <- NULL
    Pstart <- NULL
    for (i in 1:nsubj){
      if(m[i] == 0){
        P0current <- 1
        Pstartcurrent <- 1
      }
      else{
        P0current <- trans.mat.prior(m=m[i], n=nsubject.vec[i], a=a, b=b)
        Pstartcurrent <- trans.mat.prior(m=m[i], n=nsubject.vec[i], a=.9, b=.1)
      }
      P0data <- c(P0data, P0current)
      Pstart <- c(Pstart, Pstartcurrent)
    }

    ## starting values
    ols <- lm(Y~X-1)
    beta.start <- coef(ols)
    sigma2.start <- summary(ols)$sigma^2
    deltastart  <- NULL
    Ytilde <- Y - X%*%beta.start
    deltaformula <- Ytilde ~ 1 ## without intercept
    for (i in 1:nsubj){
      deltacurrent <- rep(as.vector(coef(lm(Ytilde ~ 1))), m[i] + 1)
      deltastart <- c(deltastart, deltacurrent)
    }

    ## Storage
    totalstates0 <- sum(m+1)
    betadraws0 <- matrix(0, nstore, K)
    deltadraws0 <- matrix(data=0, nstore, totalstates0)
    sigmadraws0 <- matrix(data=0, nstore, totalstates0)
    statedraws0 <- matrix(data=0, nstore, totalstates0)

    ## call C++ code to draw sample
    posterior <- .C("cHMMpanelFE",
                    deltadraws = as.double(deltadraws0),
                    sigmadraws = as.double(sigmadraws0),
                    statedraws = as.double(statedraws0),

                    betadraws = as.double(betadraws0),
                    betarow = as.integer(nrow(betadraws0)),
                    betacol = as.integer(ncol(betadraws0)),
                    totalstates = as.integer(totalstates0),

                    nsubj = as.integer(nsubj),
                    ntime = as.integer(ntime),
                    nobs = as.integer(N),
                    subjectid = as.integer(subject.id),

                    m = as.integer(m),
                    mmax = as.integer(m.max),
                    mmin = as.integer(m.min),

                    Ydata = as.double(Y),
                    Yrow = as.integer(nrow(Y)),
                    Ycol = as.integer(ncol(Y)),

                    Xdata = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),

                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    verbose = as.integer(verbose),

                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),

                    betastartdata = as.double(beta.start),
                    sigma2start = as.double(sigma2.start),
                    deltastartdata = as.double(deltastart),
                    deltastartrow = as.integer(length(deltastart)),

                    b0data = as.double(b0),
                    B0data = as.double(B0),

                    delta0data = as.double(delta0),
                    Delta0data = as.double(Delta0),

                    c0 = as.double(c0),
                    d0 = as.double(d0),

                    P0data = as.double(P0data),
                    P0row = as.integer(length(P0data)),
                    Pstartdata = as.double(Pstart),

                    subject_groupinfodata = as.double(subject.groupinfo),
                    PACKAGE="MCMCpack")

    ## pull together matrix and build MCMC object to return
    betadraws <- matrix(posterior$betadraws,
                        posterior$betarow,
                        posterior$betacol)

    sigma.samp <- as.mcmc(matrix(posterior$sigmadraws, nstore, totalstates0))
    delta.samp <- as.mcmc(matrix(posterior$deltadraws, nstore, totalstates0))
    state.samp <- as.mcmc(matrix(posterior$statedraws, nstore, totalstates0))

    ## output <- mcmc(betadraws, start=burnin+1, end=burnin+mcmc, thin=thin)
    output <- as.mcmc(betadraws[, -1])## drop the intercept
    attr(output, "title") <- "HMMpanelFE Sample"
    attr(output, "m")  <- m
    attr(output, "sigma") <- sigma.samp
    attr(output, "state") <- state.samp
    attr(output, "delta") <- delta.samp
    return(output)
  }
