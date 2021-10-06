####################################################################
## HMM Gaussian Panel Random Effects Model
## y_it = x_it'*beta + w_it*bi + e_it,
## e_it ~ N(0, sigma2)
##
## bi ~ N(0, D) : random effects coefficient
## D ~ IW(r0, R0) : covariance matrix for multiple random effects
## beta ~ N(b0, B0) : fixed effect coefficient
## sigma2 ~ IG(c0/2, d0/2) : random error
##
## written by Jong Hee Park 03/2009
## modified and integrated with other codes on 09/2011
######################################################################

#' Markov Chain Monte Carlo for the Hidden Markov Random-effects Model
#'
#' HMMpanelRE generates a sample from the posterior distribution of
#' the hidden Markov random-effects model discussed in Park (2011).
#' The code works for panel data with the same starting point.  The
#' sampling of panel parameters is based on Algorithm 2 of Chib and
#' Carlin (1999). This model uses a multivariate Normal prior for the
#' fixed effects parameters and varying individual effects, an
#' Inverse-Wishart prior on the random-effects parameters, an
#' Inverse-Gamma prior on the residual error variance, and Beta prior
#' for transition probabilities.  The user supplies data and priors,
#' and a sample from the posterior distribution is returned as an mcmc
#' object, which can be subsequently analyzed with functions provided
#' in the coda package.
#'
#' \code{HMMpanelRE} simulates from the random-effect hidden Markov
#' panel model introduced by Park (2011).
#'
#' The model takes the following form: \deqn{y_i = X_i \beta_m + W_i
#' b_i + \varepsilon_i\;\; m = 1, \ldots, M}{y_i = X_i * beta_m + W_i
#' * b_i + epsilon_i, m = 1,..., M.} Where each group \eqn{i} have
#' \eqn{k_i} observations.  Random-effects parameters are assumed to
#' be time-varying at the system level: \deqn{b_i \sim
#' \mathcal{N}_q(0, D_m)} \deqn{\varepsilon_i \sim \mathcal{N}(0,
#' \sigma^2_m I_{k_i})}
#'
#' And the errors: We assume standard, conjugate priors: \deqn{\beta
#' \sim \mathcal{N}_p(b0, B0)} And: \deqn{\sigma^{2} \sim
#' \mathcal{IG}amma(c0/2, d0/2)} And: \deqn{D \sim
#' \mathcal{IW}ishart(r0, R0)} See Chib and Carlin (1999) for more
#' details.
#'
#' And: \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}
#' Where \eqn{M} is the number of states.
#'
#' \emph{NOTE:} We do not provide default parameters for the priors on
#' the precision matrix for the random effects. When fitting one of
#' these models, it is of utmost importance to choose a prior that
#' reflects your prior beliefs about the random effects. Using the
#' \code{dwish} and \code{rwish} functions might be useful in choosing
#' these values.
#'
#' @param subject.id A numeric vector indicating the group number. It
#'   should start from 1.
#'
#' @param time.id A numeric vector indicating the time unit. It should
#'   start from 1.
#'
#' @param y The dependent variable
#'
#' @param X The model matrix of the fixed-effects
#'
#' @param W The model matrix of the random-effects. W should be a
#'   subset of X.
#'
#' @param m The number of changepoints.
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
#' @param r0 The shape parameter for the Inverse-Wishart prior on
#'   variance matrix for the random effects. Set r=q for an
#'   uninformative prior where q is the number of random effects
#'
#' @param R0 The scale matrix for the Inverse-Wishart prior on
#'   variance matrix for the random effects. This must be a square
#'   q-dimension matrix. Use plausible variance regarding random
#'   effects for the diagonal of R.
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
#' @param beta.start The starting values for the beta vector. This can
#'   either be a scalar or a column vector with dimension equal to the
#'   number of betas.  The default value of NA will use draws from the
#'   Uniform distribution with the same boundary with the data as the
#'   starting value. If this is a scalar, that value will serve as the
#'   starting value mean for all of the betas. When there is no
#'   covariate, the log value of means should be used.
#'
#' @param sigma2.start The starting values for \eqn{\sigma^2}. This
#'   can either be a scalar or a column vector with dimension equal to
#'   the number of states.
#'
#' @param D.start The starting values for the beta vector. This can
#'   either be a scalar or a column vector with dimension equal to the
#'   number of betas. The default value of NA will use draws from the
#'   Uniform distribution with the same boundary with the data as the
#'   starting value. If this is a scalar, that value will serve as the
#'   starting value mean for all of the betas. When there is no
#'   covariate, the log value of means should be used.
#'
#' @param P.start The starting values for the transition matrix. A
#'   user should provide a square matrix with dimension equal to the
#'   number of states. By default, draws from the \code{Beta(0.9,
#'   0.1)} are used to construct a proper transition matrix for each
#'   raw except the last raw.
#'
#' @param marginal.likelihood How should the marginal likelihood be
#'   calculated?  Options are: \code{none} in which case the marginal
#'   likelihood will not be calculated and \code{Chib95} in which case
#'   the method of Chib (1995) is used.
#'
#' @param ... further arguments to be passed
#'
#' @export
#'
#' @return An mcmc object that contains the posterior sample. This
#'   object can be summarized by functions provided by the coda
#'   package.  The object contains an attribute \code{prob.state}
#'   storage matrix that contains the probability of \eqn{state_i} for
#'   each period, and the log-marginal likelihood of the model
#'   (\code{logmarglike}).
#' 
#' @references Jong Hee Park, 2012. ``Unified Method for Dynamic and
#'   Cross-Sectional Heterogeneity: Introducing Hidden Markov Panel
#'   Models.''  \emph{American Journal of Political Science}.56:
#'   1040-1054. <doi: 10.1111/j.1540-5907.2012.00590.x>
#'
#'   Siddhartha Chib. 1998. ``Estimation and comparison of multiple
#'   change-point models.'' \emph{Journal of Econometrics}. 86:
#'   221-241. <doi: 10.1016/S0304-4076(97)00115-2>
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
#' ``MCMCpack: Markov Chain Monte Carlo in R.'', \emph{Journal of
#' Statistical Software}.  42(9): 1-21. \doi{10.18637/jss.v042.i09}.
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#'   ## data generating
#'   set.seed(1977)
#'   Q <- 3
#'   true.beta1   <-  c(1, 1, 1) ; true.beta2   <-  c(-1, -1, -1)
#'   true.sigma2 <-  c(2, 5); true.D1 <- diag(.5, Q); true.D2 <- diag(2.5, Q)
#'   N=30; T=100;
#'   NT <- N*T
#'   x1 <- runif(NT, 1, 2)
#'   x2 <- runif(NT, 1, 2)
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
#'   ## model fitting
#'   G <- 100
#'   b0  <- rep(0, K) ; B0  <- solve(diag(100, K))
#'   c0  <- 2; d0  <- 2
#'   r0  <- 5; R0  <- diag(c(1, 0.1, 0.1))
#'   subject.id <- c(rep(1:N, each=T))
#'   time.id <- c(rep(1:T, N))
#'   out1 <- HMMpanelRE(subject.id, time.id, y, X, W, m=1,
#'                      mcmc=G, burnin=G, thin=1, verbose=G,
#'                      b0=b0, B0=B0, c0=c0, d0=d0, r0=r0, R0=R0)
#'
#'   ## latent state changes
#'   plotState(out1)
#'
#'   ## print mcmc output
#'   summary(out1)
#'
#'
#'
#' }
#'
"HMMpanelRE" <-
  function(subject.id, time.id, y, X, W, m=1,
           mcmc=1000, burnin=1000, thin=1, verbose=0,
           b0=0, B0=0.001, c0 = 0.001, d0 = 0.001, r0, R0, a = NULL, b = NULL,
           seed = NA, beta.start = NA, sigma2.start = NA, D.start= NA, P.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...){

    cl <- match.call()
    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    ns <- m + 1
    Y <- as.matrix(y)
    X <- as.matrix(X)
    W <- as.matrix(W)
    formula <- Y ~ X-1
    ols <- lm(formula)

    N <-  nrow(Y)
    K <-  ncol(X)
    Q <-  ncol(W)

    nobs <- nrow(X)

    ## Sort Data based on time.id
    oldTSCS <- cbind(time.id, subject.id, y, X, W)
    newTSCS <- oldTSCS[order(oldTSCS[,1]),]
    YT <- as.matrix(newTSCS[,3])
    XT <- as.matrix(newTSCS[,4:(4+K-1)])
    WT <- as.matrix(newTSCS[,(4+K):(4+K+Q-1)])

    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    R0 <- as.matrix(R0)

    nstore <- mcmc/thin
    nsubj <- length(unique(subject.id))
    if (unique(subject.id)[1] != 1){
      stop("subject.id should start 1!")
    }
    ## subject.offset is the obs number from which a new subject unit starts
    subject.offset <- c(0, which(diff(sort(subject.id))==1)[-nsubj])
    ## col1: subj ID, col2: offset (C indexing), col3: #time periods in each subject
    nsubject.vec <- rep(NA, nsubj)
    for (i in 1:nsubj){
      nsubject.vec[i] <- sum(subject.id==unique(subject.id)[i])
    }
    subject.groupinfo <- cbind(unique(subject.id), subject.offset, nsubject.vec)


    ## time.groupinfo
    ## col1: time ID, col2: offset (C indexing), col3: # subjects in each time
    if(unique(time.id)[1] != 1){
      time.id <- time.id - unique(time.id)[1] + 1
      cat("time.id does not start from 1. So it is modified by subtracting the first unit of time.")
    }
    ntime <- max(nsubject.vec)## maximum time length
    ntime.vec <- rep(NA, ntime)
    for (i in 1:ntime){
      ntime.vec[i] <- sum(time.id==unique(time.id)[i])
    }
    ## time.offset is the obs number from which a new time unit starts when we stack data by time.id
    time.offset <- c(0, which(diff(sort(time.id))==1)[-ntime])
    time.groupinfo <- cbind(unique(time.id), time.offset, ntime.vec)


    ## prior inputs
    if (m > 0){
      P0 <- trans.mat.prior(m=m, n=ntime, a=a, b=b)
      ## initial values
      Pstart  <-  check.P(P.start, m, a=a, b=b)
    }
    else {
      Pstart <- P0 <- matrix(1, 1, 1)
    }
    if (is.na(beta.start[1])) {
      betastart  <- coef(ols)
    }
    else{
      betastart <- beta.start
    }
    if (is.na(sigma2.start[1])) {
      sigma2start <- summary(ols)$sigma^2
    }
    else{
      sigma2start <- sigma2.start
    }

    betadraws <- matrix(data=0, nstore, ns*K)
    sigmadraws <- matrix(data=0, nstore, ns)
    Ddraws <- matrix(data=0, nstore, ns*Q*Q)
    psdraws <- matrix(data=0, ntime, ns)
    sdraws <- matrix(data=0, nstore, ntime*ns)

    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)

    ## following MCMCregress, set chib as binary
    logmarglike <- loglik <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    ## call C++ code to draw sample
    posterior <- .C("cHMMpanelRE",
                    betadata = as.double(betadraws),
                    betarow = as.integer(nrow(betadraws)),
                    betacol = as.integer(ncol(betadraws)),
                    sigmadata = as.double(sigmadraws),
                    Ddata = as.double(Ddraws),
                    psout = as.double(psdraws),
                    sout = as.double(sdraws),

                    nsubj = as.integer(nsubj), ntime = as.integer(ntime), m = as.integer(m),
                    nobs = as.integer(nobs),
                    subjectid = as.integer(subject.id),
                    timeid = as.integer(time.id),

                    Ydata = as.double(Y), Yrow = as.integer(nrow(Y)), Ycol = as.integer(ncol(Y)),
                    Xdata = as.double(X), Xrow = as.integer(nrow(X)), Xcol = as.integer(ncol(X)),
                    Wdata = as.double(W), Wrow = as.integer(nrow(W)), Wcol = as.integer(ncol(W)),
                    YTdata = as.double(YT), XTdata = as.double(XT), WTdata = as.double(WT),
                    burnin = as.integer(burnin), mcmc = as.integer(mcmc),
                    thin = as.integer(thin), verbose = as.integer(verbose),

                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),

                    betastartdata = as.double(betastart),
                    sigma2start = as.double(sigma2start),
                    Pstart = as.double(Pstart),

                    b0data = as.double(b0), B0data = as.double(B0),
                    c0 = as.double(c0), d0 = as.double(d0),
                    r0 = as.integer(r0), R0data = as.double(R0),

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
    D.samp <- matrix(posterior$Ddata,
                     posterior$betarow,
                     Q*Q*ns)
    xnames <-  sapply(c(1:K), function(i){paste("beta", i, sep = "")})
    Dnames <-  sapply(c(1:(Q*Q)), function(i){paste("D", i, sep = "")})
    if (m == 0){
      output <- as.mcmc(cbind(beta.samp, sigma.samp, D.samp))
      names <- c(xnames, "sigma2", Dnames)
      varnames(output) <- as.list(names)
    }
    else{
      output1 <- mcmc(data=beta.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
      output2 <- mcmc(data=sigma.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
      output3 <- mcmc(data=D.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output1)  <- sapply(c(1:ns),
                                 function(i){
                                   paste(xnames, "_regime", i, sep = "")
                                 })
      varnames(output2)  <- sapply(c(1:ns),
                                   function(i){
                                     paste("sigma2_regime", i, sep = "")
                                   })
      varnames(output3)  <- sapply(c(1:ns),
                                   function(i){
                                     paste(Dnames, "_regime", i, sep = "")
                                   })

      output <- as.mcmc(cbind(output1, output2, output3))
      ps.holder   <- matrix(posterior$psout, ntime, ns)
      s.holder    <- matrix(posterior$sout, nstore, )
    }

    attr(output, "title") <- "HMMpanelRE Posterior Sample"
    attr(output, "call")   <- cl
    attr(output, "y")       <- y[1:ntime]
    attr(output, "X")       <- X[1:ntime, ]
    attr(output, "m")       <- m
    attr(output, "nsubj")   <- nsubj
    attr(output, "ntime")   <- ntime
    if (m > 0){
      attr(output, "s.store") <- s.holder
      attr(output, "prob.state") <- ps.holder/nstore
    }
    attr(output, "logmarglike") <- posterior$logmarglikeholder
    attr(output, "loglike") <- posterior$loglikeholder

    return(output)
  }
