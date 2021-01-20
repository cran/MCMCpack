####################################################################
## test subject-level breaks from panel residuals
##
## written by Jong Hee Park 03/2009
## modified and integrated with other codes by JHP 07/2011
## fixed a starting.id and ending.id
######################################################################

#' A Test for the Subject-level Break using a Unitivariate Linear Regression
#' Model with Breaks
#'
#' testpanelSubjectBreak fits a unitivariate linear regression model with
#' parametric breaks using panel residuals to test the existence of
#' subject-level breaks in panel residuals. The details are discussed in Park
#' (2011).
#'
#'
#' \code{testpanelSubjectBreak} fits a univariate linear regression model for
#' subject-level residuals from a panel model.  The details are discussed in
#' Park (2011).
#'
#' The model takes the following form:
#'
#' \deqn{e_{it} = \alpha_{im} + \varepsilon_{it}\;\; m = 1, \ldots, M}
#'
#' The errors are assumed to be time-varying at the subject level:
#'
#' \deqn{\varepsilon_{it} \sim \mathcal{N}(0, \sigma^2_{im})}
#'
#' We assume standard, semi-conjugate priors:
#'
#' \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}
#'
#' And:
#'
#' \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}
#'
#' Where \eqn{\beta} and \eqn{\sigma^{-2}} are assumed \emph{a priori}
#' independent.
#'
#' And:
#'
#' \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}
#'
#' Where \eqn{M} is the number of states.
#'
#' OLS estimates are used for starting values.
#'
#' @param subject.id A numeric vector indicating the group number. It should
#' start from 1.
#'
#' @param time.id A numeric vector indicating the time unit. It should start
#' from 1.
#'
#' @param resid A vector of panel residuals.
#'
#' @param max.break An upper bound of break numbers for the test.
#'
#' @param minimum A minimum length of time series for the test. The test will
#' skip a subject with a time series shorter than this.
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
#' @param Time Times of the observations. This will be used to find the time of
#' the first observations in panel residuals.
#'
#' @param ps.out If ps.out == TRUE, state probabilities are exported. If the
#' number of panel subjects is huge, users can turn it off to save memory.
#'
#' @param ... further arguments to be passed
#'
#' @return The returned object is a matrix containing log marginal likelihoods
#' for all HMMs.  The dimension of the returned object is the number of panel
#' subjects by max.break + 1.  If psout == TRUE, the returned object has an
#' array attribute \code{psout} containing state probabilities for all HMMs.
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
#'   G <- 1000
#'   BF <- testpanelSubjectBreak(subject.id=id, time.id=time.id,
#'          resid= resid.all, max.break=3, minimum = 10,
#'          mcmc=G, burnin = G, thin=1, verbose=G,
#'          b0=0, B0=1/100, c0=2, d0=2, Time = time.id)
#'
#'   ## estimated break numbers
#'   ## thresho
#'   estimated.breaks <- make.breaklist(BF, threshold=3)
#'
#'   ## print all posterior model probabilities
#'   print(attr(BF, "model.prob"))
#' }
#'
"testpanelSubjectBreak" <-
  function(subject.id, time.id, resid, max.break=2, minimum = 10,
           mcmc=1000, burnin=1000, thin=1, verbose=0,
           b0, B0, c0, d0, a = NULL, b = NULL, seed = NA,
           Time = NULL, ps.out = FALSE, ...){
    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    N <- length(subject.id)

    ## groupinfo matrix
    ## col1: subj ID, col2: offset (first time C indexing), col3: #time periods
    if (min(subject.id) != 1){
      stop("subject.id should start 1!")
    }
    if (min(time.id) != 1){
      stop("time.id should start 1!")
    }
    if (is.null(Time)){
      Time <- rep(N, 1)
    }
    NC <- length(unique(subject.id))
    time.list <- as.numeric(table(subject.id))

    ## Make a residula list
    resid.list <- as.list(rep(NA, NC))
    start <- 1; end <- 0
    for (i in 1:NC){
      end <- start + time.list[i] - 1
      resid.list[[i]] <- ts(resid[start:end], start=Time[start])
      start <- end + 1
    }
    ## Do the break analysis
    BFout <- matrix(NA, NC, max.break + 1)
    if (ps.out ==TRUE){
      psout <- NULL
    }
    else {
      psout <- array(NA, c(max(time.list), sum(2:(max.break+1)), NC))
    }
    for (i in 1:NC){
      residual <- resid.list[[i]]
      nk <- length(residual)
      out <- as.list(rep(NA, max.break))
      if(nk > minimum){
        for (k in 0:max.break){
          out[[k+1]] <- MCMCresidualBreakAnalysis(residual, m=k,
                                                  b0=b0, B0=B0, c0=c0, d0=d0, a=a, b=b,
                                                  burnin=burnin, mcmc=mcmc, thin=thin, verbose=verbose,
                                                  marginal.likelihood="Chib95")
          if (ps.out ==TRUE&k>0){
            if(k==1){
              start <- 1
            }
            else{
              start <- sum(2:k)+1
            }
            probstate <- attr(out[[k+1]], "prob.state")
            psout[1:length(probstate[,1]), start:(start+k), i] <- probstate
          }
          ## if no convergence diagnostic
          BFout[i, k+1] <- attr(out[[k+1]], "logmarglike")
        }
      }
      if (verbose > 0){
        cat("\n ------------------------------------------------------------- ")
        cat("\n Break analysis for subject=", i, "is just finished! \n")
      }
    }
    if (ps.out ==TRUE){
      attr(BFout, "psout") <- psout
    }
    model.prob.mat <- matrix(NA, NC, max.break + 1)
    for (i in 1:NC){
      model.prob <- exp(BFout[i, ])/sum(exp(BFout[i, ]))
      winner <- which.max(model.prob)
      if (verbose > 0){
        cat("\nPr(no residual break) for subject", i, "=",
            model.prob[1])
      }
      model.prob.mat[i,] <-  model.prob
    }
    attr(BFout, "model.prob") <- model.prob.mat
    return(BFout)
  }
