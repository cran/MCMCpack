#########################################################
## sample from the posterior distribution
## of a linear Gaussian model with multiple changepoints
## using linked C++ code in Scythe
##
## JHP 07/01/2007
## JHP 03/03/2009
#########################################################

#' Markov Chain Monte Carlo for a linear Gaussian Multiple Changepoint Model
#'
#' This function generates a sample from the posterior distribution of a linear
#' Gaussian model with multiple changepoints. The function uses the Markov
#' chain Monte Carlo method of Chib (1998).  The user supplies data and priors,
#' and a sample from the posterior distribution is returned as an mcmc object,
#' which can be subsequently analyzed with functions provided in the coda
#' package.
#'
#' \code{MCMCregressChange} simulates from the posterior distribution of the
#' linear regression model with multiple changepoints.
#'
#' The model takes the following form:
#'
#' \deqn{y_t=x_t ' \beta_i + I(s_t=i)\varepsilon_{t},\;\; i=1, \ldots, k}
#'
#' Where \eqn{k} is the number of states and \eqn{I(s_t=i)} is an
#' indicator function that becomes 1 when a state at \eqn{t} is
#' \eqn{i} and otherwise 0.
#'
#' The errors are assumed to be Gaussian in each regime:
#'
#' \deqn{I(s_t=i)\varepsilon_{t} \sim \mathcal{N}(0, \sigma^2_i)}
#'
#' We assume standard, semi-conjugate priors:
#'
#' \deqn{\beta_i \sim \mathcal{N}(b_0,B_0^{-1}),\;\; i=1, \ldots, k}
#'
#' And:
#'
#' \deqn{\sigma^{-2}_i \sim \mathcal{G}amma(c_0/2, d_0/2),\;\; i=1, \ldots, k}
#'
#' Where \eqn{\beta_i} and \eqn{\sigma^{-2}_i} are assumed \emph{a
#' priori} independent.
#'
#' The simulation proper is done in compiled C++ code to maximize efficiency.
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param m The number of changepoints.
#'
#' @param b0 The prior mean of \eqn{\beta}.  This can either be a scalar
#' or a column vector with dimension equal to the number of betas. If this
#' takes a scalar value, then that value will serve as the prior mean for all
#' of the betas.
#'
#' @param B0 The prior precision of \eqn{\beta}.  This can either be a
#' scalar or a square matrix with dimensions equal to the number of betas.  If
#' this takes a scalar value, then that value times an identity matrix serves
#' as the prior precision of beta. Default value of 0 is equivalent to an
#' improper uniform prior for beta.
#'
#' @param c0 \eqn{c_0/2} is the shape parameter for the inverse Gamma
#' prior on \eqn{\sigma^2} (the variance of the disturbances). The
#' amount of information in the inverse Gamma prior is something like that from
#' \eqn{c_0} pseudo-observations.
#'
#' @param d0 \eqn{d_0/2} is the scale parameter for the inverse Gamma
#' prior on \eqn{\sigma^2} (the variance of the disturbances). In
#' constructing the inverse Gamma prior, \eqn{d_0} acts like the sum of
#' squared errors from the \eqn{c_0} pseudo-observations.
#'
#' @param sigma.mu The mean of the inverse Gamma prior on
#' \eqn{\sigma^2}.  \eqn{sigma.mu} and
#' \eqn{sigma.var} allow users to choose the inverse Gamma prior by
#' choosing its mean and variance.
#'
#' @param sigma.var The variacne of the inverse Gamma prior on
#' \eqn{\sigma^2}.  \eqn{sigma.mu} and
#' \eqn{sigma.var} allow users to choose the inverse Gamma prior by
#' choosing its mean and variance.
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
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burnin.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the \eqn{\beta} vector, and the error variance are
#' printed to the screen every \code{verbose}th iteration.
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
#' @param beta.start The starting values for the \eqn{\beta} vector.
#' This can either be a scalar or a column vector with dimension equal to the
#' number of betas.  The default value of of NA will use the MLE estimate of
#' \eqn{\beta} as the starting value.  If this is a scalar, that value
#' will serve as the starting value mean for all of the betas.
#'
#' @param P.start The starting values for the transition matrix.  A user should
#' provide a square matrix with dimension equal to the number of states.  By
#' default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper
#' transition matrix for each raw except the last raw.
#'
#' @param random.perturb If TRUE, randomly sample hidden states whenever
#' regularly sampled hidden states have at least one single observation state
#' (SOS). SOS is a sign of overfitting in non-ergodic hidden Markov models.
#'
#' @param WAIC Compute the Widely Applicable Information Criterion (Watanabe
#' 2010).
#'
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated, and \code{Chib95} in which case the method of Chib (1995) is
#' used.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.  The object
#' contains an attribute \code{prob.state} storage matrix that contains the
#' probability of \eqn{state_i} for each period, the log-likelihood of
#' the model (\code{loglike}), and the log-marginal likelihood of the model
#' (\code{logmarglike}).
#'
#' @export
#'
#' @seealso \code{\link{plotState}}, \code{\link{plotChangepoint}}
#'
#' @references Jong Hee Park, 2012. ``Unified Method for Dynamic and
#'   Cross-Sectional Heterogeneity: Introducing Hidden Markov Panel
#'   Models.''  \emph{American Journal of Political Science}.56:
#'   1040-1054. <doi: 10.1111/j.1540-5907.2012.00590.x>
#' 
#' Sumio Watanabe. 2010. "Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory"
#' \emph{Journal of Machine Learning Research}. 11: 3571-3594.
#'
#' Siddhartha Chib. 1995. "Marginal Likelihood from the Gibbs Output."
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.
#' <doi: 10.1016/S0304-4076(97)00115-2>
#'
#' Siddhartha Chib. 1998. "Estimation and comparison of multiple change-point
#' models."  \emph{Journal of Econometrics}. 86: 221-241.
#' <doi: 10.1080/01621459.1995.10476635>
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
#' set.seed(1119)
#' n <- 100
#' x1 <- runif(n)
#' true.beta1 <- c(2, -2)
#' true.beta2 <- c(0,  2)
#' true.Sigma <- c(1, 2)
#' true.s <- rep(1:2, each=n/2)
#'
#' mu1 <- cbind(1, x1[true.s==1])%*%true.beta1
#' mu2 <- cbind(1, x1[true.s==2])%*%true.beta2
#'
#' y <- as.ts(c(rnorm(n/2, mu1, sd=sqrt(true.Sigma[1])), rnorm(n/2, mu2, sd=sqrt(true.Sigma[2]))))
#' formula=y ~ x1
#'
#' ols1 <- lm(y[true.s==1] ~x1[true.s==1])
#' ols2 <- lm(y[true.s==2] ~x1[true.s==2])
#'
#' ## prior
#' b0 <- 0
#' B0 <- 0.1
#' sigma.mu=sd(y)
#' sigma.var=var(y)
#'
#' ## models
#' model0 <-  MCMCregressChange(formula, m=0, b0=b0, B0=B0, mcmc=100, burnin=100,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model1 <-  MCMCregressChange(formula, m=1, b0=b0, B0=B0, mcmc=100, burnin=100,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model2 <-  MCMCregressChange(formula, m=2, b0=b0, B0=B0, mcmc=100, burnin=100,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model3 <-  MCMCregressChange(formula, m=3, b0=b0, B0=B0, mcmc=100, burnin=100,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model4 <-  MCMCregressChange(formula, m=4, b0=b0, B0=B0, mcmc=100, burnin=100,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model5 <-  MCMCregressChange(formula, m=5, b0=b0, B0=B0, mcmc=100, burnin=100,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#'
#' print(BayesFactor(model0, model1, model2, model3, model4, model5))
#' plotState(model1)
#' plotChangepoint(model1)
#'
#' }
#'
"MCMCregressChange"<-
  function(formula, data=parent.frame(), m = 1,
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, sigma.mu = NA, sigma.var = NA,
           a = NULL, b = NULL,
           mcmc = 1000, burnin = 1000,  thin = 1, verbose = 0,
           seed = NA, beta.start = NA, P.start = NA,  random.perturb = FALSE,
           WAIC = FALSE, marginal.likelihood = c("none", "Chib95"), ...){

    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    k <- ncol(X)    # number of covariates
    n <- length(y)
    ns <- m + 1 # number of states

    ## check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    cl <- match.call()
    nstore <- mcmc/thin

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    if(!is.na(seed)) set.seed(seed)

    ## prior
    mvn.prior <- form.mvn.prior(b0, B0, k)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    if (is.na(sigma.mu)|is.na(sigma.var)) {
      check.ig.prior(c0, d0)
    }
    else {
      c0 <- 4 + 2 *(sigma.mu^2/sigma.var)
      d0 <- 2*sigma.mu *(c0/2 - 1)
    }

    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)

    ## following MCMCregress, set chib as binary
    logmarglike <- loglik <- NULL
    marginalrun <- ifelse(marginal.likelihood == "Chib95", 1, 0)
    sos <- ifelse(random.perturb, 1, 0)

    if (m == 0){
      output <- MCMCregress(formula, data=data, burnin = burnin, mcmc = mcmc,
                            thin = thin, verbose = verbose,
                            b0 = b0, B0 = B0, c0 =c0, d0=d0,
                            marginal.likelihood = "Chib95")
      attr(output, "y") <- y
      attr(output, "m") <- m
    }
    else{
      if(k == 1){
        output <- MCMCresidualBreakAnalysis(y, data=data,  m = m,
                                            b0 = b0, B0 = B0, c0 = c0, d0 = d0,
                                            a = a, b = b,
                                            burnin = burnin, mcmc = mcmc, thin = thin, verbose = verbose,
                                            seed = seed, beta.start = beta.start, P.start = P.start,
                                            marginal.likelihood = marginal.likelihood,
                                            random.perturb = random.perturb)
      }
      else{
        ## initial values
        Pstart  <-  check.P(P.start, m, a=a, b=b)
        A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)
        betastart  <- beta.change.start(beta.start, ns, k, formula, family=gaussian, data)
        ols <- lm(y~X-1)
        Sigmastart <- rep(summary(ols)$sigma^2, ns)
        statestart <- sort(sample(1:ns, n, replace=T))

        ## call C++ code to draw sample
        posterior <- .C("cMCMCregressChange",
                        betaout = as.double(rep(0.0, nstore*ns*k)),
                        Sigmaout = as.double(rep(0.0, nstore*ns)),
                        ## Pout = as.double(rep(0.0, nstore*ns*ns)),
                        psout = as.double(rep(0.0, n*ns)),
                        sout = as.double(rep(0.0, nstore*n)),
                        yloglike = as.double(rep(0.0, nstore*n)),


                        Ydata = as.double(y),
                        Yrow = as.integer(nrow(y)),
                        Ycol = as.integer(ncol(y)),
                        Xdata = as.double(X),
                        Xrow = as.integer(nrow(X)),
                        Xcol = as.integer(ncol(X)),

                        m = as.integer(m),
                        burnin = as.integer(burnin),
                        mcmc = as.integer(mcmc),
                        thin = as.integer(thin),
                        verbose = as.integer(verbose),

                        lecuyer=as.integer(lecuyer),
                        seedarray=as.integer(seed.array),
                        lecuyerstream=as.integer(lecuyer.stream),

                        betastart = as.double(betastart),
                        Sigmastart = as.double(Sigmastart),
                        Pstart = as.double(Pstart),
                        statestart = as.integer(statestart),

                        a = as.double(a),
                        b = as.double(b),
                        b0data = as.double(b0),
                        B0data = as.double(B0),
                        c0 = as.double(c0),
                        d0 = as.double(d0),
                        A0data = as.double(A0),
                        logmarglikeholder = as.double(0.0),
                        loglikeholder = as.double(0.0),
                        marginalrun = as.integer(marginalrun),
                        sos = as.integer(sos))

        ## get marginal likelihood if Chib95
        if (marginal.likelihood == "Chib95"){
            logmarglike <- posterior$logmarglikeholder
            loglike <- posterior$loglikeholder
        }
        Waic.out <- NA
        if(WAIC){
            Y.loglike.mat <- matrix(posterior$yloglike, nstore, n)
            Waic.out <- waic(Y.loglike.mat)$total
            rm(Y.loglike.mat)
            ## cat("    Waic: ", Waic.out[1], "\n")
        }


        ## pull together matrix and build MCMC object to return
        beta.holder <- matrix(posterior$betaout, nstore, ns*k)
        Sigma.holder <- matrix(posterior$Sigmaout, nstore, ns)
        ## P.holder    <- matrix(posterior$Pout, nstore, )
        s.holder    <- matrix(posterior$sout, nstore, )
        ps.holder   <- matrix(posterior$psout, n, )

        output1 <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
        varnames(output1)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c(xnames), "_regime", i, sep = "")
                                     })
        output2 <- mcmc(data=Sigma.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
        varnames(output2)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c("sigma2"), "_regime", i, sep = "")
                                     })
        output <- as.mcmc(cbind(output1, output2))

        attr(output, "title") <- "MCMCregressChange Posterior Sample"
        attr(output, "formula") <- formula
        attr(output, "y")       <- y
        attr(output, "X")       <- X
        attr(output, "m")       <- m
        attr(output, "call")    <- cl
        attr(output, "prob.state") <- ps.holder/nstore
        ## attr(output, "P.store") <- P.holder
        attr(output, "s.store") <- s.holder
        attr(output, "logmarglike") <- logmarglike
        attr(output, "loglik") <- loglik
        attr(output, "Waic.out") <- Waic.out
      }
    }
    return(output)

  }## end of MCMC function
