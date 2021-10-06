################################
## Poisson Changepoint Model
##
## 07/14/2009 Jong Hee Park
################################

#' Markov Chain Monte Carlo for a Poisson Regression Changepoint Model
#'
#' This function generates a sample from the posterior distribution of a
#' Poisson regression model with multiple changepoints. The function uses the
#' Markov chain Monte Carlo method of Chib (1998).  The user supplies data and
#' priors, and a sample from the posterior distribution is returned as an mcmc
#' object, which can be subsequently analyzed with functions provided in the
#' coda package.
#'
#' \code{MCMCpoissonChange} simulates from the posterior distribution of a
#' Poisson regression model with multiple changepoints using the methods of
#' Chib (1998) and Fruhwirth-Schnatter and Wagner (2006).  The details of the
#' model are discussed in Park (2010).
#'
#' The model takes the following form:
#'
#' \deqn{y_t \sim \mathcal{P}oisson(\mu_t)}
#'
#' \deqn{\mu_t = x_t ' \beta_m,\;\; m = 1, \ldots, M}
#'
#' Where
#' \eqn{M} is the number of states and \eqn{\beta_m} is paramters
#' when a state is \eqn{m} at \eqn{t}.
#'
#' We assume Gaussian distribution for prior of \eqn{\beta}:
#'
#' \deqn{\beta_m \sim \mathcal{N}(b_0,B_0^{-1}),\;\; m = 1, \ldots, M}
#'
#' And:
#'
#' \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M} Where \eqn{M} is the number of states.
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
#' @param c0 \eqn{c_0} is the shape parameter for Gamma prior on
#' \eqn{\lambda} (the mean). When there is no covariate, this should be
#' provided by users. No default value is provided.
#'
#' @param d0 \eqn{d_0} is the scale parameter for Gamma prior on
#' \eqn{\lambda} (the mean). When there is no covariate, this should be
#' provided by users. No default value is provided.
#'
#' @param lambda.mu The mean of the Gamma prior on \eqn{\lambda}.
#' \eqn{sigma.mu} and \eqn{sigma.var} allow users to
#' choose the Gamma prior by choosing its mean and variance.
#'
#' @param lambda.var The variacne of the Gamma prior on \eqn{\lambda}.
#' \eqn{sigma.mu} and \eqn{sigma.var} allow users to
#' choose the Gamma prior by choosing its mean and variance.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burn-in.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0, the
#' iteration number and the posterior density samples are printed to the screen
#' every \code{verbose}th iteration.
#'
#' @param seed The seed for the random number generator.  If NA, current R
#' system seed is used.
#'
#' @param beta.start The starting values for the beta vector. This can either
#' be a scalar or a column vector with dimension equal to the number of betas.
#' The default value of NA will use draws from the Uniform distribution with
#' the same boundary with the data as the starting value. If this is a scalar,
#' that value will serve as the starting value mean for all of the betas. When
#' there is no covariate, the log value of means should be used.
#'
#' @param P.start The starting values for the transition matrix. A user should
#' provide a square matrix with dimension equal to the number of states. By
#' default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper
#' transition matrix for each raw except the last raw.
#'
#' @param marginal.likelihood How should the marginal likelihood be calculated?
#' Options are: \code{none} in which case the marginal likelihood will not be
#' calculated, and \code{Chib95} in which case the method of Chib (1995) is
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
#' @seealso \code{\link{MCMCbinaryChange}}, \code{\link{plotState}},
#' \code{\link{plotChangepoint}}
#'
#' @references Jong Hee Park. 2010. ``Structural Change in the U.S. Presidents'
#' Use of Force Abroad.'' \emph{American Journal of Political Science} 54:
#' 766-782. <doi:10.1111/j.1540-5907.2010.00459.x>
#'
#' Sylvia Fruhwirth-Schnatter and Helga Wagner 2006. ``Auxiliary Mixture
#' Sampling for Parameter-driven Models of Time Series of Counts with
#' Applications to State Space Modelling.'' \emph{Biometrika}. 93:827--841.
#'
#' Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point
#' models.'' \emph{Journal of Econometrics}. 86: 221-241.
#' <doi: 10.1016/S0304-4076(97)00115-2>
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.
#' <doi: 10.1080/01621459.1995.10476635>
#' 
#' @keywords models
#'
#' @examples
#'
#'     \dontrun{
#'     set.seed(11119)
#'     n <- 150
#'     x1 <- runif(n, 0, 0.5)
#'     true.beta1 <- c(1,  1)
#'     true.beta2 <- c(1,  -2)
#'     true.beta3 <- c(1,  2)
#'
#'     ## set true two breaks at (50, 100)
#'     true.s <- rep(1:3, each=n/3)
#'     mu1 <- exp(1 + x1[true.s==1]*1)
#'     mu2 <- exp(1 + x1[true.s==2]*-2)
#'     mu3 <- exp(1 + x1[true.s==3]*2)
#'
#'     y <- as.ts(c(rpois(n/3, mu1), rpois(n/3, mu2), rpois(n/3, mu3)))
#'     formula = y ~ x1
#'
#'     ## fit multiple models with a varying number of breaks
#'     model0 <-  MCMCpoissonChange(formula, m=0,
#'             mcmc = 1000, burnin = 1000, verbose = 500,
#'             b0 = rep(0, 2), B0 = 1/5*diag(2), marginal.likelihood = "Chib95")
#'     model1 <-  MCMCpoissonChange(formula, m=1,
#'             mcmc = 1000, burnin = 1000, verbose = 500,
#'             b0 = rep(0, 2), B0 = 1/5*diag(2), marginal.likelihood = "Chib95")
#'     model2 <-  MCMCpoissonChange(formula, m=2,
#'             mcmc = 1000, burnin = 1000, verbose = 500,
#'             b0 = rep(0, 2), B0 = 1/5*diag(2), marginal.likelihood = "Chib95")
#'     model3 <-  MCMCpoissonChange(formula, m=3,
#'             mcmc = 1000, burnin = 1000, verbose = 500,
#'             b0 = rep(0, 2), B0 = 1/5*diag(2), marginal.likelihood = "Chib95")
#'     model4 <-  MCMCpoissonChange(formula, m=4,
#'             mcmc = 1000, burnin = 1000, verbose = 500,
#'             b0 = rep(0, 2), B0 = 1/5*diag(2), marginal.likelihood = "Chib95")
#'     model5 <-  MCMCpoissonChange(formula, m=5,
#'             mcmc = 1000, burnin = 1000, verbose = 500,
#'             b0 = rep(0, 2), B0 = 1/5*diag(2), marginal.likelihood = "Chib95")
#'
#'     ## find the most reasonable one
#'     print(BayesFactor(model0, model1, model2, model3, model4, model5))
#'
#'     ## draw plots using the "right" model
#'     par(mfrow=c(attr(model2, "m") + 1, 1), mai=c(0.4, 0.6, 0.3, 0.05))
#'     plotState(model2, legend.control = c(1, 0.6))
#'     plotChangepoint(model2, verbose = TRUE, ylab="Density", start=1, overlay=TRUE)
#'
#'     ## No covariate case
#'     model2.1 <- MCMCpoissonChange(y ~ 1, m = 2, c0 = 2, d0 = 1,
#'              mcmc = 1000, burnin = 1000, verbose = 500,
#'              marginal.likelihood = "Chib95")
#'     print(BayesFactor(model2, model2.1))
#'     }
#'
"MCMCpoissonChange"<- function(formula, data = parent.frame(), m = 1,
                               b0 = 0, B0 = 1, a = NULL, b = NULL, c0 = NA, d0 = NA,
                               lambda.mu = NA, lambda.var = NA,
                               burnin = 1000, mcmc = 1000, thin = 1, verbose = 0,
                               seed = NA, beta.start = NA, P.start = NA, ## offset = NA,
                               marginal.likelihood = c("none", "Chib95"), ...) {

    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    k <- ncol(X)
    n <- length(y)
    n.arrival<-  y + 1
    NT      <-  max(n.arrival)
    tot.comp <-  n + sum(y)
    ns      <- m + 1

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

    if (k==1){
        if (!is.na(lambda.mu) && !is.na(lambda.var)) {
            c0 <- lambda.mu^2/lambda.var
            d0 <- lambda.mu/lambda.var
        }
        if ((is.na(c0)||is.na(d0))&&((is.na(lambda.mu)||is.na(lambda.var)))){
            stop("You have to provide a prior for lambda (c0 and d0 or lambda.mu and lambda.var) when there is no covariate.\n")
        }
    } else{
        c0 <- d0 <- 0
        mvn.prior <- form.mvn.prior(b0, B0, k)
        b0 <- mvn.prior[[1]]
        B0 <- mvn.prior[[2]]
    }



    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
        chib <- 1
    }
    if (m == 0){
        if (marginal.likelihood == "Chib95"){
            if (is.na(b0)||is.na(B0)){
                stop("You have to have a prior for beta (b0 and B0) when m = 0.\n")
            } else{
                output <- MCMCpoisson(formula, burnin = burnin, mcmc = mcmc,
                                      thin = thin, verbose = verbose, b0 = b0, B0 = B0,
                                      marginal.likelihood = "Laplace")
                cat("\n Chib95 method is not yet available for m = 0. Laplace method is used instead.")
            }
        }else {
            output <- MCMCpoisson(formula, burnin = burnin, mcmc = mcmc,
                                  thin = thin, verbose = verbose, b0 = b0, B0 = B0)
        }
    }else {
        ## prior
        A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)

        ## get initial values of tau from observed y
        Pstart <- check.P(P.start, m, a=a, b=b)
        betastart  <- beta.change.start(beta.start, ns, k, formula, family=poisson, data)
        if (k == 1){
            betastart <- exp(betastart)
        }
        taustart <- tau.initial(y, tot.comp)
        componentstart  <-  round(runif(tot.comp, 1, 5))
        ## if(is.na(offset)){
        ##   logoffset <- rep(0, length(y))
        ## }
        ## else{
        ##    if(length(offset) == length(y)){
        ##      logoffset <- log(offset)
        ##   }
        ##   else{
        ##     stop("\n The length of offset is not matched with y.")
        ##   }
        ##  }
        ##  print(offset)

        ## normal mixture weights
        wr  <-  c(0.2294, 0.2590, 0.2480, 0.1525, 0.0472)
        mr  <-  c(0.0982, -1.5320, -0.7433, 0.8303, -3.1428)
        sr  <-  sqrt(c(0.2401, 1.1872, 0.3782, 0.1920, 3.2375))

        ## call C++ code to draw sample
        posterior <- .C("cMCMCpoissonChange",
                        betaout = as.double(rep(0.0, nstore*ns*k)),
                        Pout = as.double(rep(0.0, nstore*ns*ns)),
                        psout = as.double(rep(0.0, n*ns)),
                        sout = as.double(rep(0.0, nstore*n)),
                        Ydata = as.double(y),
                        Yrow = as.integer(nrow(y)),
                        Ycol = as.integer(ncol(y)),
                        Xdata = as.double(X),
                        Xrow = as.integer(nrow(X)),
                        Xcol = as.integer(ncol(X)),
                        ## logoffset = as.double(logoffset),
                        m = as.integer(m),
                        burnin = as.integer(burnin),
                        mcmc = as.integer(mcmc),
                        thin = as.integer(thin),
                        verbose = as.integer(verbose),
                        betastart = as.double(betastart),
                        Pstart = as.double(Pstart),
                        taustart = as.double(taustart),
                        componentstart = as.double(componentstart),
                        a = as.double(a),
                        b = as.double(b),
                        c0 = as.double(c0),
                        d0 = as.double(d0),
                        lecuyer = as.integer(lecuyer),
                        seedarray = as.integer(seed.array),
                        lecuyerstream = as.integer(lecuyer.stream),
                        b0data = as.double(b0),
                        B0data = as.double(B0),
                        A0data = as.double(A0),
                        logmarglikeholder = as.double(0.0),
                        loglikeholder = as.double(0.0),
                        wrin = as.double(wr),
                        mrin = as.double(mr),
                        srin = as.double(sr),
                        chib = as.integer(chib),

                        PACKAGE="MCMCpack")

        ## get marginal likelihood if Chib95
        if (marginal.likelihood == "Chib95"){
            logmarglike <- posterior$logmarglikeholder
            loglike <- posterior$loglikeholder
        }
        else {
            logmarglike <- 0
            loglike <- 0
        }


      ## pull together matrix and build MCMC object to return
      beta.holder <- matrix(posterior$betaout, nstore, )
      P.holder    <- matrix(posterior$Pout, nstore, )
      s.holder    <- matrix(posterior$sout, nstore, )
      ps.holder   <- matrix(posterior$psout, n, )

      output <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output)  <- sapply(c(1:ns), function(i) { paste(xnames, "_regime", i, sep = "")})
      attr(output, "title") <- "MCMCpoissonChange Posterior Sample"
      attr(output, "formula") <- formula
      attr(output, "y")       <- y
      attr(output, "X")       <- X
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglike") <- loglike
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "s.store") <- s.holder
      attr(output, "P.store") <- P.holder
    }
    return(output)

  }
