##########################################################################
## sample from the posterior distribution of a factor analysis model
## model in R using linked C++ code in Scythe.
##
## The model is:
##
## x*_i = \Lambda \phi_i + \epsilon_i,   \epsilon_i \sim N(0, \Psi)
##
## \lambda_{ij} \sim N(l0_{ij}, L0^{-1}_{ij})
## \phi_i \sim N(0,I)
##
## and x*_i is the latent variable formed from the observed ordinal
## variable in the usual (Albert and Chib, 1993) way and is equal to
## x_i when x_i is continuous. When x_j is ordinal \Psi_jj is assumed
## to be 1.
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
## 12/2/2003
## Revised to accommodate new spec 7/20/2004
## Minor bug fix regarding std.mean 6/25/2004
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################
#' Markov Chain Monte Carlo for Mixed Data Factor Analysis Model
#'
#' This function generates a sample from the posterior distribution of a mixed
#' data (both continuous and ordinal) factor analysis model. Normal priors are
#' assumed on the factor loadings and factor scores, improper uniform priors
#' are assumed on the cutpoints, and inverse gamma priors are assumed for the
#' error variances (uniquenesses). The user supplies data and parameters for
#' the prior distributions, and a sample from the posterior distribution is
#' returned as an mcmc object, which can be subsequently analyzed with
#' functions provided in the coda package.
#'
#' The model takes the following form:
#'
#' Let \eqn{i=1,\ldots,N} index observations and \eqn{j=1,\ldots,K}
#' index response variables within an observation. An observed
#' variable \eqn{x_{ij}} can be either ordinal with a total of
#' \eqn{C_j} categories or continuous.  The distribution of \eqn{X} is
#' governed by a \eqn{N \times K} matrix of latent variables \eqn{X^*}
#' and a series of cutpoints \eqn{\gamma}. \eqn{X^*} is assumed to be
#' generated according to:
#'
#' \deqn{x^*_i = \Lambda \phi_i + \epsilon_i}
#'
#' \deqn{\epsilon_i \sim \mathcal{N}(0,\Psi)}
#'
#' where \eqn{x^*_i} is the \eqn{k}-vector of latent variables
#' specific to observation \eqn{i}, \eqn{\Lambda} is the \eqn{k \times
#' d} matrix of factor loadings, and \eqn{\phi_i} is the
#' \eqn{d}-vector of latent factor scores. It is assumed that the
#' first element of \eqn{\phi_i} is equal to 1 for all \eqn{i}.
#'
#' If the \eqn{j}th variable is ordinal, the probability that it takes the
#' value \eqn{c} in observation \eqn{i} is:
#'
#' \deqn{\pi_{ijc} = \Phi(\gamma_{jc} - \Lambda'_j\phi_i) -
#' \Phi(\gamma_{j(c-1)} - \Lambda'_j\phi_i)}
#'
#' If the \eqn{j}th variable is continuous, it is assumed that \eqn{x^*_{ij}
#' = x_{ij}} for all \eqn{i}.
#'
#' The implementation used here assumes independent conjugate priors for each
#' element of \eqn{\Lambda} and each \eqn{\phi_i}. More
#' specifically we assume:
#'
#' \deqn{\Lambda_{ij} \sim \mathcal{N}(l_{0_{ij}}, L_{0_{ij}}^{-1}),
#' i=1,\ldots,k, j=1,\ldots,d}
#'
#' \deqn{\phi_{i(2:d)} \sim \mathcal{N}(0, I), i=1,\dots,n}
#'
#' \code{MCMCmixfactanal} simulates from the posterior distribution using a
#' Metropolis-Hastings within Gibbs sampling algorithm. The algorithm employed
#' is based on work by Cowles (1996).  Note that the first element of
#' \eqn{\phi_i} is a 1. As a result, the first column of
#' \eqn{\Lambda} can be interpretated as negative item difficulty
#' parameters.  Further, the first element \eqn{\gamma_1} is
#' normalized to zero, and thus not returned in the mcmc object.  The
#' simulation proper is done in compiled C++ code to maximize efficiency.
#' Please consult the coda documentation for a comprehensive list of functions
#' that can be used to analyze the posterior sample.
#'
#' As is the case with all measurement models, make sure that you have plenty
#' of free memory, especially when storing the scores.
#'
#' @param x A one-sided formula containing the manifest variables. Ordinal
#' (including dichotomous) variables must be coded as ordered factors. Each
#' level of these ordered factors must be present in the data passed to the
#' function.  NOTE: data input is different in \code{MCMCmixfactanal} than in
#' either \code{MCMCfactanal} or \code{MCMCordfactanal}.
#'
#' @param factors The number of factors to be fitted.
#'
#' @param lambda.constraints List of lists specifying possible equality or
#' simple inequality constraints on the factor loadings. A typical entry in the
#' list has one of three forms: \code{varname=list(d,c)} which will constrain
#' the dth loading for the variable named varname to be equal to c,
#' \code{varname=list(d,"+")} which will constrain the dth loading for the
#' variable named varname to be positive, and \code{varname=list(d, "-")} which
#' will constrain the dth loading for the variable named varname to be
#' negative. If x is a matrix without column names defaults names of ``V1",
#' ``V2", ... , etc will be used. Note that, unlike \code{MCMCfactanal}, the
#' \eqn{\Lambda} matrix used here has \code{factors}+1 columns. The
#' first column of \eqn{\Lambda} corresponds to negative item
#' difficulty parameters for ordinal manifest variables and mean parameters for
#' continuous manifest variables and should generally not be constrained
#' directly by the user.
#'
#' @param data A data frame.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' iterations must be divisible by this value.
#'
#' @param tune The tuning parameter for the Metropolis-Hastings sampling. Can
#' be either a scalar or a \eqn{k}-vector (where \eqn{k} is the number of
#' manifest variables). \code{tune} must be strictly positive.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is great than 0 the
#' iteration number and the Metropolis-Hastings acceptance rate are printed to
#' the screen every \code{verbose}th iteration.
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
#' @param lambda.start Starting values for the factor loading matrix Lambda. If
#' \code{lambda.start} is set to a scalar the starting value for all
#' unconstrained loadings will be set to that scalar. If \code{lambda.start} is
#' a matrix of the same dimensions as Lambda then the \code{lambda.start}
#' matrix is used as the starting values (except for equality-constrained
#' elements). If \code{lambda.start} is set to \code{NA} (the default) then
#' starting values for unconstrained elements in the first column of Lambda are
#' based on the observed response pattern, the remaining unconstrained elements
#' of Lambda are set to 0, and starting values for inequality constrained
#' elements are set to either 1.0 or -1.0 depending on the nature of the
#' constraints.
#'
#' @param psi.start Starting values for the error variance (uniqueness) matrix.
#' If \code{psi.start} is set to a scalar then the starting value for all
#' diagonal elements of \code{Psi} that represent error variances for
#' continuous variables are set to this value. If \code{psi.start} is a
#' \eqn{k}-vector (where \eqn{k} is the number of manifest variables)
#' then the staring value of \code{Psi} has \code{psi.start} on the main
#' diagonal with the exception that entries corresponding to error variances
#' for ordinal variables are set to 1.. If \code{psi.start} is set to \code{NA}
#' (the default) the starting values of all the continuous variable
#' uniquenesses are set to 0.5. Error variances for ordinal response variables
#' are always constrained (regardless of the value of \code{psi.start} to have
#' an error variance of 1 in order to achieve identification.
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
#' @param store.lambda A switch that determines whether or not to store the
#' factor loadings for posterior analysis. By default, the factor loadings are
#' all stored.
#'
#' @param store.scores A switch that determines whether or not to store the
#' factor scores for posterior analysis.  \emph{NOTE: This takes an enormous
#' amount of memory, so should only be used if the chain is thinned heavily, or
#' for applications with a small number of observations}.  By default, the
#' factor scores are not stored.
#'
#' @param std.mean If \code{TRUE} (the default) the continuous manifest
#' variables are rescaled to have zero mean.
#'
#' @param std.var If \code{TRUE} (the default) the continuous manifest
#' variables are rescaled to have unit variance.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' 
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}},
#' \code{\link[stats]{factanal}}, \code{\link[MCMCpack]{MCMCfactanal}},
#' \code{\link[MCMCpack]{MCMCordfactanal}}, \code{\link[MCMCpack]{MCMCirt1d}},
#' \code{\link[MCMCpack]{MCMCirtKd}}
#'
#' @references Kevin M. Quinn. 2004. ``Bayesian Factor Analysis for Mixed
#' Ordinal and Continuous Responses.'' \emph{Political Analysis}. 12: 338-353.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' M. K. Cowles. 1996. ``Accelerating Monte Carlo Markov Chain Convergence for
#' Cumulative-link Generalized Linear Models." \emph{Statistics and Computing.}
#' 6: 101-110.
#'
#' Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling."
#' Springer: New York.
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
#' \dontrun{
#' data(PErisk)
#'
#' post <- MCMCmixfactanal(~courts+barb2+prsexp2+prscorr2+gdpw2,
#'                         factors=1, data=PErisk,
#'                         lambda.constraints = list(courts=list(2,"-")),
#'                         burnin=5000, mcmc=1000000, thin=50,
#'                         verbose=500, L0=.25, store.lambda=TRUE,
#'                         store.scores=TRUE, tune=1.2)
#' plot(post)
#' summary(post)
#'
#'
#'
#'
#' library(MASS)
#' data(Cars93)
#' attach(Cars93)
#' new.cars <- data.frame(Price, MPG.city, MPG.highway,
#'                  Cylinders, EngineSize, Horsepower,
#'                  RPM, Length, Wheelbase, Width, Weight, Origin)
#' rownames(new.cars) <- paste(Manufacturer, Model)
#' detach(Cars93)
#'
#' # drop obs 57 (Mazda RX 7) b/c it has a rotary engine
#' new.cars <- new.cars[-57,]
#' # drop 3 cylinder cars
#' new.cars <- new.cars[new.cars$Cylinders!=3,]
#' # drop 5 cylinder cars
#' new.cars <- new.cars[new.cars$Cylinders!=5,]
#'
#' new.cars$log.Price <- log(new.cars$Price)
#' new.cars$log.MPG.city <- log(new.cars$MPG.city)
#' new.cars$log.MPG.highway <- log(new.cars$MPG.highway)
#' new.cars$log.EngineSize <- log(new.cars$EngineSize)
#' new.cars$log.Horsepower <- log(new.cars$Horsepower)
#'
#' new.cars$Cylinders <- ordered(new.cars$Cylinders)
#' new.cars$Origin    <- ordered(new.cars$Origin)
#'
#'
#'
#' post <- MCMCmixfactanal(~log.Price+log.MPG.city+
#'                  log.MPG.highway+Cylinders+log.EngineSize+
#'                  log.Horsepower+RPM+Length+
#'                  Wheelbase+Width+Weight+Origin, data=new.cars,
#'                  lambda.constraints=list(log.Horsepower=list(2,"+"),
#'                  log.Horsepower=c(3,0), weight=list(3,"+")),
#'                  factors=2,
#'                  burnin=5000, mcmc=500000, thin=100, verbose=500,
#'                  L0=.25, tune=3.0)
#' plot(post)
#' summary(post)
#'
#' }
#'
"MCMCmixfactanal" <-
  function(x, factors, lambda.constraints=list(),
           data=parent.frame(), burnin = 1000, mcmc = 20000,
           thin=1, tune=NA, verbose = 0, seed = NA,
           lambda.start = NA, psi.start=NA,
           l0=0, L0=0, a0=0.001, b0=0.001,
           store.lambda=TRUE, store.scores=FALSE,
           std.mean=TRUE, std.var=TRUE, ... ) {

    call <- match.call()
    echo.name <- NULL
    mt <- terms(x, data=data)
    if (attr(mt, "response") > 0)
      stop("Response not allowed in formula in MCMCmixfactanal().\n")
    if(missing(data)) data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$factors <- mf$lambda.constraints <- mf$burnin <- mf$mcmc <- NULL
    mf$thin <- mf$tune <- mf$verbose <- mf$seed <- NULL
    mf$lambda.start <- mf$l0 <- mf$L0 <- mf$a0 <- mf$b0 <- NULL
    mf$store.lambda <- mf$store.scores <- mf$std.mean <- NULL
    mf$std.var <- mf$... <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf$na.action <- 'na.pass'
    mf <- eval(mf, sys.frame(sys.parent()))
    attributes(mt)$intercept <- 0
    Xterm.length <- length(attr(mt, "variables"))
    X <- subset(mf,
                select=as.character(attr(mt, "variables"))[2:Xterm.length])

    N <- nrow(X)	      # number of observations
    K <- ncol(X)              # number of manifest variables
    ncat <- matrix(NA, K, 1)  # vector of number of categ. in each man. var.
    for (i in 1:K){
      if (is.numeric(X[,i])){
        ncat[i] <- -999
        X[is.na(X[,i]), i] <- -999
      }
   else if (is.ordered(X[, i])) {
     ncat[i] <- nlevels(X[, i])
     temp <- as.integer(X[,i])
     temp <- ifelse(is.na(X[,i]) | (X[,i] == "<NA>"), -999, temp)
     X[, i] <- temp
   }
   else {
        stop("Manifest variable ", dimnames(X)[[2]][i],
             " neither ordered factor nor numeric variable.\n")
      }
    }

    X <- as.matrix(X)
    xvars <- dimnames(X)[[2]] # X variable names
    xobs <- dimnames(X)[[1]]  # observation names

    if (is.null(xobs)){
      xobs <- 1:N
    }

    # standardize X
    if (std.mean){
      for (i in 1:K){
        if (ncat[i] == -999){
          X[,i] <- X[,i]-mean(X[,i])
        }
      }
    }
    if (std.var){
      for (i in 1:K){
        if (ncat[i] == -999){
          X[,i] <- (X[,i] - mean(X[,i]))/sd(X[,i]) + mean(X[,i])
        }
      }
    }

    n.ord.ge3 <- 0
    for (i in 1:K)
      if (ncat[i] >= 3) n.ord.ge3 <- n.ord.ge3 + 1

    check.mcmc.parameters(burnin, mcmc, thin)

    ## setup constraints on Lambda
    holder <- build.factor.constraints(lambda.constraints, X, K, factors+1)
    Lambda.eq.constraints <- holder[[1]]
    Lambda.ineq.constraints <- holder[[2]]
    X.names <- holder[[3]]

    ## if subtracting out the mean of continuous X then constrain
    ## the mean parameter to 0
    for (i in 1:K){
      if (ncat[i] < 2 && std.mean==TRUE){
        if ((Lambda.eq.constraints[i,1] == -999 ||
             Lambda.eq.constraints[i,1] == 0.0) &&
            Lambda.ineq.constraints[i,1] == 0.0){
          Lambda.eq.constraints[i,1] <- 0.0
        }
        else {
          cat("Constraints on Lambda are logically\ninconsistent with std.mean==TRUE.\n")
          stop("Please respecify and call MCMCmixfactanal() again\n")
        }
      }
    }


    ## setup and check prior on Psi
    holder <- form.ig.diagmat.prior(a0, b0, K)
    a0 <- holder[[1]]
    b0 <- holder[[2]]

    ## setup prior on Lambda
    holder <- form.factload.norm.prior(l0, L0, K, factors+1, X.names)
    Lambda.prior.mean <- holder[[1]]
    Lambda.prior.prec <- holder[[2]]

    # seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    # Starting values for Lambda
    Lambda <- matrix(0, K, factors+1)
    if (is.na(lambda.start)){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              if (j==1){
                if (ncat[i] < 2){
                  Lambda[i,j] <- mean(X[,i]!=-999)
                }
                if (ncat[i] == 2){
                  probit.out <- glm(as.factor(X[X[,i]!=-999,i])~1,
                                    family=binomial(link="probit"))
                  probit.beta <- coef(probit.out)
                  Lambda[i,j] <- probit.beta[1]
                }
                if (ncat[i] > 2){
                  polr.out <- polr(ordered(X[X[,i]!=-999,i])~1)
                  Lambda[i,j] <- -polr.out$zeta[1]*.588
                }
              }
            }
            if(Lambda.ineq.constraints[i,j]>0){
              Lambda[i,j] <- 1.0
            }
            if(Lambda.ineq.constraints[i,j]<0){
              Lambda[i,j] <- -1.0
            }
          }
          else Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }
    }
    else if (is.matrix(lambda.start)){
      if (nrow(lambda.start)==K && ncol(lambda.start)==(factors+1))
        Lambda  <- lambda.start
      else {
        cat("Starting values not of correct size for model specification.\n")
        stop("Please respecify and call ", echo.name, "() again\n")
      }
    }
    else if (length(lambda.start)==1 && is.numeric(lambda.start)){
      Lambda <- matrix(lambda.start, K, factors+1)
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }
    }
    else {
      cat("Starting values neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call ", echo.name, "() again\n")
    }

    # check MH tuning parameter
    if (is.na(tune)){
      tune <- matrix(NA, K, 1)
      for (i in 1:K){
        tune[i] <- abs(0.05/ncat[i])
      }
    }
    else if (is.double(tune)){
      tune <- matrix(abs(tune/ncat), K, 1)
    }

    # starting values for gamma (note: not changeable by user)
    if (max(ncat) <= 2){
      gamma <- matrix(0, 3, K)
    }
    else {
      gamma <- matrix(0, max(ncat)+1, K)
    }
    for (i in 1:K){
      if (ncat[i]<=2){
        gamma[1,i] <- -300
        gamma[2,i] <- 0
        gamma[3,i] <- 300
      }
      if(ncat[i] > 2) {
        polr.out <- polr(ordered(X[X[,i]!=-999,i])~1)
        gamma[1,i] <- -300
        gamma[2,i] <- 0
        gamma[3:ncat[i],i] <- (polr.out$zeta[2:(ncat[i]-1)] -
                               polr.out$zeta[1])*.588

        gamma[ncat[i]+1,i] <- 300
      }
    }

    ## starting values for Psi
    Psi <- factuniqueness.start(psi.start, X)
    for (i in 1:K){
      if (ncat[i] >= 2){
        Psi[i,i] <- 1.0
      }
    }

    # define holder for posterior sample
    if (store.scores == FALSE && store.lambda == FALSE){
      sample <- matrix(data=0, mcmc/thin, length(gamma)+K)
    }
    else if (store.scores == TRUE && store.lambda == FALSE){
      sample <- matrix(data=0, mcmc/thin, (factors+1)*N + length(gamma)+K)
    }
    else if(store.scores == FALSE && store.lambda == TRUE) {
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+length(gamma)+K)
    }
    else { # store.scores==TRUE && store.lambda==TRUE
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+(factors+1)*N +
                       length(gamma)+K)
    }

    accepts <- matrix(0, K, 1)

    # Call the C++ code to do the real work
    posterior <- NULL
    posterior <- .C("mixfactanalpost",
                    samdata = as.double(sample),
                    samrow = as.integer(nrow(sample)),
                    samcol = as.integer(ncol(sample)),
                    X = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    tune = as.double(tune),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    Lambda = as.double(Lambda),
                    Lambdarow = as.integer(nrow(Lambda)),
                    Lambdacol = as.integer(ncol(Lambda)),
                    gamma = as.double(gamma),
                    gammarow = as.integer(nrow(gamma)),
                    gammacol = as.integer(ncol(gamma)),
                    Psi = as.double(Psi),
                    Psirow = as.integer(nrow(Psi)),
                    Psicol = as.integer(ncol(Psi)),
                    ncat = as.integer(ncat),
                    ncatrow = as.integer(nrow(ncat)),
                    ncatcol = as.integer(ncol(ncat)),
                    Lameq = as.double(Lambda.eq.constraints),
                    Lameqrow = as.integer(nrow(Lambda.eq.constraints)),
                    Lameqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lamineq = as.double(Lambda.ineq.constraints),
                    Lamineqrow = as.integer(nrow(Lambda.ineq.constraints)),
                    Lamineqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lampmean = as.double(Lambda.prior.mean),
                    Lampmeanrow = as.integer(nrow(Lambda.prior.mean)),
                    Lampmeancol = as.integer(ncol(Lambda.prior.prec)),
                    Lampprec = as.double(Lambda.prior.prec),
                    Lampprecrow = as.integer(nrow(Lambda.prior.prec)),
                    Lamppreccol = as.integer(ncol(Lambda.prior.prec)),
                    a0 = as.double(a0),
                    a0row = as.integer(nrow(a0)),
                    a0col = as.integer(ncol(a0)),
                    b0 = as.double(b0),
                    b0row = as.integer(nrow(b0)),
                    b0col = as.integer(ncol(b0)),
                    storelambda = as.integer(store.lambda),
                    storescores = as.integer(store.scores),
                    accepts = as.integer(accepts),
                    acceptsrow = as.integer(nrow(accepts)),
                    acceptscol = as.integer(ncol(accepts)),
                    PACKAGE="MCMCpack"
                    )

    accepts <- matrix(posterior$accepts, posterior$acceptsrow,
                      posterior$acceptscol, byrow=TRUE)
    rownames(accepts) <- X.names
    colnames(accepts) <- ""
    cat("\n\nAcceptance rates:\n")
    print(t(accepts) / (posterior$burnin+posterior$mcmc), digits=2)

    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$samdata, posterior$samrow, posterior$samcol,
                     byrow=FALSE)
    output <- mcmc(data=sample,start=1, end=mcmc, thin=thin)

    par.names <- NULL
    if (store.lambda==TRUE){
      Lambda.names <- paste(paste("Lambda",
                                  rep(X.names,
                                      each=(factors+1)), sep=""),
                            rep(1:(factors+1),K), sep=".")
      par.names <- c(par.names, Lambda.names)
    }

    gamma.names <- paste(paste("gamma",
                               rep(0:(nrow(gamma)-1),
                                   each=K), sep=""),
                         rep(X.names,  nrow(gamma)), sep=".")
    par.names <- c(par.names, gamma.names)

    if (store.scores==TRUE){
      phi.names <- paste(paste("phi",
                               rep(xobs, each=(factors+1)), sep="."),
                         rep(1:(factors+1),(factors+1)), sep=".")
      par.names <- c(par.names, phi.names)
    }

    Psi.names <- paste("Psi", X.names, sep=".")
    par.names <- c(par.names, Psi.names)

    varnames(output) <- par.names

    # get rid of columns for constrained parameters
    output.df <- as.data.frame(as.matrix(output))
    output.var <- diag(var(output.df))
    output.df <- output.df[,output.var != 0]
    output <- mcmc(as.matrix(output.df), start=burnin+1, end=burnin+mcmc,
                   thin=thin)

    # add constraint info so this isn't lost
    attr(output, "constraints") <- lambda.constraints
    attr(output, "n.manifest") <- K
    attr(output, "n.factors") <- factors
    attr(output, "accept.rates") <- t(accepts) / (posterior$burnin+posterior$mcmc)
      attr(output,"title") <-
        "MCMCpack Mixed Data Factor Analysis Posterior Sample"

    return(output)

  }
