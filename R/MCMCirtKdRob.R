##########################################################################
## sample from a K-dimensional two-parameter item response model with
## logit link that has been rescaled so that the inverse link is:
##
## \delta0 + (1 - \delta0 - \delta1)*\Phi(.)
##
## where \delta0 \in (0, k0) and \delta1 \in (0, k1)
##
## priors for deltas are rescaled beta with parameters c0, d0, and c1, d1
##
##
## datamatrix is assumed to be nsubjects by nitems
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
## Feb. 17, 2005
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


#' Markov Chain Monte Carlo for Robust K-Dimensional Item Response Theory Model
#'
#' This function generates a posterior sample from a Robust K-dimensional item
#' response theory (IRT) model with logistic link, independent standard normal
#' priors on the subject abilities (ideal points), and independent normal
#' priors on the item parameters.  The user supplies data and priors, and a
#' sample from the posterior distribution is returned as an mcmc object, which
#' can be subsequently analyzed with functions provided in the coda package.
#'
#' \code{MCMCirtKdRob} simulates from the posterior using the slice sampling
#' algorithm of Neal (2003).  The simulation proper is done in compiled C++
#' code to maximize efficiency.  Please consult the coda documentation for a
#' comprehensive list of functions that can be used to analyze the posterior
#' sample.
#'
#' The model takes the following form.  We assume that each subject has an
#' subject ability (ideal point) denoted \eqn{\theta_j} \eqn{(K \times
#' 1)}, and that each item has a scalar difficulty parameter
#' \eqn{\alpha_i} and discrimination parameter \eqn{\beta_i}
#' \eqn{(K \times 1)}.  The observed choice by subject \eqn{j} on
#' item \eqn{i} is the observed data matrix which is \eqn{(I \times J)}.
#'
#' The probability that subject \eqn{j} answers item \eqn{i} correctly is
#' assumed to be:
#'
#' \deqn{\pi_{ij} = \delta_0 + (1 - \delta_0 - \delta_1) / (1+\exp(\alpha_i - \beta_i \theta_j))}
#'
#' This model was discussed in Bafumi et al. (2005).
#'
#' We assume the following priors.  For the subject abilities (ideal points) we
#' assume independent standard Normal priors:
#'
#' \deqn{\theta_{j,k} \sim \mathcal{N}(0,1)}
#'
#' These cannot be changed by the user.  For each item parameter, we
#' assume independent Normal priors:
#'
#' \deqn{\left[\alpha_i, \beta_i \right]' \sim \mathcal{N}_{(K+1)} (b_{0,i},B_{0,i})}
#'
#' Where \eqn{B_{0,i}} is a diagonal matrix.  One can specify a
#' separate prior mean and precision for each item parameter. We also
#' assume \eqn{\delta_0 / k_0 \sim }\eqn{
#' \mathcal{B}eta(c_0, d_0)} and
#' \eqn{\delta_1 / k_1 \sim }\eqn{
#' \mathcal{B}eta(c_1, d_1)}.
#'
#' The model is identified by constraints on the item parameters and / or
#' ability parameters. See Rivers (2004) for a discussion of identification of
#' IRT models.
#'
#' As is the case with all measurement models, make sure that you have plenty
#' of free memory, especially when storing the item parameters.
#'
#' @param datamatrix The matrix of data.  Must be 0, 1, or missing values.  It
#' is of dimensionality subjects by items.
#'
#' @param dimensions The number of dimensions in the latent space.
#'
#' @param item.constraints List of lists specifying possible equality
#'   or simple inequality constraints on the item parameters. A
#'   typical entry in the list has one of three forms:
#'   \code{rowname=list(d,c)} which will constrain the dth item
#'   parameter for the item named rowname to be equal to c,
#'   \code{rowname=list(d,"+")} which will constrain the dth item
#'   parameter for the item named rowname to be positive, and
#'   \code{rowname=list(d, "-")} which will constrain the dth item
#'   parameter for the item named rowname to be negative. If
#'   datamatrix is a matrix without row names defaults names of
#'   ``V1", ``V2", ... , etc will be used. In a \eqn{K}-dimensional
#'   model, the first item parameter for item \eqn{i} is the
#'   difficulty parameter (\eqn{\alpha_i}), the second item parameter
#'   is the discrimation parameter on dimension 1 (\eqn{\beta_{i,1}}),
#'   the third item parameter is the discrimation parameter on
#'   dimension 2 (\eqn{\beta_{i,2}}), ..., and the \eqn{(K+1)}th item
#'   parameter is the discrimation parameter on dimension \eqn{K}
#'   (\eqn{\beta_{i,K}}).  The item difficulty parameters
#'   (\eqn{\alpha}) should generally not be constrained.
#'
#' @param ability.constraints List of lists specifying possible equality or
#' simple inequality constraints on the ability parameters. A typical entry in
#' the list has one of three forms: \code{colname=list(d,c)} which will
#' constrain the dth ability parameter for the subject named colname to be
#' equal to c, \code{colname=list(d,"+")} which will constrain the dth ability
#' parameter for the subject named colname to be positive, and
#' \code{colname=list(d, "-")} which will constrain the dth ability parameter
#' for the subject named colname to be negative. If datamatrix is a matrix
#' without column names defaults names of ``V1", ``V2", ... , etc will be used.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of iterations for the sampler after burn-in.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' iterations must be divisible by this value.
#'
#' @param interval.method Method for finding the slicing interval. Can be equal
#' to either \code{step} in which case the stepping out algorithm of Neal
#' (2003) is used or \code{doubling} in which case the doubling procedure of
#' Neal (2003) is used. The stepping out method tends to be faster on a
#' per-iteration basis as it typically requires few function calls. The
#' doubling method expands the initial interval more quickly which makes the
#' Markov chain mix somewhat more quickly-- at least in some situations.
#'
#' @param theta.w The initial width of the slice sampling interval for each
#' ability parameter (the elements of \eqn{\theta})
#'
#' @param theta.mp The parameter governing the maximum possible width of the
#' slice interval for each ability parameter (the elements of
#' \eqn{\theta}). If \code{interval.method="step"} then the maximum
#' width is \code{theta.w * theta.mp}.
#'
#' If \code{interval.method="doubling"} then the maximum width is \code{theta.w
#' * 2^theta.mp}.
#'
#' @param alphabeta.w The initial width of the slice sampling interval for each
#' item parameter (the elements of \eqn{\alpha} and \eqn{\beta})
#'
#' @param alphabeta.mp The parameter governing the maximum possible width of
#' the slice interval for each item parameters (the elements of
#' \eqn{\alpha} and \eqn{\beta}). If \code{interval.method="step"}
#' then the maximum width is \code{alphabeta.w * alphabeta.mp}.
#'
#' If \code{interval.method="doubling"} then the maximum width is
#' \code{alphabeta.w * 2^alphabeta.mp}.
#'
#' @param delta0.w The initial width of the slice sampling interval for
#' \eqn{\delta_0}
#'
#' @param delta0.mp The parameter governing the maximum possible width of the
#' slice interval for \eqn{\delta_0}. If \code{interval.method="step"}
#' then the maximum width is \code{delta0.w * delta0.mp}. If
#' \code{interval.method="doubling"} then the maximum width is \code{delta0.w *
#' 2^delta0.mp}.
#'
#' @param delta1.w The initial width of the slice sampling interval for
#' \eqn{\delta_1}
#'
#' @param delta1.mp The parameter governing the maximum possible width of the
#' slice interval for \eqn{\delta_1}. If \code{interval.method="step"}
#' then the maximum width is \code{delta1.w * delta1.mp}. If
#' \code{interval.method="doubling"} then the maximum width is \code{delta1.w *
#' 2^delta1.mp}.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If verbose > 0, the iteration number with
#' be printed to the screen every verbose'th iteration.
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
#' @param theta.start The starting values for the ability parameters
#' \eqn{\theta}. Can be either a scalar or a matrix with number of rows
#' equal to the number of subjects and number of columns equal to the dimension
#' \eqn{K} of the latent space. If \code{theta.start=NA} then starting
#' values will be chosen that are 0 for unconstrained subjects, -0.5 for
#' subjects with negative inequality constraints and 0.5 for subjects with
#' positive inequality constraints.
#'
#' @param alphabeta.start The starting values for the \eqn{\alpha} and
#' \eqn{\beta} difficulty and discrimination parameters. If
#' \code{alphabeta.start} is set to a scalar the starting value for all
#' unconstrained item parameters will be set to that scalar. If
#' \code{alphabeta.start} is a matrix of dimension \eqn{(K+1) \times
#' items} then the \code{alphabeta.start} matrix is used as the
#' starting values (except for equality-constrained elements). If
#' \code{alphabeta.start} is set to \code{NA} (the default) then starting
#' values for unconstrained elements are set to values generated from a series
#' of proportional odds logistic regression fits, and starting values for
#' inequality constrained elements are set to either 1.0 or -1.0 depending on
#' the nature of the constraints.
#'
#' @param delta0.start The starting value for the \eqn{\delta_0}
#' parameter.
#'
#' @param delta1.start The starting value for the \eqn{\delta_1}
#' parameter.
#'
#' @param b0 The prior means of the \eqn{\alpha} and \eqn{\beta}
#' difficulty and discrimination parameters, stacked for all items.  If a
#' scalar is passed, it is used as the prior mean for all items.
#'
#' @param B0 The prior precisions (inverse variances) of the independent Normal
#' prior on the item parameters.  Can be either a scalar or a matrix of
#' dimension \eqn{(K+1) \times items}.
#'
#' @param k0 \eqn{\delta_0} is constrained to lie in the interval
#' between 0 and \code{k0}.
#'
#' @param k1 \eqn{\delta_1} is constrained to lie in the interval
#' between 0 and \code{k1}.
#'
#' @param c0 Parameter governing the prior for \eqn{\delta_0}.
#' \eqn{\delta_0} divided by \code{k0} is assumed to be follow a beta
#' distribution with first parameter \code{c0}.
#'
#' @param d0 Parameter governing the prior for \eqn{\delta_0}.
#' \eqn{\delta_0} divided by \code{k0} is assumed to be follow a beta
#' distribution with second parameter \code{d0}.
#'
#' @param c1 Parameter governing the prior for \eqn{\delta_1}.
#' \eqn{\delta_1} divided by \code{k1} is assumed to be follow a beta
#' distribution with first parameter \code{c1}.
#'
#' @param d1 Parameter governing the prior for \eqn{\delta_1}.
#' \eqn{\delta_1} divided by \code{k1} is assumed to be follow a beta
#' distribution with second parameter \code{d1}.
#'
#' @param store.item A switch that determines whether or not to store the item
#' parameters for posterior analysis.  \emph{NOTE: This typically takes an
#' enormous amount of memory, so should only be used if the chain is thinned
#' heavily, or for applications with a small number of items}.  By default, the
#' item parameters are not stored.
#'
#' @param store.ability A switch that determines whether or not to store the
#' subject abilities for posterior analysis.  By default, the item parameters
#' are all stored.
#'
#' @param drop.constant.items A switch that determines whether or not items
#' that have no variation should be deleted before fitting the model. Default =
#' TRUE.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MCMCpack]{MCMCirt1d}}, \code{\link[MCMCpack]{MCMCirtKd}}
#'
#' @references James H. Albert. 1992. ``Bayesian Estimation of Normal Ogive
#' Item Response Curves Using Gibbs Sampling." \emph{Journal of Educational
#' Statistics}.  17: 251-269.
#'
#' Joseph Bafumi, Andrew Gelman, David K. Park, and Noah Kaplan. 2005.
#' ``Practical Issues in Implementing and Understanding Bayesian Ideal Point
#' Estimation.'' \emph{Political Analysis}.
#'
#' Joshua Clinton, Simon Jackman, and Douglas Rivers. 2004. ``The Statistical
#' Analysis of Roll Call Data."  \emph{American Political Science Review}.  98:
#' 355-370.
#'
#' Simon Jackman. 2001. ``Multidimensional Analysis of Roll Call Data via
#' Bayesian Simulation.'' \emph{Political Analysis.} 9: 227-241.
#'
#' Valen E. Johnson and James H. Albert. 1999. \emph{Ordinal Data Modeling}.
#' Springer: New York.
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
#' Douglas Rivers.  2004.  ``Identification of Multidimensional Item-Response
#' Models."  Stanford University, typescript.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
#'    ## Court example with ability (ideal point) and
#'    ##  item (case) constraints
#'    data(SupremeCourt)
#'    post1 <- MCMCirtKdRob(t(SupremeCourt), dimensions=1,
#'                    burnin=500, mcmc=5000, thin=1,
#'                    B0=.25, store.item=TRUE, store.ability=TRUE,
#'                    ability.constraints=list("Thomas"=list(1,"+"),
#'                    "Stevens"=list(1,-4)),
#'                    item.constraints=list("1"=list(2,"-")),
#'                    verbose=50)
#'    plot(post1)
#'    summary(post1)
#'
#'    ## Senate example with ability (ideal point) constraints
#'    data(Senate)
#'    namestring <- as.character(Senate$member)
#'    namestring[78] <- "CHAFEE1"
#'    namestring[79] <- "CHAFEE2"
#'    namestring[59] <- "SMITH.NH"
#'    namestring[74] <- "SMITH.OR"
#'    rownames(Senate) <- namestring
#'    post2 <- MCMCirtKdRob(Senate[,6:677], dimensions=1,
#'                          burnin=1000, mcmc=5000, thin=1,
#'                          ability.constraints=list("KENNEDY"=list(1,-4),
#'                                   "HELMS"=list(1, 4), "ASHCROFT"=list(1,"+"),
#'                                   "BOXER"=list(1,"-"), "KERRY"=list(1,"-"),
#'                                   "HATCH"=list(1,"+")),
#'                          B0=0.1, store.ability=TRUE, store.item=FALSE,
#'                          verbose=5, k0=0.15, k1=0.15,
#'                          delta0.start=0.13, delta1.start=0.13)
#'
#'    plot(post2)
#'    summary(post2)
#'    }
#'
"MCMCirtKdRob" <-
  function(datamatrix, dimensions, item.constraints=list(),
           ability.constraints=list(),
           burnin = 500, mcmc = 5000, thin=1, interval.method="step",
           theta.w=0.5, theta.mp=4, alphabeta.w=1.0, alphabeta.mp=4,
           delta0.w=NA, delta0.mp=3, delta1.w=NA, delta1.mp=3,
           verbose = FALSE, seed = NA, theta.start=NA,
           alphabeta.start = NA, delta0.start = NA,
           delta1.start = NA, b0=0, B0=0,
           k0=.1, k1=.1, c0=1, d0=1, c1=1, d1=1,
           store.item=TRUE, store.ability=FALSE,
           drop.constant.items=TRUE, ... ) {

    ## set X up
    if (drop.constant.items==TRUE){
      x.col.var <- apply(datamatrix, 2, var, na.rm=TRUE)
      keep.inds <- x.col.var>0
      keep.inds[is.na(keep.inds)] <- FALSE
      datamatrix <- datamatrix[,keep.inds]
    }
    X <- as.data.frame(datamatrix)
    xvars <- dimnames(X)[[2]]
    xobs <- dimnames(X)[[1]]
    N <- nrow(X)    # number of subjects
    K <- ncol(X)    # number of items
    for (i in 1:K){
      X[is.na(X[,i]), i] <- -999
    }
    if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) != (N*K)) {
      cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    X <- as.matrix(X)

    ## take care of the case where X has no row names
    if (is.null(xobs)){
      xobs <- 1:N
    }

    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    ## check slice sampling parameters
    if (!(interval.method %in% c("step", "doubling"))){
      cat("Error: interval.method not equal to 'step' or 'doubling'.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    method.step <- 0
    if (interval.method == "step"){
      method.step <- 1
    }
    if (theta.w <= 0 ){
      cat("Error: theta.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (theta.mp < 1 ){
      cat("Error: theta.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (alphabeta.w <= 0 ){
      cat("Error: alphabeta.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (alphabeta.mp < 1 ){
      cat("Error: alphabeta.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (is.na(delta0.w)){
      delta0.w <- 0.25*k0
    }
    if (delta0.w <= 0 ){
      cat("Error: delta0.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (delta0.mp < 1 ){
      cat("Error: delta0.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (is.na(delta1.w)){
      delta1.w <- 0.25*k1
    }
    if (delta1.w <= 0 ){
      cat("Error: delta1.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (delta1.mp < 1 ){
      cat("Error: delta1.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }

    ## error check the prior parameters for delta
    if (k0 < 0 | k0 > 0.5){
      cat("Error: k0 not in (0, 0.5).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (k1 < 0 | k1 > 0.5){
      cat("Error: k1 not in (0, 0.5).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (c0 < 0){
      cat("Error: c0 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (c1 < 0){
      cat("Error: c1 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (d0 < 0){
      cat("Error: d0 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (d1 < 0){
      cat("Error: d1 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }

    ## setup constraints on Lambda = (alpha, beta)
    holder <- build.factor.constraints(item.constraints, X, K, dimensions+1)
    Lambda.eq.constraints <- holder[[1]]
    Lambda.ineq.constraints <- holder[[2]]
    X.names <- holder[[3]]

    ## setup constraints on theta
    holder <- build.factor.constraints(ability.constraints, t(X),
                                       N, dimensions)
    theta.eq.constraints <- holder[[1]]
    theta.ineq.constraints <- holder[[2]]

    ## setup prior on Lambda
    holder <- form.factload.norm.prior(b0, B0, K, dimensions+1, X.names)
    Lambda.prior.mean <- holder[[1]]
    Lambda.prior.prec <- holder[[2]]

    # seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Starting values for delta0 and delta1
    if (is.na(delta0.start)){
      delta0.start <- 0.5 * k0;
    }
    if (is.na(delta1.start)){
      delta1.start <- 0.5 * k1;
    }
    if (delta0.start < 0 | delta0.start > k0){
      cat("Error: delta0 not in (0, k0).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    if (delta1.start < 0 | delta1.start > k1){
      cat("Error: delta1 not in (0, k1).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }



    ## Starting values for Lambda
    Lambda <- matrix(0, K, dimensions+1)
    if (is.na(alphabeta.start)){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:(dimensions+1)){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              if (j==1){
                  probit.out <- glm(as.factor(X[X[,i]!=-999,i])~1,
                                    family=binomial(link="logit"))
                  probit.beta <- coef(probit.out)
                  Lambda[i,j] <- -1 * probit.beta[1]
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
    else if (is.matrix(alphabeta.start)){
      if (nrow(alphabeta.start)==K && ncol(alphabeta.start)==(dimensions+1))
        Lambda  <- alphabeta.start
      else {
        cat("Starting values not of correct size for model specification.\n")
        stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)
      }
    }
    else if (length(alphabeta.start)==1 && is.numeric(alphabeta.start)){
      Lambda <- matrix(alphabeta.start, K, dimensions+1)
      for (i in 1:K){
        for (j in 1:(dimensions+1)){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }
    }
    else {
      cat("Starting values for alpha & beta neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)
    }
    for (i in 1:K){
      lam.sqdist <- sum(Lambda[i,]^2)
      while (lam.sqdist > 100){
        Lambda[i,] <- Lambda[i,] * 0.95
        lam.sqdist <- sum(Lambda[i,]^2)
      }
    }



    ## Starting values for theta
    if (is.na(theta.start)){
      theta <- matrix(0, N, dimensions)
    }
    else if(is.null(dim(theta.start)) & is.numeric(theta.start)){
      theta <- matrix(theta.start, N, dimensions)
    }
    else if(nrow(theta.start)==N & ncol(theta.start)==dimensions){
      theta <- theta.start
    }
    else{
      cat("Starting values for theta neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)
    }
    for (i in 1:N){
      for (j in 1:dimensions){
        if (theta.eq.constraints[i,j]==-999){
          if(theta.ineq.constraints[i,j]>0){
            theta[i,j] <- 0.5
          }
          if(theta.ineq.constraints[i,j]<0){
            theta[i,j] <- -0.5
          }
        }
        else theta[i,j] <- theta.eq.constraints[i,j]
      }
    }





    ## define holder for posterior sample
    if (store.ability == FALSE && store.item == FALSE){
      cat("You need to store either the ability or item parameters.\n")
      stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)
    }
    else if (store.ability == TRUE && store.item == FALSE){
      sample <- matrix(data=0, mcmc/thin, (dimensions+1)*N+2)
    }
    else if(store.ability == FALSE && store.item == TRUE) {
      sample <- matrix(data=0, mcmc/thin, K*(dimensions+1)+2)
    }
    else { # store.ability==TRUE && store.item==TRUE
      sample <- matrix(data=0, mcmc/thin, K*(dimensions+1)+(dimensions+1)*N+2)
    }


    ## Call the C++ code to do the real work
    posterior <- .C("irtKdRobpost",
                    samdata = as.double(sample),
                    samrow = as.integer(nrow(sample)),
                    samcol = as.integer(ncol(sample)),
                    X = as.integer(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    method.step = as.integer(method.step),
                    theta.w = as.double(theta.w),
                    theta.mp = as.integer(theta.mp),
                    ab.w = as.double(alphabeta.w),
                    ab.mp = as.integer(alphabeta.mp),
                    delta0.w = as.double(delta0.w),
                    delta0.mp = as.integer(delta0.mp),
                    delta1.w = as.double(delta1.w),
                    delta1.mp = as.integer(delta1.mp),
                    delta0 = as.double(delta0.start),
                    delta1 = as.double(delta1.start),
                    Lambda = as.double(Lambda),
                    Lambdarow = as.integer(nrow(Lambda)),
                    Lambdacol = as.integer(ncol(Lambda)),
                    theta = as.double(theta),
                    thetarow = as.integer(nrow(theta)),
                    thetacol = as.integer(ncol(theta)),
                    Lameq = as.double(Lambda.eq.constraints),
                    Lameqrow = as.integer(nrow(Lambda.eq.constraints)),
                    Lameqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lamineq = as.double(Lambda.ineq.constraints),
                    Lamineqrow = as.integer(nrow(Lambda.ineq.constraints)),
                    Lamineqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    theteq = as.double(theta.eq.constraints),
                    theteqrow = as.integer(nrow(theta.eq.constraints)),
                    theteqcol = as.integer(ncol(theta.ineq.constraints)),
                    thetineq = as.double(theta.ineq.constraints),
                    thetineqrow = as.integer(nrow(theta.ineq.constraints)),
                    thetineqcol = as.integer(ncol(theta.ineq.constraints)),
                    Lampmean = as.double(Lambda.prior.mean),
                    Lampmeanrow = as.integer(nrow(Lambda.prior.mean)),
                    Lampmeancol = as.integer(ncol(Lambda.prior.prec)),
                    Lampprec = as.double(Lambda.prior.prec),
                    Lampprecrow = as.integer(nrow(Lambda.prior.prec)),
                    Lamppreccol = as.integer(ncol(Lambda.prior.prec)),
                    k0 = as.double(k0),
                    k1 = as.double(k1),
                    c0 = as.double(c0),
                    c1 = as.double(c1),
                    d0 = as.double(d0),
                    d1 = as.double(d1),
                    storeitem = as.integer(store.item),
                    storesability = as.integer(store.ability),
                    PACKAGE="MCMCpack"
                    )



    ## put together matrix and build MCMC object to return
    sample <- matrix(posterior$samdata, posterior$samrow, posterior$samcol,
                     byrow=FALSE)
    output <- mcmc(data=sample,start=1, end=mcmc, thin=thin)

    par.names <- NULL
    if (store.item==TRUE){
      alpha.hold <- paste("alpha", X.names, sep=".")
      beta.hold <- paste("beta", X.names, sep = ".")
      beta.hold <- rep(beta.hold, dimensions, each=dimensions)
      beta.hold <- paste(beta.hold, 1:dimensions, sep=".")

      Lambda.names <- t(cbind(matrix(alpha.hold, K, 1),
                              matrix(beta.hold,K,dimensions,byrow=TRUE)))
      dim(Lambda.names) <- NULL
      par.names <- c(par.names, Lambda.names)
    }

    if (store.ability==TRUE){
      phi.names <- paste(paste("theta",
                               rep(xobs, each=(dimensions+1)), sep="."),
                         rep(0:dimensions,(dimensions+1)), sep=".")
      par.names <- c(par.names, phi.names)
    }

    par.names <- c("delta0", "delta1", par.names)

    varnames(output) <- par.names

    ## get rid of columns for constrained parameters
    output.df <- as.data.frame(as.matrix(output))
    output.var <- diag(var(output.df))
    output.df <- output.df[,output.var != 0]
    output <- mcmc(as.matrix(output.df), start=1, end=mcmc, thin=thin)

    ## add constraint info so this isn't lost
    attr(output, "constraints") <- item.constraints
    attr(output, "n.items") <- K
    attr(output, "n.dimensions") <- dimensions
    attr(output,"title") <-
      "MCMCpack Robust K-Dimensional Item Response Theory Model Posterior Sample"

return(output)


  }
