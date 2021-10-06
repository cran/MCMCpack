##########################################################################
## sample from the posterior distribution of a 2-dimensional pairwise
##    comparisons model with Dirichlet process prior in R using
##    linked C++ code in Scythe.
##
## Model is:
##
##  i = 1,...,n.resp  (resondents)
##  j = 1,...,n.cand  (candidates)
##
## Y_{ijj'} = 1 if i chooses j over j'
## Y_{ijj'} = 0 if i chooses j' over j
## Y_{ijj'} = NA if i chooses neither;
##
## Pr(Y_{ijj'} = 1) = \Phi( z_{i}^{T} [\theta_{j} - \theta_{ j'} ] )
## z_{i}=[cos(gamma_{i}) sin(gamma_{i})]^{T}
## 
## gamma_i has a Dirichlet Process prior.
## We use a finite stick breaking process to approximate the infinite
## Dirichlet Process.
##
## theta_j \overset{ind}{\sim} N_{2}(0, I_{2})
##            (some theta_js truncated above or below 0, or fixed to constants)
##
##
## candidate IDs in columns 2 to 4 need to begin with a letter
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Based on 1-d code from KQ 3/17/2015
## Initial 2-d DP code QY 11/14/2019
## Various improvements and modifications KQ and QY 2019 to 2021
## Added to MCMCpack KQ 6/26/2021
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


#' Markov Chain Monte Carlo for the Two-Dimensional Pairwise Comparisons
#' Model with Dirichlet Process Prior in Yu and Quinn (2021)
#'
#' This function generates a sample from the posterior distribution of a
#' model for pairwise comparisons data with a probit link. Unlike standard
#' models for pairwise comparisons data, in this model the latent attribute
#' of each item being compared is a vector in two-dimensional Euclidean space.
#'
#' \code{MCMCpaircompare2d} uses the data augmentation approach of Albert and
#' Chib (1993) in conjunction with Gibbs and Metropolis-within-Gibbs steps
#' to fit the model. The user supplies data and a sample from the
#' posterior is returned as an \code{mcmc} object, which can be subsequently
#' analyzed in the \code{coda} package.
#'
#' The simulation is done in compiled C++ code to maximize efficiency.
#'
#' Please consult the \code{coda} package documentation for a comprehensive
#' list of functions that can be used to analyze the posterior sample.
#'
#' The model takes the following form:
#' 
#' \deqn{i = 1,...,I \ \ \ \   (raters) }
#' \deqn{j = 1,...,J \ \ \ \   (items)  }
#' \deqn{Y_{ijj'} = 1 \ \  if  \ \  i \ \  chooses \ \  j \ \ over \ \ j'}
#' \deqn{Y_{ijj'} = 0 \ \ if \ \ i \ \ chooses \ \ j' \ \ over \ \ j}
#' \deqn{Y_{ijj'} = NA \ \  if \ \ i \ \ chooses \ \ neither}
#'
#' \deqn{\Pr(Y_{ijj'} = 1) = \Phi( \mathbf{z}_{i}' [\boldsymbol{\theta}_{j} -
#' \boldsymbol{\theta}_{ j'} ])}
#' \deqn{\mathbf{z}_{i}=[\cos(\gamma_{i}), \  \sin(\gamma_{i})]' }
#' 
#' The following priors are assumed:
#' \deqn{\gamma_i \sim G}
#' \deqn{G \sim \mathcal{DP}(\alpha G_0)}
#' \deqn{G_0 = \mathcal{U}nif(0,   \pi/2)}
#' \deqn{\alpha \sim \mathcal{G}amma(a_0, b_0)} 
#' \deqn{\boldsymbol{\theta}_j \sim
#' \mathcal{N}_{2}(\mathbf{0}, \mathbf{I}_{2})}
#' For identification, some \eqn{\boldsymbol{\theta}_j}s are truncated
#' above or below 0, or fixed to constants.
#'
#'
#' @param pwc.data A data.frame containing the pairwise comparisons data.
#' Each row of \code{pwc.data} corresponds to a single pairwise comparison.
#' \code{pwc.data} needs to have exactly four columns. The first column
#' contains a unique identifier for the rater. Column two contains the unique
#' identifier for the first item being compared. Column three contains the
#' unique identifier for the second item being compared. Column four contains
#' the unique identifier of the item selected from the two items being
#' compared. If a tie occurred, the entry in the fourth column should be NA.
#' \strong{The identifiers in columns 2 through 4 must start with a letter.
#' Examples are provided below.} 
#'
#' @param theta.constraints A list specifying possible simple equality or
#' inequality constraints on the item parameters. A
#'   typical entry in the list has one of three forms:
#'   \code{itemname=list(d,c)} which will constrain the dth dimension of
#'   theta for the item named \code{itemname} to be equal to c,
#'   \code{itemname=list(d,"+")} which will constrain the dth dimension of
#'   theta for the item named \code{itemname} to be positive, and
#'   \code{itemname=list(d, "-")} which will constrain the dth dimension of
#'   theta for the item named \code{itemname} to be negative. 
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of Gibbs iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' Gibbs iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 
#' output is printed to the screen every
#' \code{verbose}th iteration.
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
#' @param gamma.start The starting value for the gamma vector.  This
#' can either be a scalar or a column vector with dimension equal to the number
#' of raters.  If this takes a scalar value, then that value will serve as the
#' starting value for all of the gammas. The default value of NA will set the
#' starting value of each gamma parameter to \eqn{\pi/4}.
#'
#'
#' @param theta.start Starting values for the theta. Can be either a numeric
#' scalar, a J by 2 matrix (where J is the number of items compared), or NA.
#' If a scalar, all theta values are set to that value (except elements already
#' specified via theta.contraints. If NA, then non constrained elements of
#' theta are set equal to 0, elements constrained to be positive are set equal
#' to 0.5, elements constrained to be negative are set equal to -0.5 and
#' elements with equality constraints are set to satisfy those constraints. 
#'
#' @param store.theta Should the theta draws be returned? Default is TRUE. 
#'
#' @param store.gamma Should the gamma draws be returned? Default is TRUE.
#'
#' @param tune Tuning parameter for the random walk Metropolis proposal for
#' each gamma_i. \code{tune} is the width of the uniform proposal centered at
#' the current value of gamma_i. Must be a positive scalar.
#'
#' @param procrustes Should the theta and gamma draws be post-processed with
#' a Procrustes transformation? Default is FALSE. The Procrustes target matrix
#' is derived from the constrained elements of theta. Each row of theta that
#' has both theta values constrained is part of the of the target matrix.
#' Elements with equality constraints are set to those values. Elements
#' constrained to be positive are set to 1. Elements constrained to be negative
#' are set to -1. If \code{procrustes} is set to \code{TRUE} theta.constraints
#' must be set so that there are at least three rows of theta that have both
#' elements of theta constrained. 
#'
#' @param alpha.start The starting value for the DP concentration parameter
#' alpha. Must be a positive scalar. Defaults to 1. If \code{alpha.fixed} is
#' set equal to TRUE, then alpha is held fixed at alpha.start.
#'
#' @param cluster.max The maximum number of clusters allowed in the
#' approximation to the DP prior for gamma. Defaults to 100. Must be a
#' positive integer. 
#'
#' @param cluster.mcmc The number of additional MCMC iterations that are done
#' to sample each cluster-specific gamma value within one main MCMC iteration.
#' Must be a positive integer. Defaults to 500. Setting this to a lower value
#' speeds runtime at the cost of (possibly) worse mixing. 
#'
#' @param alpha.fixed Logical value indicating whether the DP concentration
#' parameter alpha be held fixed (\code{TRUE}) or estimated (\code{FALSE}). 
#'
#' @param a0 The shape parameter of the gamma prior for alpha. This is the
#' same parameterization of the gamma distribution as R's internal
#' \code{rgamma()} function. Only relevant if \code{alpha.fixed} is set equal
#' to \code{FALSE}. Defaults to 1.
#'
#' @param b0 The rate parameter of the gamma prior for alpha.  This is the
#' same parameterization of the gamma distribution as R's internal
#' \code{rgamma()} function. Only relevant if \code{alpha.fixed} is set equal
#' to \code{FALSE}. Defaults to 1.
#' 
#' @param ... further arguments to be passed
#'
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package. Most of the column
#' names of the mcmc object are self explanatory. Note however that the columns
#' with names of the form "cluster.[raterID]" give the cluster membership of
#' each rater at each stored MCMC iteration. Because of the possibility of
#' label switching, the particular values of these cluster membership variables
#' are not meaningful. What is meaningful is whether two raters share the same
#' cluster membership value at a particular MCMC iteration. This indicates
#' that those two raters were clustered together during that iteration.
#' Finally, note that the "n.clusters" column gives the number of distinct
#' gamma values at each iteration, i.e. the number of clusters at that
#' iteration.
#'
#' @author Qiushi Yu <yuqiushi@umich.edu> and
#' Kevin M. Quinn <kmq@umich.edu>
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MCMCpack]{MCMCpaircompare}},
#' \code{\link[MCMCpack]{MCMCpaircompare2dDP}}
#'
#' @references Albert, J. H. and S. Chib. 1993. ``Bayesian Analysis of Binary
#' and Polychotomous Response Data.'' \emph{J. Amer. Statist. Assoc.} 88,
#' 669-679
#'
#' Yu, Qiushi and Kevin M. Quinn. 2021. ``A Multidimensional Pairwise
#' Comparison Model for Heterogeneous Perceptions with an Application to
#' Modeling the Perceived Truthfulness of Public Statements on COVID-19.''
#' University of Michigan Working Paper.
#' 
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
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
#' \dontrun{
#' ## a synthetic data example
#' set.seed(123)
#' 
#' I <- 65  ## number of raters
#' J <- 50 ## number of items to be compared
#' 
#' ## 3 clusters:
#' ## raters 1 to 5 put most weight on dimension 1
#' ## raters 6 to 10 put most weight on dimension 2
#' ## raters 11 to I put substantial weight on both dimensions
#' gamma.true <- c(rep(0.05, 5),
#'              rep(1.50, 5),
#'              rep(0.7, I-10) )
#' theta1.true <- rnorm(J, m=0, s=1)
#' theta2.true <- rnorm(J, m=0, s=1)
#' theta1.true[1] <- 2
#' theta2.true[1] <- 2
#' theta1.true[2] <- -2
#' theta2.true[2] <- -2
#' theta1.true[3] <-  2
#' theta2.true[3] <- -2
#' 
#' 
#' 
#' n.comparisons <- 125 ## number of pairwise comparisons for each rater
#' 
#' ## generate synthetic data according to the assumed model
#' rater.id <- NULL
#' item.1.id <- NULL
#' item.2.id <- NULL
#' choice.id <- NULL
#' for (i in 1:I){
#'     for (c in 1:n.comparisons){
#'         rater.id <- c(rater.id, i+100)
#'         item.numbers <- sample(1:J, size=2, replace=FALSE)
#'         item.1 <- item.numbers[1]
#'         item.2 <- item.numbers[2]
#'         item.1.id <- c(item.1.id, item.1)
#'         item.2.id <- c(item.2.id, item.2)
#'         z <- c(cos(gamma.true[i]), sin(gamma.true[i]))
#'         eta <- z[1] * (theta1.true[item.1] - theta1.true[item.2])  +
#'             z[2] * (theta2.true[item.1] - theta2.true[item.2])
#'         prob.item.1.chosen <- pnorm(eta)
#'         u <- runif(1)
#'         if (u <= prob.item.1.chosen){
#'             choice.id <- c(choice.id, item.1)
#'         }
#'         else{
#'             choice.id <- c(choice.id, item.2)
#'         }
#'     }
#' }
#' item.1.id <- paste("item", item.1.id+100, sep=".")
#' item.2.id <- paste("item", item.2.id+100, sep=".")
#' choice.id <- paste("item", choice.id+100, sep=".")
#' 
#' sim.data <- data.frame(rater.id, item.1.id, item.2.id, choice.id)
#' 
#' 
#' ## fit the model (should be run for more than 10500 iterations)
#' posterior <- MCMCpaircompare2dDP(pwc.data=sim.data,
#'                                  theta.constraints=list(item.101=list(1,2),
#'                                                         item.101=list(2,2),
#'                                                         item.102=list(1,-2),
#'                                                         item.102=list(2,-2),
#'                                                         item.103=list(1,"+"),
#'                                                         item.103=list(2,"-")),
#'                                  verbose=100,
#'                                  burnin=500, mcmc=10000, thin=5,
#'                                  cluster.mcmc=10,
#'                                  store.theta=TRUE, store.gamma=TRUE,
#'                                  tune=0.1)
#' 
#' 
#' 
#' 
#' 
#' theta1.draws <- posterior[, grep("theta1", colnames(posterior))]
#' theta2.draws <- posterior[, grep("theta2", colnames(posterior))]
#' gamma.draws <- posterior[, grep("gamma", colnames(posterior))]
#' 
#' theta1.post.med <- apply(theta1.draws, 2, median)
#' theta2.post.med <- apply(theta2.draws, 2, median)
#' gamma.post.med <- apply(gamma.draws, 2, median)
#' 
#' theta1.post.025 <- apply(theta1.draws, 2, quantile, prob=0.025)
#' theta1.post.975 <- apply(theta1.draws, 2, quantile, prob=0.975)
#' theta2.post.025 <- apply(theta2.draws, 2, quantile, prob=0.025)
#' theta2.post.975 <- apply(theta2.draws, 2, quantile, prob=0.975)
#' gamma.post.025 <- apply(gamma.draws, 2, quantile, prob=0.025)
#' gamma.post.975 <- apply(gamma.draws, 2, quantile, prob=0.975)
#' 
#' 
#' 
#' ## compare estimates to truth
#' par(mfrow=c(2,2))
#' plot(theta1.true, theta1.post.med, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5),
#'      col=rgb(0,0,0,0.3))
#' segments(x0=theta1.true, x1=theta1.true,
#'          y0=theta1.post.025, y1=theta1.post.975,
#'          col=rgb(0,0,0,0.3)) 
#' abline(0, 1, col=rgb(1,0,0,0.5))
#' 
#' plot(theta2.true, theta2.post.med, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5),
#'      col=rgb(0,0,0,0.3))
#' segments(x0=theta2.true, x1=theta2.true,
#'          y0=theta2.post.025, y1=theta2.post.975,
#'          col=rgb(0,0,0,0.3)) 
#' abline(0, 1, col=rgb(1,0,0,0.5))
#' 
#' plot(gamma.true, gamma.post.med, xlim=c(0, 1.6), ylim=c(0, 1.6),
#'      col=rgb(0,0,0,0.3))
#' segments(x0=gamma.true, x1=gamma.true,
#'          y0=gamma.post.025, y1=gamma.post.975,
#'          col=rgb(0,0,0,0.3)) 
#' abline(0, 1, col=rgb(1,0,0,0.5))
#' 
#' 
#' ## plot point estimates 
#' plot(theta1.post.med, theta2.post.med,
#'      xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5),
#'      col=rgb(0,0,0,0.3))
#' for (i in 1:length(gamma.post.med)){
#'     arrows(x0=0, y0=0,
#'            x1=cos(gamma.post.med[i]),
#'            y1=sin(gamma.post.med[i]),
#'            col=rgb(1,0,0,0.2), len=0.05, lwd=0.5)
#' }
#' 
#' }
#' @export


"MCMCpaircompare2dDP" <- function(pwc.data, theta.constraints=list(),
                                burnin=1000, mcmc=20000, thin=1,
                                verbose=0, seed=NA,
                                gamma.start=NA,
                                theta.start=NA,
                                store.theta=TRUE,
                                store.gamma=FALSE,
                                tune=0.3, procrustes=FALSE,
                                alpha.start=1, cluster.max=100,
                                cluster.mcmc=500,
                                alpha.fixed=TRUE, 
                                a0=1, b0=1,
                                ...){
    
    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    if (!(procrustes %in% c(TRUE, FALSE))){
        cat("procrustes cannot take a value other than TRUE or FALSE.\n")
        stop("Please check function arguments and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    
    ## convert all columns to character data
    pwc.data[,1] <- as.character(pwc.data[,1])
    pwc.data[,2] <- as.character(pwc.data[,2])
    pwc.data[,3] <- as.character(pwc.data[,3])
    pwc.data[,4] <- as.character(pwc.data[,4])
    
    ## check input data
    if (ncol(pwc.data) != 4){
        cat("pwc.data must have 4 columns. The specified pwc.data does not have 4 columns.\n")
        stop("Please check data and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    for (i in 1:nrow(pwc.data)){
        if (!(pwc.data[i,4] %in% c(NA, pwc.data[i,2], pwc.data[i,3]))){
            cat("pwc.data[", i, ",4] is not in {NA, pwc.data[", i, ",2:3]}.\n", sep="")
            stop("Please check data and try MCMCpaircompare2dDP() again.\n",
                 call.=FALSE)      
        }
    }
    
    
    
    ## extract key constants from pwc.data
    n <- nrow(pwc.data)
    n.resp <- length(unique(pwc.data[,1]))
    n.cand <- length(unique( c(pwc.data[,2], pwc.data[,3])))
    
    
    ## convert pwc.data into purely numeric matrix
    resp.codes <- sort(unique(pwc.data[,1]))
    cand.codes <- sort(unique( c(pwc.data[,2], pwc.data[,3]) ))
    pwc.data.numeric <- matrix(-999, nrow(pwc.data), 4)
    for (p in 1:n){
        if (pwc.data[p,1] > 0){
            pwc.data.numeric[p,1] <- which(pwc.data[p,1] == resp.codes)
        }
        if (pwc.data[p,2] > 0){
            pwc.data.numeric[p,2] <- which(pwc.data[p,2] == cand.codes)
        }
        if (pwc.data[p,3] > 0){
            pwc.data.numeric[p,3] <- which(pwc.data[p,3] == cand.codes)
        }
        if (pwc.data[p,4] > 0){
            pwc.data.numeric[p,4] <- which(pwc.data[p,4] == cand.codes)
        }
    }
    

    if (length(cluster.max) > 1){
        cat("cluster.max is not a scalar in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (!is.numeric(cluster.max)){
        cat("cluster.max is non-numeric in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (cluster.max <= 0){
        cat("cluster.max takes a value <= 0 in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }

    if (length(cluster.mcmc) > 1){
        cat("cluster.mcmc is not a scalar in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (!is.numeric(cluster.mcmc)){
        cat("cluster.mcmc is non-numeric in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (cluster.mcmc <= 0){
        cat("cluster.mcmc takes a value <= 0 in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }


    if (length(a0) > 1){
        cat("a0 is not a scalar in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (!is.numeric(a0)){
        cat("a0 is non-numeric in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (a0 <= 0){
        cat("a0 takes a value <= 0 in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }


    if (length(b0) > 1){
        cat("b0 is not a scalar in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (!is.numeric(b0)){
        cat("b0 is non-numeric in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (b0 <= 0){
        cat("b0 takes a value <= 0 in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    
    ## set up constraints on theta
    holder <- build.pairwise.theta.constraints(theta.constraints,
                                               cand.codes, n.cand, 2)
    theta.eq.constraints <- holder[[1]]
    theta.ineq.constraints <- holder[[2]]
    
    ## starting values for theta
    theta <- pairwise.theta.start(theta.start, n.cand, 2,
                                  theta.eq.constraints,
                                  theta.ineq.constraints)

    ## starting values for alpha
    if (length(alpha.start) > 1){
        cat("alpha.start is not a scalar in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (!is.numeric(alpha.start)){
        cat("alpha.start is non-numeric in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (alpha.start <= 0){
        cat("alpha.start takes a value <= 0 in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    
    
    ## starting values for cluster.gamma
    cluster.gamma<-gamma.start

    if (all(is.na(gamma.start))){
       cluster.gamma <- gamma.start <- rep(pi/4, cluster.max)
    }
    if (length(gamma.start) < cluster.max){
       cluster.gamma <- rep(gamma.start, length.out=cluster.max)
    }
    if (!is.numeric(gamma.start)){
        cat("gamma.start is non-numeric in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (min(gamma.start) < 0){
        cat("gamma.start takes a value < 0 in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    if (max(gamma.start) > pi/2){
        cat("gamma.start takes a value > pi/2 in MCMCpaircompare2dDP().\n")
        stop("Please check specification and try MCMCpaircompare2dDP() again.\n",
             call.=FALSE)
    }
    
    ##starting values for judge membership
    judge.cluster.membership<-rep(NA, n.resp)
    ##starting values gamma_i, which depends on judge i's cluster membership
    gamma<-rep(NA,n.resp)
    ## I randomly assign judges to 1:cluster.max clusters with equal likelihood
    for(i in 1:n.resp){
      judge.cluster.membership[i]<-sample(1:cluster.max, 1)
      gamma[i]<-cluster.gamma[judge.cluster.membership[i]]
    }

    
    ## define holder for posterior sample
    if(store.gamma == FALSE & store.theta == TRUE) {
        sample <- matrix(data=0, mcmc/thin, 2*n.cand)
    }
    else if (store.gamma == TRUE & store.theta == FALSE){
      #we store judge i's gamma, judge i's cluster membership, and the # of unique clusters in one iteration
        sample <- matrix(data=0, mcmc/thin, n.resp*2+1)
    }
    else if (store.gamma == TRUE & store.theta == TRUE){
        sample <- matrix(data=0, mcmc/thin, 2*n.cand + n.resp*2+1)
    }
    else{
        cat("Error: store.gamma == FALSE & store.theta == FALSE.\n")
        stop("Please respecify and call MCMCpaircompare() again.\n",
             call.=FALSE)      
    }
    #if alpha.fixed=FALSE, we store alpha posterior sample
    if(alpha.fixed==FALSE){
      sample<-cbind(sample,0)
    }
    
    
    ## define holder for the acceptance rate of each gamma    
    gamma_accept_rate<-rep(0, n.resp)
    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    

    ## create theta.sub.target and associated row indicator
    theta.eq.const.holder <- theta.eq.constraints
    theta.ineq.const.holder <- theta.ineq.constraints
    theta.eq.const.holder[theta.eq.const.holder[,1] == -999 &
                          theta.ineq.const.holder[,1] !=  0, 1] <-
        theta.ineq.const.holder[theta.eq.const.holder[,1] == -999 &
                                theta.ineq.const.holder[,1] !=  0, 1]
    theta.eq.const.holder[theta.eq.const.holder[,2] == -999 &
                          theta.ineq.const.holder[,2] !=  0, 2] <-
        theta.ineq.const.holder[theta.eq.const.holder[,2] == -999 &
                                theta.ineq.const.holder[,2] !=  0, 2]

    theta.sub.target.indic <- ((theta.eq.const.holder[,1] != -999) &
        (theta.eq.const.holder[,2] != -999))
    
    theta.sub.target <- theta.eq.const.holder[theta.sub.target.indic,]
        
    if (procrustes == TRUE){
        if (nrow(theta.sub.target) < 3){
        cat("Error: procrustes == TRUE but theta.constraints has < 3 rows.\n")
        stop("Please respecify and call MCMCpaircompare() again.\n",
             call.=FALSE)      
        }
        else{
            theta.eq.constraints <- 0 * theta.eq.constraints - 999
            theta.ineq.constraints <- 0 * theta.ineq.constraints
        }        
    }
    
    
    
    ## call C++ code to draw sample
    posterior <- .C("cMCMCpaircompare2dDP",
                    sampledata = as.double(sample),
                    samplerow = as.integer(nrow(sample)),
                    samplecol = as.integer(ncol(sample)),
                    pwc.datanumericdata = as.integer(pwc.data.numeric-1),#minus one because respondents and items are indexed from zero in C++
                    pwc.datanumericrow = as.integer(nrow(pwc.data.numeric)),
                    pwc.datanumericcol = as.integer(ncol(pwc.data.numeric)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    clustermcmc = as.integer(cluster.mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    thetastartdata = as.double(theta),
                    thetastartrow = as.integer(nrow(theta)),
                    thetastartcol = as.integer(ncol(theta)),
                    gammastartdata = as.double(gamma),
                    gammastartrow = as.integer(length(gamma)),
                    gammastartcol = as.integer(1),
                    clustergammastartdata=as.double(cluster.gamma), 
                    clustergammastartrow= as.integer(length(cluster.gamma)), 
                    clustergammastartcol= as.integer(1),
                    judgeclustermembershipstartdata=as.integer(judge.cluster.membership-1), #minus one because clusters are indexed from zero in C++
                    judgeclustermembershipstartrow=as.integer(length(judge.cluster.membership)), 
                    judgeclustermembershipstartcol= as.integer(1),
                    tunevalue=as.double(tune),
                    thetaeqdata = as.double(theta.eq.constraints),
                    thetaeqrow = as.integer(nrow(theta.eq.constraints)),
                    thetaeqcol = as.integer(ncol(theta.eq.constraints)),
                    thetaineqdata = as.double(theta.ineq.constraints),
                    thetaineqrow = as.integer(nrow(theta.ineq.constraints)),
                    thetaineqcol = as.integer(ncol(theta.ineq.constraints)),
                    storegamma = as.integer(store.gamma),
                    storetheta = as.integer(store.theta),
                    gammaacceptrate=as.double(gamma_accept_rate),
                    alpha=as.double(alpha.start),
                    clustermax=as.integer(cluster.max),
                    alphafixed=as.integer(alpha.fixed), 
                    a=as.double(a0),
                    b=as.double(b0),
                    PACKAGE="MCMCpack"
                    )
    
    ## undo the C++ indexing by 0
    posterior$pwc.datanumericdata <- posterior$pwc.datanumericdata + 1
    
    theta.names <- c(paste("theta1.", cand.codes, sep = ""), paste("theta2.", cand.codes, sep = ""))
    ## I store theta's in the following way:
    ## if we have 1,2,3,4,5 cand.codes,
    ## then theta.names is "theta1.1" "theta1.2" "theta1.3" "theta1.4" "theta1.5"
    ##                     "theta2.1" "theta2.2" "theta2.3" "theta2.4" "theta2.5"
    gamma.names <- paste("gamma.", resp.codes, sep = "")
    
    cluster.names <- paste("cluster.", resp.codes, sep = "")
    
    ## put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol,
                     byrow=FALSE)
    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    
    names <- NULL
    if(store.theta == TRUE) {
        names <- c(names, theta.names)
    }
    if (store.gamma == TRUE){
        names <- c(names, gamma.names, cluster.names, "n.clusters")
    }
    if (alpha.fixed == FALSE){
      names <- c(names, "alpha")
    }
    varnames(output) <- names


    if (procrustes){
        gammas <- output[, grep("gamma", colnames(output))]
        thetas <- output[, grep("theta", colnames(output))]
        thetas1 <- thetas[,1:n.cand]
        thetas2 <- thetas[,(n.cand+1):ncol(thetas)]
        for (iter in 1:nrow(thetas)){
            theta.sub <- cbind(thetas1[iter, theta.sub.target.indic],
                               thetas2[iter,  theta.sub.target.indic])
            procrust.out <- procrustes(theta.sub, theta.sub.target,
                                       translation=FALSE, dilation=FALSE)
            R <- procrust.out$R
            
            theta.mat <- cbind(thetas1[iter,], thetas2[iter,])
            theta.mat <- theta.mat %*% R
            thetas1[iter, ] <- theta.mat[,1]
            thetas2[iter, ] <- theta.mat[,2]
            
            z <- cbind(cos(gammas[iter,]), sin(gammas[iter,]))
            z <- t(R) %*% t(z)
            gammas[iter,] <- atan2(z[2,], z[1,])            
        }
        # output <- mcmc(data=cbind(thetas1, thetas2, gammas),
        #                start=burnin+1, end=burnin+mcmc, thin=thin)
        output[,1:(2*n.cand)]<-cbind(thetas1, thetas2)
        output[,(1+2*n.cand):(n.resp+2*n.cand)]<-gammas
    }


    
    attr(output,"title") <-
        "MCMCpaircompare2dDP Posterior Sample"
    gamma_accept_rate <- posterior$gammaacceptrate
    attr(output, "gamma_accept_rate") <- gamma_accept_rate
    attr(output, "procrustes") <- procrustes
    return(output)
    
    
} ## end MCMCpaircompare





