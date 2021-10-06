##########################################################################
## sample from the posterior distribution of a pairwise comparisons model
## in R using linked C++ code in Scythe.
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
## Pr(Y_{ijj'} = 1) = \Phi( \alpha_{i} [\theta_{j} - \theta_{ j'} ] )
##
## alpha_i \overset{iid}{\sim} N(a, A^{-1})
## theta_j \overset{ind}{\sim} N(0, 1)
##            (some theta_js truncated above or below 0, or fixed to constants)
##
##
## candidate IDs in columns 2 to 4 need to begin with a letter
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Original code KQ 3/17/2015
## Added to MCMCpack KQ 6/24/2021
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for a Pairwise Comparisons Model with Probit Link
#'
#' This function generates a sample from the posterior distribution of a
#' model for pairwise comparisons data with a probit link. Thurstone's model
#' is a special case of this model when the \eqn{\alpha} parameter is fixed at
#' 1.
#'
#' \code{MCMCpaircompare} uses the data augmentation approach of Albert and
#' Chib (1993). The user supplies data and priors, and a sample from the
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
#' \deqn{Pr(Y_{ijj'} = 1) = \Phi( \alpha_{i} [\theta_{j} - \theta_{ j'} ] ) }
#'
#' The following Gaussian priors are assumed:
#' \deqn{\alpha_i \sim \mathcal{N}(a, A^{-1})}
#' \deqn{\theta_j \sim \mathcal{N}(0, 1)}
#' For identification, some \eqn{\theta_j}s are truncated above or below 0,
#' or fixed to constants.
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
#' For applications without raters (such as sports competitions) all entries
#' in the first column should be set to a single value and \code{alpha.fixed}
#' (see below) should be set to \code{TRUE}. \strong{The identifiers in
#' columns 2 through 4 must start with a letter. Examples are provided below.} 
#'
#' @param theta.constraints A list specifying possible simple equality or
#' inequality constraints on the item parameters. A typical entry in the
#' list has one of three forms: \code{itemname=c} which will constrain the
#' item parameter for the item named \code{itemname} to be equal to c,
#' \code{itemname="+"} which will constrain the item parameter for the
#' item named \code{itemname} to be positive, and \code{itemname="-"} which
#' will constrain the item parameter for the item named \code{itemname} to
#' be negative. 
#'
#' @param alpha.fixed Should alpha be fixed to a constant value of 1 for all
#' raters? Default is FALSE. If set to FALSE, an alpha value is estimated for
#' each rater. 
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
#' @param alpha.start The starting value for the alpha vector.  This
#' can either be a scalar or a column vector with dimension equal to the number
#' of alphas.  If this takes a scalar value, then that value will serve as the
#' starting value for all of the alphas. The default value of NA will set the
#' starting value of each alpha parameter to 1.
#'
#' @param a The prior mean of alpha. Must be a scalar. Default is 0. 
#'
#' @param A The prior precision of alpha. Must be a positive scalar.
#' Default is 0.25 (prior variance is 4). 
#'
#' @param store.theta Should the theta draws be returned? Default is TRUE. 
#'
#' @param store.alpha Should the alpha draws be returned? Default is FALSE.
#'
#' @param ... further arguments to be passed
#'
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MCMCpack]{MCMCpaircompare2d}},
#' \code{\link[MCMCpack]{MCMCpaircompare2dDP}}
#'
#' @references Albert, J. H. and S. Chib. 1993. ``Bayesian Analysis of Binary
#' and Polychotomous Response Data.'' \emph{J. Amer. Statist. Assoc.} 88,
#' 669-679
#'
#' Yu, Qiushi and Kevin M. Quinn. 2021. ``A Multidimensional Pairwise
#' Comparison Model for Heterogeneous Perception with an Application to
#' Modeling the Perceived Truthfulness of Public Statements on COVID-19.''
#' University of Michigan Working Paper.
#' 
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21. \doi{10.18637/jss.v042.i09}.
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
#'   \dontrun{
#'   ## Euro 2016 example
#'   data(Euro2016)
#'
#' posterior1 <- MCMCpaircompare(pwc.data=Euro2016,
#'                               theta.constraints=list(Ukraine="-",
#'                                                      Portugal="+"),
#'                               alpha.fixed=TRUE,
#'                               verbose=10000,
#'                               burnin=10000, mcmc=500000, thin=100,
#'                               store.theta=TRUE, store.alpha=FALSE)
#'
#' ## alternative identification constraints
#' posterior2 <- MCMCpaircompare(pwc.data=Euro2016,
#'                               theta.constraints=list(Ukraine="-",
#'                                                      Portugal=1),
#'                               alpha.fixed=TRUE,
#'                               verbose=10000,
#'                               burnin=10000, mcmc=500000, thin=100,
#'                               store.theta=TRUE, store.alpha=FALSE)
#'
#'
#'
#'
#'
#'
#'
#'
#' ## a synthetic data example with estimated rater-specific parameters
#' set.seed(123)
#' 
#' I <- 65  ## number of raters
#' J <- 50 ## number of items to be compared
#' 
#' 
#' ## raters 1 to 5 have less sensitivity to stimuli than raters 6 through I
#' alpha.true <- c(rnorm(5, m=0.2, s=0.05), rnorm(I - 5, m=1, s=0.1))
#' theta.true <- sort(rnorm(J, m=0, s=1))
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
#'         eta <- alpha.true[i] * (theta.true[item.1] - theta.true[item.2])
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
#' ## fit the model
#' posterior <- MCMCpaircompare(pwc.data=sim.data,
#'                              theta.constraints=list(item.101=-2,
#'                                                     item.150=2),
#'                              alpha.fixed=FALSE,
#'                              verbose=10000,
#'                              a=0, A=0.5,
#'                              burnin=10000, mcmc=200000, thin=100,
#'                              store.theta=TRUE, store.alpha=TRUE)
#' 
#' theta.draws <- posterior[, grep("theta", colnames(posterior))]
#' alpha.draws <- posterior[, grep("alpha", colnames(posterior))]
#' 
#' theta.post.med <- apply(theta.draws, 2, median)
#' alpha.post.med <- apply(alpha.draws, 2, median)
#' 
#' theta.post.025 <- apply(theta.draws, 2, quantile, prob=0.025)
#' theta.post.975 <- apply(theta.draws, 2, quantile, prob=0.975)
#' alpha.post.025 <- apply(alpha.draws, 2, quantile, prob=0.025)
#' alpha.post.975 <- apply(alpha.draws, 2, quantile, prob=0.975)
#' 
#' ## compare estimates to truth
#' par(mfrow=c(1,2))
#' plot(theta.true, theta.post.med, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5),
#'      col=rgb(0,0,0,0.3))
#' segments(x0=theta.true, x1=theta.true,
#'          y0=theta.post.025, y1=theta.post.975,
#'          col=rgb(0,0,0,0.3)) 
#' abline(0, 1, col=rgb(1,0,0,0.5))
#' 
#' plot(alpha.true, alpha.post.med, xlim=c(0, 1.2), ylim=c(0, 3),
#'      col=rgb(0,0,0,0.3))
#' segments(x0=alpha.true, x1=alpha.true,
#'          y0=alpha.post.025, y1=alpha.post.975,
#'          col=rgb(0,0,0,0.3)) 
#' abline(0, 1, col=rgb(1,0,0,0.5))
#' 
#' }
#' 
#' @export


"MCMCpaircompare" <- function(pwc.data, theta.constraints=list(),
                              alpha.fixed=FALSE,
                              burnin=1000, mcmc=20000, thin=1,
                              verbose=0, seed=NA,
                              alpha.start=NA,
                              a=0, A=0.25,
                              store.theta=TRUE,
                              store.alpha=FALSE,
                              ...){

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    if (!is.logical(alpha.fixed)){
        cat("alpha.fixed must be a logical value.\n")
        stop("Please check data and try MCMCpaircompare() again.\n",
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
        stop("Please check data and try MCMCpaircompare() again.\n",
             call.=FALSE)
    }
    for (i in 1:nrow(pwc.data)){
        if (!(pwc.data[i,4] %in% c(NA, pwc.data[i,2], pwc.data[i,3]))){
            cat("pwc.data[", i, ",4] is not in {NA, pwc.data[", i, ",2:3]}.\n", sep="")
            stop("Please check data and try MCMCpaircompare() again.\n",
                 call.=FALSE)      
        }
    }
    if (!is.numeric(a) | length(a) !=1){
        cat("a must be a scalar.\n")
        stop("Please check specification and try MCMCpaircompare() again.\n",
             call.=FALSE)          
    }
    if (!is.numeric(A) | length(A) !=1 | A <= 0){
        cat("A must be a positive scalar.\n")
        stop("Please check specification and try MCMCpaircompare() again.\n",
             call.=FALSE)          
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
        if (!is.na(pwc.data[p,1])){
            pwc.data.numeric[p,1] <- which(pwc.data[p,1] == resp.codes)
        }
        if (!is.na(pwc.data[p,2])){
            pwc.data.numeric[p,2] <- which(pwc.data[p,2] == cand.codes)
        }
        if (!is.na(pwc.data[p,3])){
            pwc.data.numeric[p,3] <- which(pwc.data[p,3] == cand.codes)
        }
        if (!is.na(pwc.data[p,4])){
            pwc.data.numeric[p,4] <- which(pwc.data[p,4] == cand.codes)
        }
    }
    
    
    ## set up constraints on theta
    if(length(theta.constraints) != 0) {
        for (i in 1:length(theta.constraints)){
            theta.constraints[[i]] <-
                list(as.integer(1), theta.constraints[[i]][1])
        }
    }
    theta.eq.constraints <- matrix(NA, n.cand)
    theta.ineq.constraints <- matrix(0, n.cand)
    rownames(theta.eq.constraints) <- cand.codes
    rownames(theta.ineq.constraints) <- cand.codes
    if (length(theta.constraints) != 0){
        constraint.names <- names(theta.constraints)  
        for (i in 1:length(constraint.names)){
            name.i <- constraint.names[i]
            theta.constraints.i <- theta.constraints[[i]]
            cand.index <- theta.constraints.i[[1]]
            replace.element <- theta.constraints.i[[2]]
            if (is.numeric(replace.element)){
                theta.eq.constraints[rownames(theta.eq.constraints)==name.i,
                                     cand.index] <- replace.element
            }
            if (replace.element=="+"){
                theta.ineq.constraints[rownames(theta.ineq.constraints)==name.i,
                                       cand.index] <- 1 
            }
            if (replace.element=="-"){
                theta.ineq.constraints[rownames(theta.ineq.constraints)==name.i,
                                       cand.index] <- -1
            }
        }
    }
    
    testmat <- theta.ineq.constraints * theta.eq.constraints
    
    if (min(is.na(testmat))==0){
        if ( min(testmat[!is.na(testmat)]) < 0){
            cat("Constraints on theta are logically inconsistent.\n")
            stop("Please respecify and call ", calling.function(), " again.\n")
        }
    }
    theta.eq.constraints[is.na(theta.eq.constraints)] <- -999
    
    
    ## starting values for theta
    theta.start <- rep(0, n.cand)
    
    for (j in 1:n.cand){
        cand.code.j <- cand.codes[j]
        if (theta.eq.constraints[cand.code.j,1] != -999){
            theta.start[j] <- theta.eq.constraints[cand.code.j,1]
        }
        if (theta.ineq.constraints[cand.code.j,1] != 0){
            theta.start[j] <- theta.ineq.constraints[cand.code.j,1]
        }    
    }
    
    
    ## starting values for alpha
    if (is.na(alpha.start)){
        alpha.start <- rep(1, n.resp)
    }
    if (length(alpha.start) < n.resp){
        alpha.start <- rep(alpha.start, length.out=n.resp)
    }
    if (!is.numeric(alpha.start)){
        cat("alpha.start is non-numeric in MCMCpaircompare().\n")
        stop("Please check specification and try MCMCpaircompare() again.\n",
             call.=FALSE)          
    }
    
    if (alpha.fixed){
        alpha.start <- rep(1, n.resp)
    }
    
    
    ## define holder for posterior sample
    if(store.alpha == FALSE & store.theta == TRUE) {
        sample <- matrix(data=0, mcmc/thin, n.cand)
    }
    else if (store.alpha == TRUE & store.theta == FALSE){
        sample <- matrix(data=0, mcmc/thin, n.resp)
    }
    else if (store.alpha == TRUE & store.theta == TRUE){
        sample <- matrix(data=0, mcmc/thin, n.cand + n.resp)
    }
    else{
        cat("Error: store.alpha == FALSE & store.theta == FALSE.\n")
        stop("Please respecify and call MCMCpaircompare() again.\n",
             call.=FALSE)      
    }
    
    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    
    
    ## call C++ code to draw sample
    posterior <- .C("cMCMCpaircompare",
                    sampledata = as.double(sample),
                    samplerow = as.integer(nrow(sample)),
                    samplecol = as.integer(ncol(sample)),
                    pwc.datanumericdata = as.integer(pwc.data.numeric-1),
                    pwc.datanumericrow = as.integer(nrow(pwc.data.numeric)),
                    pwc.datanumericcol = as.integer(ncol(pwc.data.numeric)),
                    alphafixed = as.integer(alpha.fixed),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    thetastartdata = as.double(theta.start),
                    thetastartrow = as.integer(length(theta.start)),
                    thetastartcol = as.integer(1),
                    astartdata = as.double(alpha.start),
                    astartrow = as.integer(length(alpha.start)),
                    astartcol = as.integer(1),
                    a=as.double(a),
                    A=as.double(A),
                    thetaeqdata = as.double(theta.eq.constraints),
                    thetaeqrow = as.integer(nrow(theta.eq.constraints)),
                    thetaeqcol = as.integer(ncol(theta.eq.constraints)),
                    thetaineqdata = as.double(theta.ineq.constraints),
                    thetaineqrow = as.integer(nrow(theta.ineq.constraints)),
                    thetaineqcol = as.integer(ncol(theta.ineq.constraints)),
                    storealpha = as.integer(store.alpha),
                    storetheta = as.integer(store.theta),
                    PACKAGE="MCMCpack"
                    )
    
    ## undo the C++ indexing by 0
    posterior$pwc.datanumericdata <- posterior$pwc.datanumericdata + 1
    
    theta.names <- paste("theta.", cand.codes, sep = "")
    alpha.names <- paste("alpha.", resp.codes, sep = "")
    
    ## put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol,
                     byrow=FALSE)
    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    
    names <- NULL
    if(store.theta == TRUE) {
        names <- c(names, theta.names)
    }
    if (store.alpha == TRUE){
        names <- c(names, alpha.names)
    }
    varnames(output) <- names
    attr(output,"title") <-
        "MCMCpaircompare Posterior Sample"
    return(output)
    
    
} ## end MCMCpaircompare





