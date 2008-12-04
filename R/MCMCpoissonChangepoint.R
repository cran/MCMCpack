##########################################################################
## sample from the posterior distribution
## of a Poisson model with multiple changepoints
## using linked C++ code in Scythe 1.0
##
## JHP 07/01/2007
##
## Revised on 09/12/2007 JHP	  
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

"MCMCpoissonChangepoint"<-
    function(data,  m = 1, c0 = NA, d0 = NA, a = NULL, b = NULL,
            burnin = 10000, mcmc = 10000, thin = 1, verbose = 0,
            seed = NA, lambda.start = NA, P.start = NA,
            marginal.likelihood = c("none", "Chib95"), ...) {

    ## check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    cl <- match.call()

    ## ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
   
    ## sample size
    y <- as.matrix(data)
    n <- nrow(y)
    ns <- m+1

    ## prior 
    A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)
    if (is.na(c0)||is.na(d0))
        stop("Please specify prior for lambda (c0 and d0) and call MCMCpoissonChangepoint again.\n")
    
    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)

    ## following MCMCregress, set chib as binary
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    Pstart <- check.P(P.start, m=m, n=n, a=a, b=b)
    lambdastart <- check.theta(lambda.start, ns, y, min=range(y)[1], max=range(y)[2])

    nstore <- mcmc/thin

    ## call C++ code to draw sample
    posterior <- .C("MCMCpoissonChangepoint",
                    lambdaout = as.double(rep(0.0, nstore*ns)),
                    Pout = as.double(rep(0.0, nstore*ns*ns)),
                    psout = as.double(rep(0.0, n*ns)),
                    sout = as.double(rep(0.0, nstore*n)),
                    Ydata = as.double(y),
                    Yrow = as.integer(nrow(y)),
                    Ycol = as.integer(ncol(y)),
                    m = as.integer(m),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
					lecuyer=as.integer(lecuyer), 
					seedarray=as.integer(seed.array),
					lecuyerstream=as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    lambdastart = as.double(lambdastart),
                    Pstart = as.double(Pstart),
                    a = as.double(a),
                    b = as.double(b),
                    c0 = as.double(c0),
                    d0 = as.double(d0),
                    A0data = as.double(A0),
                    logmarglikeholder = as.double(0.0),
                    chib = as.integer(chib))

    ## get marginal likelihood if Chib95
    if (marginal.likelihood == "Chib95"){
      logmarglike <- posterior$logmarglikeholder
	  ##loglike <- posterior$loglikeholder
    }

    ## pull together matrix and build MCMC object to return
    lambda.holder <- matrix(posterior$lambdaout, mcmc/thin)
    P.holder    <- matrix(posterior$Pout, mcmc/thin)
    s.holder    <- matrix(posterior$sout, mcmc/thin)
    ps.holder   <- matrix(posterior$psout, n)

    output <- mcmc(data=lambda.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
    varnames(output)  <- paste("lambda.", 1:ns, sep = "")
    attr(output,"title") <- "MCMCpoissonChangepoint Posterior Sample"
    attr(output, "y")    <- y
    attr(output, "m")    <- m
    attr(output, "call") <- cl
    attr(output, "logmarglike") <- logmarglike
    attr(output, "prob.state") <- ps.holder/(mcmc/thin)
    attr(output, "s.store") <- s.holder
    return(output)

 }## end of MCMC function

