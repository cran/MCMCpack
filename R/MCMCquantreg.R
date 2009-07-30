"MCMCquantreg" <-
  function(formula, p=0.5, data=NULL, burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001,
           ...) {
    
    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    cl <- match.call()
    if (p<=0 || p>=1){
	stop("p must be in (0,1).\n Please respecify and call again.\n")
	}
    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]    
    K <- ncol(X)  # number of covariates
        
    ## starting values and priors
    ols.fit <- lm(formula)
    defaults <- matrix(coef(ols.fit),K,1)
    defaults[1] <- defaults[1]+summary(ols.fit)$sigma*qnorm(p)
    beta.start <- coef.start(beta.start, K, formula, family=gaussian, data, defaults=defaults)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    check.ig.prior(c0, d0)


    B0.eigenvalues <- eigen(B0)$values
    if (min(B0.eigenvalues) < 0){
      stop("B0 is not positive semi-definite.\nPlease respecify and call again.\n")
    }
    

    
    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, K+1)
    posterior <- NULL 
    
    if (p==0.5) {
     ## call C++/Scythe function "MCMCmedreg" to draw samples
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCmedreg", 
                     sample.nonconst=sample, Y=Y, X=X,                      burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, c0=as.double(c0), d0=as.double(d0), package="MCMCpack")
}

    else{    

    ## call C++/Scythe function "MCMCquantreg" to draw samples
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCquantreg", 
                     sample.nonconst=sample, p=as.double(p), Y=Y, X=X, 
                     burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, c0=as.double(c0), d0=as.double(d0), package="MCMCpack")
}
    
## pull together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior,
                               names=c(xnames,"sigma"),
                               title="MCMCquantreg Posterior Sample",
                               y=Y, call=cl)
    
    return(output)
}
