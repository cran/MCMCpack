## sample from the posterior distribution of a logistic regression
## model in R using linked C++ code in Scythe
##
## KQ 1/23/2003
##
## Modified to meet new developer specification 7/15/2004 KQ
## Modified for new Scythe and rngs 7/25/2004 KQ
## note: B0 is now a precision

"MCMClogit" <-
  function(formula, data = parent.frame(), burnin = 1000, mcmc = 10000,
           thin=1, tune=1.1, verbose = FALSE, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, ...) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]    
    K <- ncol(X)  # number of covariates
        
    ## starting values and priors
    beta.start <- coef.start(beta.start, K, formula, family=binomial, data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    ## form the tuning parameter
    tune <- vector.tune(tune, K)
    V <- vcov(glm(formula=formula, data=data, family=binomial))
      
    ## y \in {0, 1} error checking
    if (sum(Y!=0 & Y!=1) > 0) {
      cat("Elements of Y equal to something other than 0 or 1.\n")
      stop("Check data and call MCMClogit() again. \n") 
    }
   
   
    ## define holder for posterior density sample
    sample <- matrix(data=0, mcmc/thin, dim(X)[2] )

    ## call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMClogit",
                     sample=sample, Y=Y, X=X, burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     tune=tune, lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, V=V) 
    
    ## put together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior, names=xnames,
                               title="MCMClogit Posterior Density Sample")
    return(output)    
  }

##########################################################################


