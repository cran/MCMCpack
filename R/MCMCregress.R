# MCMCregress.R samples from the posterior distribution of a Gaussian
# linear regression model in R using linked C++ code in Scythe
#
# Original written by ADM and KQ 5/21/2002
# Updated with helper functions ADM 5/28/2004
# Modified to meet new developer specification 6/18/2004 KQ
# Modified for new Scythe and rngs 7/22/2004 ADM

"MCMCregress" <-
  function(formula, data=parent.frame(), burnin = 1000, mcmc = 10000,
   thin=1, verbose = FALSE, seed = NA, beta.start = NA,
   b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, ...) {
    
    # checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    
    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    # form response and model matrices
    holder <- parse.formula(formula, data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]    
    K <- ncol(X)  # number of covariates
        
    # starting values and priors
    beta.start <- coef.start(beta.start, K, formula, family=gaussian, data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    check.ig.prior(c0, d0)
   
    # define holder for posterior density sample
    sample <- matrix(data=0, mcmc/thin, K+1)

    # call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCregress", 
                     sample=sample, Y=Y, X=X, burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, c0=as.double(c0), d0=as.double(d0))
                  
    # pull together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior,
                               names=c(xnames, "sigma2"),
                               title="MCMCregress Posterior Density Sample")
    return(output)
  }
