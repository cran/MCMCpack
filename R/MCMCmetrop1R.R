## samples from a user-written posterior code in R using a
## random walk Metropolis algorithm
##
## KQ 6/24/2004
##

"MCMCmetrop1R" <- function(fun, theta.init,
                           burnin=500, mcmc=20000, thin=1,
                           tune=1, verbose=0, seed=NA, logfun=TRUE,
                           optim.trace=0, optim.REPORT=10,
                           optim.maxit=500, ...){
  
  ## error checking here
  check.offset(list(...))
  check.mcmc.parameters(burnin, mcmc, thin)
    
  ## form the tuning vector
  tune <- vector.tune(tune, length(theta.init))
  
  ## form seed
  seeds <- form.seeds(seed) 
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
  
  ## setup the environment so that fun can see the things passed as ...
  ## it should be the case that users can specify arguments with the same
  ## names as variables defined in MCMCmetrop1R without causing problems.
  my.env <- new.env()
  environment(fun) <- my.env
  dots <- list(...)
  dotnames <- names(dots)
  ndots <- length(dots)
  if (ndots >= 1){
    for (i in 1:ndots){
      assign(x=dotnames[[i]], value=dots[[i]], inherits=FALSE, envir=my.env)
    }
  }
  
  ## find approx mode and Hessian using optim()
  opt.out <- optim(theta.init, fun,
                   control=list(fnscale=-1, trace=optim.trace,
                     REPORT=optim.REPORT, maxit=optim.maxit),
                   method="BFGS", hessian=TRUE)
  if(opt.out$convergence!=0){
    warning("Mode and Hessian were not found with call to optim().\nSampling proceeded anyway. \n") 
  }
  
  propvar <- tune %*% solve(-1*opt.out$hessian) %*% tune
  
  ## Call the C++ function to do the MCMC sampling 
  sample <- .Call("MCMCmetrop1R_cc", fun, as.double(theta.init),
                  my.env, as.integer(burnin), as.integer(mcmc),
                  as.integer(thin),
                  as.logical(verbose),
                  lecuyer=as.integer(lecuyer), 
                  seedarray=as.integer(seed.array),
                  lecuyerstream=as.integer(lecuyer.stream),
                  as.logical(logfun),
                  as.matrix(propvar),
                  PACKAGE="MCMCpack")

  ## turn sample into an mcmc object
  sample <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
  return(sample)
}
 
