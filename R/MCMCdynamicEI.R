# sample from the posterior of Quinn's dynamic ecological inference model
# in R using linked C++ code in Scythe
#
# KQ 10/25/2002

"MCMCdynamicEI" <-
  function(r0, r1, c0, c1, burnin=5000, mcmc=200000,
           thin=200, tune=2.65316, verbose=FALSE, seed=NA,
           W=0, a0=1, b0=0.03, a1=1,
           b1=0.03, ...){
    
    # Error checking
    if (length(r0) != length(r1)){
      cat("length(r0) != length(r1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r0) != length(c0)){
      cat("length(r0) != length(c0).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r0) != length(c1)){
      cat("length(r0) != length(c1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r1) != length(c0)){
      cat("length(r1) != length(c0).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r1) != length(c1)){
      cat("length(r1) != length(c1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(c0) != length(c1)){
      cat("length(c0) != length(c1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (min((r0+r1) == (c0+c1))==0){
      cat("Rows and columns do not sum to same thing.\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    check.mcmc.parameters(burnin, mcmc, thin)
    tune <- scalar.tune(tune)

    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    if (a0 <= 0 ){
      cat("Parameter a0 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }
    
    if (b0 <= 0 ){
      cat("Parameter b0 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }
    
    if (a1 <= 0 ){
      cat("Parameter a1 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }
    
    if (b1 <= 0 ){
      cat("Parameter b1 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }
    
    ntables = length(r0)
    
    if (W==0){ # construct weight matrix for a simple random walk assuming
               # tables are temporally ordered and 1 time unit apart
      W <- matrix(0, ntables, ntables)
      for (i in 2:(ntables)){
        W[i,i-1] <- 1
        W[i-1,i] <- 1
      }
    }

    # setup matrix to hold output from sampling
    sample <- matrix(0, mcmc/thin, ntables*2+2)

    # call C++ code to do the sampling
    C.sample <- .C("dynamicEI",
                   samdata = as.double(sample),
                   samrow = as.integer(nrow(sample)),
                   samcol = as.integer(ncol(sample)),
                   r0 = as.double(r0),
                   r1 = as.double(r1),
                   c0 = as.double(c0),
                   c1 = as.double(c1),
                   ntables = as.integer(ntables),
                   burnin = as.integer(burnin),
                   mcmc = as.integer(mcmc),
                   thin = as.integer(thin),
                   W = as.double(W),
                   a0 = as.double(a0),
                   b0 = as.double(b0),
                   a1 = as.double(a1),
                   b1 = as.double(b1),
                   verbose = as.integer(verbose),
                   tune = as.double(tune),
                   lecuyer = as.integer(lecuyer),
                   seedarray = as.integer(seed.array),
                   lecuyerstream = as.integer(lecuyer.stream),
                   accepts = as.integer(0),
                   PACKAGE="MCMCpack"
                   )

    cat(" Overall acceptance rate = ",
        C.sample$accepts / (C.sample$burnin+C.sample$mcmc) / ntables, "\n")
    
    sample <- matrix(C.sample$samdata, C.sample$samrow, C.sample$samcol,
                     byrow=TRUE)
    output <- mcmc(data=sample, start=1, end=mcmc, thin=thin)
    p0names <- paste("p0table", 1:ntables, sep="")
    p1names <- paste("p1table", 1:ntables, sep="")
    varnames(output) <- c(p0names, p1names, "sigma^2_0", "sigma^2_1")
    
    attr(output, "title") <- "MCMCpack Quinn's Dynamic EI Model Posterior Density Sample" 
    
    
    return(output)
    
  }
