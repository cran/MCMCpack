# sample from the posterior distribution of Wakefield's baseline model
# for ecological inference in R using linked C++ code in Scythe
#
# KQ 10/22/2002

"MCMChierEI" <-
  function(r0, r1, c0, c1, burnin=5000, mcmc=50000, thin=1,
           verbose=0, seed=NA,
           m0=0, M0=2.287656,
           m1=0, M1=2.287656,
           a0=0.825, b0=0.0105,
           a1=0.825, b1=0.0105, ...){
    
    # Error checking
    if (length(r0) != length(r1)){
      cat("length(r0) != length(r1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    if (length(r0) != length(c0)){
      cat("length(r0) != length(c0).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }
    
    if (length(r0) != length(c1)){
      cat("length(r0) != length(c1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }
    
    if (length(r1) != length(c0)){
      cat("length(r1) != length(c0).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }
    
    if (length(r1) != length(c1)){
      cat("length(r1) != length(c1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }
    
    if (length(c0) != length(c1)){
      cat("length(c0) != length(c1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }
    
    if (min((r0+r1) == (c0+c1))==0){
      cat("Rows and columns do not sum to same thing.\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    check.mcmc.parameters(burnin, mcmc, thin)

    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
   
    
    if (M0 <= 0 ){
      cat("Parameter M0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (M1 <= 0 ){
      cat("Parameter M1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (a0 <= 0 ){
      cat("Parameter a0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (a1 <= 0 ){
      cat("Parameter a1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (b0 <= 0 ){
      cat("Parameter b0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (b1 <= 0 ){
      cat("Parameter b1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
   
    # setup matrix to hold output from sampling
    ntables = length(r0)
    sample <- matrix(0, mcmc/thin, ntables*2+4)

    # call C++ code to do the sampling
    C.sample <- .C("hierEI",
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
                   mu0priormean = as.double(m0),
                   mu0priorvar = as.double(M0),
                   mu1priormean = as.double(m1),
                   mu1priorvar = as.double(M1),
                   a0 = as.double(a0),
                   b0 = as.double(b0),
                   a1 = as.double(a1),
                   b1 = as.double(b1),
                   verbose = as.integer(verbose),
                   lecuyer = as.integer(lecuyer),
                   seedarray = as.integer(seed.array),
                   lecuyerstream = as.integer(lecuyer.stream),
                   PACKAGE="MCMCpack"
                   )

    sample <- matrix(C.sample$samdata, C.sample$samrow, C.sample$samcol,
                     byrow=TRUE)
    
    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    p0names <- paste("p0table", 1:ntables, sep="")
    p1names <- paste("p1table", 1:ntables, sep="")
    varnames(output) <- c(p0names, p1names, "mu0", "mu1", "sigma^2.0",
                          "sigma^2.1")
    
    attr(output, "title") <- "MCMCpack Wakefield's Hierarchical EI Model Posterior Density Sample"
        
    return(output)
    
  }
