# sample from the posterior distribution of Wakefield's baseline model
# for ecological inference in R using linked C++ code in Scythe
#
# KQ 10/22/2002

"MCMChierEI" <-
  function(r0, r1, c0, c1, burnin=1000, mcmc=50000, thin=1,
           m0=0, M0=10,
           m1=0, M1=10,
           nu0=1.0, delta0=0.5,
           nu1=1.0, delta1=0.5,
           verbose=FALSE, tune=2.65316, seed=0, ...){
    
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

    check.parameters(burnin, mcmc, thin, "MCMChierEI", tune)
    
    if (M0 <= 0 ){
      cat("Parameter M0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (M1 <= 0 ){
      cat("Parameter M1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (nu0 <= 0 ){
      cat("Parameter nu0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (nu1 <= 0 ){
      cat("Parameter nu1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (delta0 <= 0 ){
      cat("Parameter delta0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }
    
    if (delta1 <= 0 ){
      cat("Parameter delta1 <= 0.\n")
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
                   nu0 = as.double(nu0),
                   delta0 = as.double(delta0),
                   nu1 = as.double(nu1),
                   delta1 = as.double(delta1),
                   verbose = as.integer(verbose),
                   tune = as.double(tune),
                   seed = as.integer(seed),
                   accepts = as.integer(0),
                   PACKAGE="MCMCpack"
                   )

    cat(" Overall acceptance rate = ",
        C.sample$accepts / (C.sample$burnin+C.sample$mcmc) / ntables, "\n")
    
    sample <- matrix(C.sample$samdata, C.sample$samrow, C.sample$samcol,
                     byrow=TRUE)
    
    output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
    p0names <- paste("p0table", 1:ntables, sep="")
    p1names <- paste("p1table", 1:ntables, sep="")
    varnames(output) <- c(p0names, p1names, "mu0", "mu1", "sigma^2_0",
                          "sigma^2_1")
    
    attr(output, "title") <- "MCMCpack Wakefield's Hierarchical EI Model Posterior Density Sample"
        
    return(output)
    
  }
