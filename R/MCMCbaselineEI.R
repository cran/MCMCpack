# sample from the posterior distribution of Wakefield's baseline model
# for ecological inference in R using linked C++ code in Scythe
#
# KQ 10/22/2002

"MCMCbaselineEI" <-
  function(r0, r1, c0, c1, burnin=1000, mcmc=50000, thin=10,
           tune=2.65316, verbose=FALSE, seed=NA,
           a0=1, b0=1, a1=1, b1=1, method="NA",
           ...){


    # Error checking
    if (length(r0) != length(r1)){
      cat("length(r0) != length(r1).\n")
      stop("Please check data and try MCMCbaselineEI() again.\n")
    }

    if (length(r0) != length(c0)){
      cat("length(r0) != length(c0).\n")
      stop("Please check data and try MCMCbaselineEI() again.\n")  
     }

    if (length(r0) != length(c1)){
      cat("length(r0) != length(c1).\n")
      stop("Please check data and try MCMCbaselineEI() again.\n")
    }

    if (length(r1) != length(c0)){
      cat("length(r1) != length(c0).\n")
      stop("Please check data and try MCMCbaselineEI() again.\n")
    }
    if (length(r1) != length(c1)){
      cat("length(r1) != length(c1).\n")
      stop("Please check data and try MCMCbaselineEI() again.\n")
    }
    
    if (length(c0) != length(c1)){
      cat("length(c0) != length(c1).\n")
      stop("Please check data and try MCMCbaselineEI() again.\n")
    }

    if (min((r0+r1) == (c0+c1))==0){
      cat("Rows and columns do not sum to same thing.\n")
      stop("Please check data and try MCMCbaselineEI() again.\n")
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
      stop("Please respecify and try MCMCbaselineEI() again.\n")
    }

    if (b0 <= 0 ){
      cat("Parameter b0 <= 0.\n")
      stop("Please respecify and try MCMCbaselineEI() again.\n")
    }

    if (a1 <= 0 ){
      cat("Parameter a1 <= 0.\n")
      stop("Please respecify and try MCMCbaselineEI() again.\n")
    }

    if (b1 <= 0 ){
      cat("Parameter b1 <= 0.\n")
      stop("Please respecify and try  MCMCbaselineEI() again.\n")
    }

    if (!(method %in% c("DA", "NA")) ){
      cat("Parameter method not DA or NA.\n")
      stop("Please respecify and try MCMCbaselineEI() again.\n")
    }

    # setup matrix to hold output from sampling
    ntables = length(r0)
    sample <- matrix(0, mcmc/thin, ntables*4)

    # call C++ code to do the sampling
    if (method=="NA"){
      C.sample <- .C("baselineNA",
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
    }
    else if (method=="DA"){
      C.sample <- .C("baselineDA",
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
                     a0 = as.double(a0),
                     b0 = as.double(b0),
                     a1 = as.double(a1),
                     b1 = as.double(b1),
                     verbose = as.integer(verbose),
                     lecuyer = as.integer(lecuyer),
                     seedarray = as.integer(seed.arry),
                     lecuyerstream = as.integer(lecuyer.stream),
                     PACKAGE="MCMCpack"
                     )
    }
    
    sample <- matrix(C.sample$samdata, C.sample$samrow, C.sample$samcol,
                     byrow=TRUE)
    if (method=="NA"){
      sample <- sample[,1:(ntables*2)]
      output <- mcmc(data=sample, start=1, end=mcmc, thin=thin)
      p0names <- paste("p0table", 1:ntables, sep="")
      p1names <- paste("p1table", 1:ntables, sep="")
      varnames(output) <- c(p0names, p1names)
      cat(" Overall acceptance rate = ",
          C.sample$accepts / (C.sample$burnin+C.sample$mcmc) / ntables, "\n")
    }
    else{
      output <- mcmc(data=sample, start=1, end=mcmc, thin=thin)
      p0names <- paste("p0table", 1:ntables, sep="")
      p1names <- paste("p1table", 1:ntables, sep="")
      y0names <- paste("y0table", 1:ntables, sep="")
      y1names <- paste("y1table", 1:ntables, sep="")
      varnames(output) <- c(p0names, p1names, y0names, y1names)
    }
    
    attr(output, "title") <- paste("MCMCpack Wakefield's Baseline EI Model Posterior Density Sample, Method =", method)

    return(output)
    
  }
  
