# sample from the posterior distribution of a one-dimensional item
# response theory model in R using linked C++ code in Scythe.
#
# ADM and KQ 1/23/2003
 
"MCMCirt1d" <-
  function(datamatrix, theta.fixed = 1, burnin = 500, mcmc = 1000,
           thin=5, verbose = FALSE, seed = 0, theta.start = NA, 
           alpha.start = NA, beta.start = NA, t0 = 0, T0 = 1,
           b0.alpha = 0, b0.beta = 0, B0.alpha = 1, B0.beta = 1,
           B0.corr = 0, store.item = FALSE, ... ) {

    set.seed(83829)
    
    # burnin / mcmc / thin error checking
    check.parameters(burnin, mcmc, thin, "MCMCirt1d")

    # check vote matrix and convert to work with C++ code 
    datamatrix <- as.matrix(datamatrix)   
    K <- nrow(datamatrix)   # cases
    J <- ncol(datamatrix)   # justices
    if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) != (J*K)) {
      cat("Data matrix contains elements other than 0, 1 or NA.\n")
      stop("Please check data and try MCMCirt1d() again.\n")
    }
    datamatrix[is.na(datamatrix)] <- 9   

    # starting values for theta error checking
    if (is.na(theta.start)) {
      theta.start <- as.matrix(rnorm(J,1))
    }
    else if(is.null(dim(theta.start))) {
      theta.start <- theta.start * matrix(1,J,1)  
    }
    else if((dim(theta.start)[1] != J) || (dim(theta.start)[2] != 1)) {
      cat("Starting value for theta not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }

    # starting values for alpha and beta error checking
    if (is.na(alpha.start)) {
      alpha.start <- as.matrix(rnorm(K,1))
    }
    else if(is.null(dim(alpha.start))) {
      alpha.start <- alpha.start * matrix(1,K,1)  
    }
    else if((dim(alpha.start)[1] != K) || (dim(alpha.start)[2] != 1)) {
      cat("Starting value for alpha not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }      
    if (is.na(beta.start)) {
      beta.start <- as.matrix(rnorm(K,1))
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,K,1)  
    }
    else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
      cat("Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }    

    # theta fixed error checking / this is the index of the person
    # to constrain to be negative, which identifies the model
    theta.fixed <- as.integer(theta.fixed)
    if(theta.fixed < 1 || theta.fixed > J) {
      cat("Index of actor to fix is outside range.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")      
    }

    # prior for theta error checking (mean)
    if(is.null(dim(t0))) {
      t0 <- t0 * matrix(1,J,1)  
    }
    if((dim(t0)[1] != J) || (dim(t0)[2] != 1)) {
      cat("Vector of prior means t0 not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }
   
    # prior for theta error checking (variance)
    if(is.null(dim(T0))) {
      T0 <- T0 * matrix(1,J,1)  
    }
    if((dim(T0)[1] != J) || (dim(T0)[2] != 1)) {
      cat("Vector of prior variances T0 not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }
    if(sum(T0 > 0) != J) {
      cat("Some Elements of Vector of prior variances T0 not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")         
    }

    # prior for alpha and beta error checking (mean)
    if(is.null(dim(b0.alpha))) {
      b0.alpha <- b0.alpha * matrix(1,K,1)  
    }
    if((dim(b0.alpha)[1] != K) || (dim(b0.alpha)[2] != 1)) {
      cat("Vector of prior means b0.alpha not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }
    if(is.null(dim(b0.beta))) {
      b0.beta <- b0.beta * matrix(1,K,1)  
    }
    if((dim(b0.beta)[1] != K) || (dim(b0.beta)[2] != 1)) {
      cat("Vector of prior means b0.beta not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }
         
    # prior for alpha and beta error checking (variance)
    if(is.null(dim(B0.alpha))) {
      B0.alpha <- B0.alpha * matrix(1,K,1)  
    }
    if((dim(B0.alpha)[1] != K) || (dim(B0.alpha)[2] != 1)) {
      cat("Vector of prior variances B0.alpha not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }
    if(sum(B0.alpha > 0) != K) {
      cat("Elements of prior variances B0.alpha not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")         
    }
    if(is.null(dim(B0.beta))) {
      B0.beta <- B0.beta * matrix(1,K,1)  
    }
    if((dim(B0.beta)[1] != K) || (dim(B0.beta)[2] != 1)) {
      cat("Elements of prior variances B0.beta not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }
    if(sum(B0.beta > 0) != K) {
      cat("Elements of prior variances B0.beta negative.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")         
    }
   
    # prior for alpha and beta error checking (correlation)
    if(is.null(dim(B0.corr))) {
      B0.corr <- B0.corr * matrix(1,K,1)  
    }
    if((dim(B0.corr)[1] != K) || (dim(B0.corr)[2] != 1)) {
      cat("Vector of covariances alpha.b0.beta not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")
    }
    if(sum(B0.corr >= -1.0 & B0.corr<=1.0) != K) {
      cat("Elements of prior correlations B0.corr not defined.\n")
      stop("Please respecify and call MCMCirt1d() again.\n")       
    }
    alpha.beta.cov <- B0.corr / (sqrt(B0.alpha) * sqrt(B0.beta))   
         
    # define holder for posterior density sample
    if(store.item == FALSE) {
      sample <- matrix(data=0, mcmc/thin, J)
    }
    else {
      sample <- matrix(data=0, mcmc/thin, J + 2 * K)
    }
   
    # call C++ code to draw sample
    inv.obj <- .C("irt1dpost",
                  samdata = as.double(sample),
                  samrow = as.integer(nrow(sample)),
                  samcol = as.integer(ncol(sample)),
                  votedata = as.double(datamatrix),
                  voterow = as.integer(nrow(datamatrix)),
                  votecol = as.integer(ncol(datamatrix)),    
                  burnin = as.integer(burnin),
                  gibbs = as.integer(mcmc),
                  thin = as.integer(thin),
                  seed = as.integer(seed),
                  verbose = as.integer(verbose),
                  tstartdata = as.double(theta.start),
                  tstartrow = as.integer(nrow(theta.start)),
                  tstartcol = as.integer(ncol(theta.start)),
                  thetafixed = as.integer(theta.fixed),
                  astartdata = as.double(alpha.start),
                  astartrow = as.integer(nrow(alpha.start)),
                  astartcol = as.integer(ncol(alpha.start)),
                  bstartdata = as.double(beta.start),
                  bstartrow = as.integer(nrow(beta.start)),
                  bstartcol = as.integer(ncol(beta.start)),
                  t0data = as.double(t0),
                  t0row = as.integer(nrow(t0)),
                  t0col = as.integer(ncol(t0)),   
                  T0data = as.double(T0),
                  T0row = as.integer(nrow(T0)),
                  T0col = as.integer(ncol(T0)),
                  ameandata = as.double(b0.alpha),
                  ameanrow = as.integer(nrow(b0.alpha)),
                  ameancol = as.integer(ncol(b0.alpha)),   
                  bmeandata = as.double(b0.beta),
                  bmeanrow = as.integer(nrow(b0.beta)),
                  bmeancol = as.integer(ncol(b0.beta)), 
                  avardata = as.double(B0.alpha),
                  avarrow = as.integer(nrow(B0.alpha)),
                  avarcol = as.integer(ncol(B0.alpha)), 
                  bvardata = as.double(B0.beta),
                  bvarrow = as.integer(nrow(B0.beta)),
                  bvarcol = as.integer(ncol(B0.beta)), 
                  abcovdata = as.double(alpha.beta.cov),
                  abcovrow = as.integer(nrow(alpha.beta.cov)),
                  abcovcol = as.integer(ncol(alpha.beta.cov)),
                  store = as.integer(store.item),
                  PACKAGE="MCMCpack"
                  )
    
    theta.names <- paste("theta", 1:J, sep = "")
    alpha.beta.names <- paste(rep(c("alpha","beta"), K), rep(1:K, each = 2),
                              sep = "")
   
    # put together matrix and build MCMC object to return
    sample <- matrix(inv.obj$samdata, inv.obj$samrow, inv.obj$samcol,
                     byrow=TRUE)
    output <- mcmc2(data=sample,start=1, end=mcmc, thin=thin)
   
    if(store.item == FALSE) {
      names <- theta.names
    }
    else {
      names <- c(theta.names, alpha.beta.names)
    }
    varnames(output) <- names
    attr(output,"title") <-
      "MCMCpack One Dimensional IRT Model Posterior Density Sample"
    return(output)
    
  }
