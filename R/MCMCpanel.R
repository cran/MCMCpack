# sample from the posterior distribution of general linear panel
# model in R using linked C++ code in Scythe
#
# NOTE: this does not take a data argument because only
#       matrices are passed.  This should probably be fixed.  Also,
#       re-implementing using something likes Bates' syntax
#       would be nice.  The helper functions could also be used 
#       all over the place here. This is another good project for a grad
#       student.
#
#
# ADM and KQ 8/1/2002
# updated with Ben Goodrich's feedback and new spec ADM 7/28/2004


"MCMCpanel" <-
  function(obs, Y, X, W, burnin = 1000, mcmc = 10000, thin = 5, 
           verbose = 0, seed = NA, sigma2.start = NA,
           D.start = NA, b0 = 0, B0 = 1, eta0, R0, nu0 = 0.001,
           delta0 = 0.001, ...) {

    # DESCRIPTION:
    #
    #   MCMCpanel fits a general linear panel model using Algorithm 2 of
    #   Chib and Carlin (1999).  The program calls a compiled C++ shared
    #   library to perform the actual sampling. The model takes the
    #   following form:
    #
    #      y_i = X_i \beta + W_i b_i + \varepsilon_i
    #
    #      b_i \sim N_q(0,D)
    #
    #      \varepsilon_i \sim N_k(0,\sigma^2 I_n)
    #
    #   With conjugate priors:
    #
    #      \beta \sim N_p(\beta_0, \B_0^-1)
    #
    #      D^-1 \sim Wishart(\nu_0^{-1} R_0, \nu_0)
    #
    #      \sigma^-2 \sim Gamma(\nu_00/2, \delta_00/2) 
    #
    #   The model is defined in terms of k (the number of responses
    #   per subject, assumed to be constant across subjects), p (the
    #   number of columns in the design matrix of covariates), and 
    #   q (the number of columns in the design matrix), and n (the
    #   number of subjects).  The components of the model are the
    #   following:
    #
    #      y_i (k \times 1) vector of responses for subject i
    #
    #      X_i (k \times p) matrix of covariates for subject i      
    #
    #      \beta (p \times 1) vector of fixed effects coefficients
    #
    #      W_i (k \times q) design matrix for random effects for subject i
    #
    #      b_i (q \times 1) vector of random effects for subject i
    #
    #      \varepsilon (k \times 1) vector of errors for subject i
    #

    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    # model parameters
    n <- length(unique(obs))
    k <- length(Y) / n
    p <- dim(X)[2]
    q <- dim(W)[2]

    # check data conformability
    obs.temp <- as.matrix(obs)
    if (any(obs.temp[,1] != obs)) {
      cat("Error: obs is not a column vector.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }   
    if(length(unique(tabulate(obs))) != 1) {
      cat("Error: Panel is not balanced [check obs vector].\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }
    Y.temp <- as.matrix(Y)
    if (any(Y.temp[,1] != Y)) {
      cat("Error: X matrix is not conformable [does not match Y].\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }
    if(dim(W)[1] != n * k) {
      cat("Error: W matrix is not conformable [does not match Y].\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }   
 
    # check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin

    # starting values for beta error checking
    beta.start <- NA
    ols.beta <- solve(t(X) %*% X) %*% t(X) %*% Y
    ols.sigma2 <-
      t(Y - X %*% ols.beta) %*% (Y - X %*% ols.beta) / (k*n - p - 1)
    ols.sigma2 <- as.double(ols.sigma2)
    if (is.na(beta.start)){ # use least squares estimates
      beta.start <- ols.beta
    }
    if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,p,1)  
    }
    if((dim(beta.start)[1] != p) || (dim(beta.start)[2] != 1)) {
      cat("Error: Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }
  
    # sigma2 starting values error checking
    if (is.na(sigma2.start)){
      sigma2.start <- ols.sigma2
    }   
    if(sigma2.start <= 0) {
      cat("Error: Starting value for sigma2 negative.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }   

    # starting values for D error checking
    if (is.na(D.start)){ # use matrix of ones
      D.start <- .5 * ols.sigma2 * diag(q)
    }   
    if(is.null(dim(D.start))) {
      D.start <- D.start * diag(q)
    }
    if((dim(D.start)[1] != q) || (dim(D.start)[2] != q)) {
      cat("Error: Starting value for D not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }

    # set up prior for beta
    if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,p,1)  
    }
    if((dim(b0)[1] != p) || (dim(b0)[2] != 1)) {
      cat("Error: N(b0,B0^-1) prior b0 not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }  
    if(is.null(dim(B0))) {
      B0 <- B0 * diag(p)
    }
    if((dim(B0)[1] != p) || (dim(B0)[2] != p)) {
      cat("Error: N(b0,B0^-1) prior B0 not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    } 
   
    # set up prior for sigma2
    if(nu0 <= 0) {
      cat("Error: G(nu0,delta0) prior nu0 less than or equal to zero.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }
    if(delta0 <= 0) {
      cat("Error: G(nu0,delta0) prior delta0 less than or equal to zero.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }
       
    # set up prior for D
    if(eta0 < q) {
      cat("Error: Wishart(eta0,R0) prior eta0 less than or equal to q.\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }   
    if(is.null(dim(R0))) {
      R0 <- R0 * diag(q)
    }
    if((dim(R0)[1] != q) || (dim(R0)[2] != q)) {
      cat("Error: Wishart(eta0,R0) prior R0 not comformable [q times q].\n")
      stop("Please respecify and call MCMCpanel() again.\n", call.=FALSE)
    }
         
    # set up big holder matrix
    sample <- matrix(0, mcmc/thin, p+q*q+1)
   
    # call C++ code to draw sample
    inv.obj <- .C("panelpost",
                  samdata = as.double(sample),
                  samrow = as.integer(nrow(sample)),
                  samcol = as.integer(ncol(sample)),
                  obsdata = as.double(obs),
                  obsrow = as.integer(length(obs)),
                  obscol = as.integer(1),   
                  ydata = as.double(Y),
                  yrow = as.integer(length(Y)),
                  ycol = as.integer(1),   
                  xdata = as.double(X),
                  xrow = as.integer(nrow(X)),
                  xcol = as.integer(ncol(X)),   
                  wdata = as.double(W),
                  wrow = as.integer(nrow(W)),
                  wcol = as.integer(ncol(W)),   
                  burnin = as.integer(burnin),
                  gibbs = as.integer(mcmc),
                  thin = as.integer(thin),
                  lecuyer = as.integer(lecuyer),
                  seedarray = as.integer(seed.array),
                  lecuyerstream = as.integer(lecuyer.stream),
                  verbose = as.integer(verbose),
                  bstartdata = as.double(beta.start),
                  bstartrow = as.integer(nrow(beta.start)),
                  bstartcol = as.integer(ncol(beta.start)),
                  sigma2start = as.double(sigma2.start),
                  Dstartdata = as.double(D.start),
                  Dstartrow = as.integer(nrow(D.start)),
                  Dstartcol = as.integer(ncol(D.start)),
                  b0data = as.double(b0),
                  b0row = as.integer(nrow(b0)),
                  b0col = as.integer(ncol(b0)),   
                  B0data = as.double(B0),
                  B0row = as.integer(nrow(B0)),
                  B0col = as.integer(ncol(B0)),  
                  nu0 = as.double(nu0),
                  delta0 = as.double(delta0),
                  eta0 = as.integer(eta0),    
                  R0data = as.double(R0),
                  R0row = as.integer(nrow(R0)),
                  R0col = as.integer(ncol(R0)),
                  n = as.integer(n),
                  k = as.integer(k),
                  p = as.integer(p),
                  q = as.integer(q),
                  PACKAGE="MCMCpack"
                  )     
    
    # put together matrix and build MCMC object to return
    sample <- matrix(inv.obj$samdata, inv.obj$samrow, inv.obj$samcol,
                     byrow=TRUE)   
  
    if (length(colnames(X))>0) {
        beta.names <- colnames(X)
    }
    else beta.names <- paste("beta", 1:p, sep = "")
    D.names <- paste("D", 1:(q*q), sep = "")
    sigma2.names <- "sigma2"
    names <- c(beta.names, D.names, sigma2.names)   

    output <- mcmc(data=sample, start=1, end=mcmc, thin=thin)
    varnames(output) <- names
    attr(output,"title") <- 
      "MCMCpack Linear Panel Model Posterior Density Sample"
    return(output)
  }

