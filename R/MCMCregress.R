# sample from the posterior distribution of a Gaussian linear regression
# model in R using linked C++ code in Scythe
#
# ADM and KQ 5/21/2002

"MCMCregress" <-
  function(formula, data = list(), burnin = 1000, mcmc = 10000,
   thin=5, verbose = FALSE, seed = 0, sigma2.start = NA,
   b0 = 0, B0 = 0, nu = 0.001, delta = 0.001, ...) {
    
    # extract X, Y, and variable names from the model formula and frame       
    call <- match.call()
    mt <- terms(formula, data=data)
    if(missing(data)) data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$seed <- mf$verbose <- mf$beta.start <- mf$sigma2.start <- NULL
    mf$burnin <- mf$mcmc <- mf$thin <- NULL
    mf$b0 <- mf$B0 <- mf$nu <- mf$delta <- mf$... <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
  
    # null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    X <- as.matrix(X)         # X matrix
    xvars <- dimnames(X)[[2]] # X variable names
    Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
    N <- nrow(X)	             # number of observations      
    K <- ncol(X)              # number of covariates
   
    # burnin / mcmc / thin error checking
    check.parameters(burnin, mcmc, thin, "MCMCregress")
    
    # starting values for beta error checking
    beta.start <- NA
    if (is.na(beta.start)){ # use MLEs
      beta.start <- matrix(coef(lm(formula, data=data)), K, 1)
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,K,1)  
    }
    else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
      cat("Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCregress() again.\n")
    }

       # sigma2 starting values error checking
    if (is.na(sigma2.start)){ # use MLE
      lm.out <- lm(formula, data=data)
      sigma2.start <- var(residuals(lm.out))
    }
    else if(sigma2.start <= 0) {
      cat("Starting value for sigma2 negative.\n")
      stop("Please respecify and call MCMCregress() again.\n")
    }   
    
    # prior for beta error checking
    if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,K,1)  
    }
    if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
      cat("N(b0,B0^-1) prior b0 not conformable.\n")
      cat("Please respecify and call MCMCregress() again.\n") 
    }  
    if(is.null(dim(B0))) {
      B0 <- B0 * diag(K)    
    }
    if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
      cat("N(b0,B0^-1) prior B0 not conformable.\n")
      stop("Please respecify and call MCMCregress() again.\n")
    }  
   
    # prior for sigma2 error checking
    if(nu <= 0) {
      cat("IG(nu/2,delta/2) prior nu less than or equal to zero.\n")
      stop("Please respecify and call MCMCregress() again.\n")
    }
    if(delta <= 0) {
      cat("IG(nu/2,delta/2) prior delta less than or equal to zero.\n")
      stop("Please respecify and call MCMCregress() again.\n")      
    }  
   
    # define holder for posterior density sample
    sample <- matrix(data=0, mcmc/thin, dim(X)[2] + 1)
  
    # call C++ code to draw sample
    inv.obj <- .C("regpost",
                  samdata = as.double(sample),
                  samrow = as.integer(nrow(sample)),
                  samcol = as.integer(ncol(sample)),
                  Xdata = as.double(X),
                  Xrow = as.integer(nrow(X)),
                  Xcol = as.integer(ncol(X)),   
                  Ydata = as.double(Y),
                  Yrow = as.integer(nrow(Y)),
                  Ycol = as.integer(ncol(Y)),   
                  burnin = as.integer(burnin),
                  gibbs = as.integer(mcmc),
                  thin = as.integer(thin),
                  seed = as.integer(seed),
                  verbose = as.integer(verbose),
                  bstartdata = as.double(beta.start),
                  bstartrow = as.integer(nrow(beta.start)),
                  bstartcol = as.integer(ncol(beta.start)),
                  sigma2start = as.double(sigma2.start),
                  b0data = as.double(b0),
                  b0row = as.integer(nrow(b0)),
                  b0col = as.integer(ncol(b0)),   
                  B0data = as.double(B0),
                  B0row = as.integer(nrow(B0)),
                  B0col = as.integer(ncol(B0)),   
                  nu = as.double(nu),
                  delta = as.double(delta),
                  PACKAGE="MCMCpack"
                  )
   
    # put together matrix and build MCMC object to return
    sample <- matrix(inv.obj$samdata, inv.obj$samrow, inv.obj$samcol, byrow=TRUE)
    output <- mcmc2(data=sample,start=1, end=mcmc, thin=thin)
    names <- c(xvars, "sigma2")
    varnames(output) <- names
    attr(output,"title") <- "MCMCregress Posterior Density Sample"
    return(output)   
  }
