# sample from the posterior distribution of an ordered probit model
# via the data augmentation approach of Cowles (1996)
#
# KQ 1/25/2003

"MCMCoprobit" <-
  function(formula, data = list(), burnin = 1000, mcmc = 10000,
           thin = 5, tune = NA, verbose = FALSE, seed = 0, beta.start = NA,
           b0 = 0, B0 = 0.001, ...) {
   
   # extract X, Y, and variable names from the model formula and frame
   call <- match.call()
   mt <- terms(formula, data=data)
   if(missing(data)) data <- sys.frame(sys.parent())
   mf <- match.call(expand.dots = FALSE)
   mf$burnin <- mf$mcmc <- mf$b0 <- mf$B0 <- NULL
   mf$thin <- mf$... <- mf$tune <- mf$verbose <- mf$seed <- NULL
   mf$beta.start <- NULL
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent()))
   vars <- as.character(attr(mt, "variables"))[-1] # y varname and x varnames
   
   # null model support
   X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)# else NULL
   X.names <- dimnames(X)[[2]]
   Y    <- model.response(mf, "numeric")
   Y    <- factor(Y, ordered=TRUE)
   ncat <- nlevels(Y)             # number of categories of y
   cat  <- levels(Y)              # values of categories of y
   N <- nrow(X)	                  # number of observations
   K <- ncol(X)	                  # number of covariates
   
   # convert data to matrices to be passed
   Y <- as.matrix(as.integer(Y))
   X <- as.matrix(X)
   
   # burnin / mcmc / thin error checking
   check.parameters(burnin, mcmc, thin, "MCMCoprobit", tune)
    
   # check tuning parameter
   if (is.na(tune)){
     tune <- 0.05/ncat
   }
 
   xint <- match("(Intercept)", colnames(X), nomatch=0)
   if (xint > 0){
     new.X <- X[, -xint, drop=FALSE]
   }
   else warning("An intercept is needed and assumed in MCMCoprobit()\n.")
   if (ncol(new.X) == 0){
     polr.out <- polr(ordered(Y)~1)
   }
   else {
     polr.out <- polr(ordered(Y)~new.X)
   }
   
   # starting values for beta error checking
   if (is.na(beta.start)){
     beta.start <- matrix(0, K, 1)
     beta.start[1] <- -.588 * polr.out$zeta[1]
     if( ncol(new.X) > 0){
       beta.start[2:K] <- .588 * coef(polr.out)
     }
   }
   else if(is.null(dim(beta.start))) {
     beta.start <- beta.start * matrix(1,K,1)  
   }
   else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
     cat("Starting value for beta not conformable.\n")
     stop("Please respecify and call MCMCoprobit() again.\n")
   }
   
   # prior for beta error checking
   if(is.null(dim(b0))) {
     b0 <- b0 * matrix(1,K,1)  
   }
   if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
     cat("N(b0,B0) prior b0 not conformable.\n")
     stop("Please respecify and call MCMCoprobit() again.\n") 
   }   
   if(is.null(dim(B0))) {
     B0 <- B0 * diag(K)    
   }
   if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
     cat("N(b0,B0) prior B0 not conformable.\n")
     stop("Please respecify and call MCMCoprobit() again.\n")
   }

   # form gamma starting values (note: not changeable)
   gamma <- matrix(NA,ncat+1,1)
   gamma[1] <- -300
   gamma[2] <- 0
   gamma[3:ncat] <- (polr.out$zeta[2:(ncat-1)] - polr.out$zeta[1])*.588
   gamma[ncat+1] <- 300
      
   # posterior density sample
   sample <- matrix(data=0, mcmc/thin, K + ncat + 1)
   
   posterior <- .C("oprobitpost",
                   samdata = as.double(sample),
                   samrow = as.integer(nrow(sample)),
                   samcol = as.integer(ncol(sample)),
                   Ydata = as.integer(Y),
                   Yrow = as.integer(nrow(Y)),
                   Ycol = as.integer(ncol(Y)),   
                   Xdata = as.double(X),
                   Xrow = as.integer(nrow(X)),
                   Xcol = as.integer(ncol(X)),   
                   burnin = as.integer(burnin),
                   mcmc = as.integer(mcmc),
                   thin = as.integer(thin),
                   tune = as.double(tune),
                   seed = as.integer(seed),
                   verbose = as.integer(verbose),
                   bstartdata = as.double(beta.start),
                   bstartrow = as.integer(nrow(beta.start)),
                   bstartcol = as.integer(ncol(beta.start)),
                   gammadata = as.double(gamma),
                   gammarow = as.integer(nrow(gamma)),
                   gammacol = as.integer(ncol(gamma)),
                   b0data = as.double(b0),
                   b0row = as.integer(nrow(b0)),
                   b0col = as.integer(ncol(b0)),   
                   B0data = as.double(B0),
                   B0row = as.integer(nrow(B0)),
                   B0col = as.integer(ncol(B0)),
                   accepts = as.integer(0),
                   PACKAGE="MCMCpack"
                   )

   cat(" Overall acceptance rate = ",
       posterior$accepts / (posterior$burnin+posterior$mcmc), "\n")
   
   
   # put together matrix and build MCMC object to return
   sample <- matrix(posterior$samdata, posterior$samrow,
                    posterior$samcol, byrow=TRUE)
   sample <- sample[,c(1:K, (K+3):(K+ncat))]

   output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
   names <- c(X.names, paste("gamma", 2:(ncat-1), sep=""))     
   varnames(output) <- names
   attr(output,"title") <-
     "MCMCpack Ordered Probit Posterior Density Sample"   
   return(output)
 }
