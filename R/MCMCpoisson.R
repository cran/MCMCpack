##########################################################################
#
# sample from the posterior distribution of a Poisson regression
# model in R using linked C++ code in Scythe
#
# Andrew D. Martin
# Washington University
#
# Kevin M. Quinn
# University of Washington
#
# January 24, 2003
#
# Bug fixed by KQ on 3/17/2003
#
##########################################################################

MCMCpoisson <- function(formula, data = list(), burnin = 1000, mcmc = 10000,
   thin=5, tune=1.1, verbose = FALSE, seed = 0, beta.start = NA,
   b0 = 0, B0 = 0.001, ...) {
  
   #
   # extract X, Y, and variable names from the model formula and frame       
   #
   call <- match.call()
   mt <- terms(formula, data=data)
   if(missing(data)) data <- sys.frame(sys.parent())
   mf <- match.call(expand.dots = FALSE)
   mf$seed <- mf$verbose <- mf$beta.start <- NULL
   mf$burnin <- mf$mcmc <- mf$thin <- mf$tune <- NULL
   mf$b0 <- mf$B0 <- mf$... <- NULL
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent()))
  
   # null model support
   X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
   X <- as.matrix(X)         # X matrix
   xvars <- dimnames(X)[[2]] # X variable names
   Y <- as.matrix(as.integer(model.response(mf, "numeric"))) # Y matrix
   N <- nrow(X)	             # number of observations      
   K <- ncol(X)              # number of covariates
     
   # test y non-negative
   if (sum(Y < 0) > 0) {
      cat("\n Elements of Y negative. ")
      stop("\n Check data and call MCMCpoisson() again. \n") 
     }
   
   # burnin / mcmc / thin error checking
   if(mcmc %% thin != 0) {
      cat("mcmc interval not evenly divisible by thinning interval.")
      stop("\n Please respecify and call MCMCpoisson() again.\n")
      }
   if(mcmc < 0) {
      cat("mcmc interval negative.")
      stop("\n Please respecify and call MCMCpoisson() again.\n")
      }
   if(burnin < 0) {
      cat("Burnin interval negative.\n")
      stop("Please respecify and call MCMCpoisson() again.\n")
      }
   if(thin < 0) {
      cat("Thinning interval negative.\n")
      stop("Please respecify and call MCMCpoisson() again.\n")
      }  
   
   # starting values for beta error checking
   library(MASS)
   glm.out <- glm(formula, data=data, family=poisson())
   m <- coef(glm.out)
   V <- vcov(glm.out)
   if (is.na(beta.start)) { # use MLEs
     beta.start <- matrix(m, K, 1)
     }
   else if(is.null(dim(beta.start))) {
     beta.start <- beta.start * matrix(1,K,1)  
     }
   else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
     cat("Starting value for beta not conformable.")
     stop("\n Please respecify and call MCMCpoisson() again.\n")
     }

   # prior for beta error checking
   if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,K,1)  
      }
   if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
     cat("N(b0,B0) prior b0 not conformable.")
     stop("\n Please respecify and call MCMCpoisson() again.\n")
     }
   if(is.null(dim(B0))) {
     B0 <- B0 * diag(K)    
     }
   if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
     cat("N(b0,B0) prior B0 not conformable.")
     stop("\n Please respecify and call MCMCpoisson() again.\n") 
     }
   
   # define holder for posterior density sample
   sample <- matrix(data=0, mcmc/thin, dim(X)[2] )
  
   # call C++ code to draw sample
   posterior <- .C("poissonpost",
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
                   mcmc = as.integer(mcmc),
                   thin = as.integer(thin),
                   seed = as.integer(seed),
                   verbose = as.integer(verbose),
                   bstartdata = as.double(beta.start),
                   bstartrow = as.integer(nrow(beta.start)),
                   bstartcol = as.integer(ncol(beta.start)),
                   b0data = as.double(b0),
                   b0row = as.integer(nrow(b0)),
                   b0col = as.integer(ncol(b0)),   
                   B0data = as.double(B0),
                   B0row = as.integer(nrow(B0)),
                   B0col = as.integer(ncol(B0)),
                   mdata = as.double(m),
                   mrow = as.integer(length(m)),
                   mcol = as.integer(1),
                   Vdata = as.double(V),
                   Vrow = as.integer(nrow(V)),
                   Vcol = as.integer(ncol(V)),
                   tune = as.double(tune),
                   accepts = as.integer(0),
                   PACKAGE="MCMCpack"
                   )
   cat("overall acceptance rate = ",
      posterior$accepts / (posterior$burnin+posterior$mcmc), "\n")
  
   # put together matrix and build MCMC object to return
   sample <- matrix(posterior$samdata, posterior$samrow,
                    posterior$samcol, byrow=TRUE)
   output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
   names <- c(xvars)
   varnames(output) <- names
   attr(output,"title") <-
     "MCMCpack Poisson Regression Posterior Density Sample"
   return(output)

}

##########################################################################
