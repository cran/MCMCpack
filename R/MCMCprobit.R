##########################################################################
# sample from the posterior distribution of a Gaussian linear regression
# model in R using linked C++ code in Scythe
#
# Andrew D. Martin
# Washington University
#
# Kevin M. Quinn
# University of Washington
#
# May 21, 2002
#
##########################################################################

MCMCprobit <- function(formula, data = list(), burnin = 1000, mcmc = 10000,
   thin = 5, verbose = FALSE, seed = 0, beta.start = NA,
   b0 = 0, B0 = 0, bayes.resid=FALSE, ...) {
  
   #
   # extract X, Y, and variable names from the model formula and frame       
   #
   call <- match.call()
   mt <- terms(formula, data=data)
   if(missing(data)) data <- sys.frame(sys.parent())
   mf <- match.call(expand.dots = FALSE)
   mf$seed <- mf$verbose <- mf$beta.start <- mf$bayes.resid <- NULL
   mf$burnin <- mf$mcmc <- mf$thin <- NULL
   mf$b0 <- mf$B0 <- mf$... <- NULL
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
      

   # residuals setup
   if (is.logical(bayes.resid) && bayes.resid==TRUE){
     resvec <- 1:N
   }
   else if (bayes.resid != FALSE){
     resvec <- bayes.resid
   }

   
   # y \in {0, 1} error checking
   if (sum(Y!=0 & Y!=1) > 0) {
     cat("Elements of Y equal to something other than 0 or 1.")
     stop("\n Check data and call MCMCprobit() again. \n") 
     }
   
   # burnin / mcmc / thin error checking
   if(mcmc %% thin != 0) {
      cat("Gibbs interval not evenly divisible by thinning interval.")
      stop("\n Please respecify and call MCMCprobit() again.\n")
      } 
   if(mcmc < 0) {
      cat("Gibbs interval negative.")
      stop("\n Please respecify and call MCMCprobit() again.\n")
      } 
   if(burnin < 0) {
      cat("Burnin interval negative.")
      stop("\n Please respecify and call MCMCprobit() again.\n")
      }       
   if(thin < 0) {
      cat("Thinning interval negative.")
      stop("\n Please respecify and call MCMCprobit() again.\n")
      }          
   
   # starting values for beta error checking
   if (is.na(beta.start)){ # use MLEs
     beta.start <- matrix(coef(glm(formula, data=data,
                                  family=binomial(link=probit))), K, 1)
   }
   else if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,K,1)  
      }
   else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
      cat("Starting value for beta not conformable.")
      stop("\n Please respecify and call MCMCprobit() again.\n")
      }

   # prior for beta error checking
   if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,K,1)  
    }
   if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
     cat("N(b0,B0^-1) prior b0 not conformable.")
     stop("\n Please respecify and call MCMCprobit() again.\n")
   }  
   
   if(is.null(dim(B0))) {
     B0 <- B0 * diag(K)    
   }
   if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
     cat("N(b0,B0^-1) prior B0 not conformable.")
     stop("\n Please respecify and call MCMCprobit() again.\n")
   }  
   


   if (bayes.resid == FALSE){
   
     # define holder for posterior density sample
     sample <- matrix(data=0, mcmc/thin, dim(X)[2] )
  
     # call C++ code to draw sample
     posterior <- .C("probitpost",
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
                     b0data = as.double(b0),
                     b0row = as.integer(nrow(b0)),
                     b0col = as.integer(ncol(b0)),   
                     B0data = as.double(B0),
                     B0row = as.integer(nrow(B0)),
                     B0col = as.integer(ncol(B0)),   
                     PACKAGE="MCMCpack"
                     )
   
     # put together matrix and build MCMC object to return
     sample <- matrix(posterior$samdata, posterior$samrow,
                      posterior$samcol, byrow=TRUE)
     output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
     names <- c(xvars)
     varnames(output) <- names
     attr(output,"title") <-
        "MCMCpack Probit Regression Posterior Density Sample"
   }
   else{
     # define holder for posterior density sample
     sample <- matrix(data=0, mcmc/thin, dim(X)[2]+length(resvec) )
  
     # call C++ code to draw sample
     posterior <- .C("probitpostres",
                     samdata = as.double(sample),
                     samrow = as.integer(nrow(sample)),
                     samcol = as.integer(ncol(sample)),
                     Xdata = as.double(X),
                     Xrow = as.integer(nrow(X)),
                     Xcol = as.integer(ncol(X)),   
                     Ydata = as.double(Y),
                     Yrow = as.integer(nrow(Y)),
                     Ycol = as.integer(ncol(Y)),
                     resvecdata = as.double(resvec),
                     resvecrow = as.integer(length(resvec)),
                     resveccol = as.integer(1),
                     burnin = as.integer(burnin),
                     gibbs = as.integer(mcmc),
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
                     PACKAGE="MCMCpack"   
                     )
     # put together matrix and build MCMC object to return
     sample <- matrix(posterior$samdata, posterior$samrow,
                      posterior$samcol, byrow=TRUE)
     output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
     names <- c(xvars, paste("epsilonstar", as.character(resvec), sep="") )
     varnames(output) <- names
     attr(output,"title") <-
        "MCMCpack Probit Regression Posterior Density Sample"
   }
     return(output)

     
}

##########################################################################
