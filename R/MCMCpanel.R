##########################################################################
# sample from the posterior distribution of general linear pane;
# model in R using linked C++ code in Scythe
#
# Andrew D. Martin
# Washington University
#
# Kevin M. Quinn
# University of Washington
#
# August 1, 2002
#
##########################################################################

MCMCpanel <- function(obs, Y, X, W, burnin = 1000, mcmc = 10000, thin = 5, 
    verbose = FALSE, seed = 0, beta.start = NA, sigma2.start = NA,
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
   # This program is free software; you can redistribute it and/or
   # modify it under the terms of the GNU General Public License
   # as published by the Free Software Foundation; either version 2
   # of the License, or (at your option) any later version.
   #
   # This program is distributed in the hope that it will be useful,
   # but WITHOUT ANY WARRANTY; without even the implied warranty of
   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   # GNU General Public License for more details.
   #  
   # For a copy of the GNU General Public License
   # write to the Free Software Foundation, Inc.,
   # 59 Temple Place - Suite 330, Boston, MA
   # 02111-1307, USA.  
   #
   # Copyright (C) 2002 Andrew D. Martin 
   #
   # Andrew D. Martin
   # Department of Political Science
   # Washington University
   # Campus Box 1063
   # One Brookings Drive
   # St. Louis, MO  63130
   # admartin@artsci.wustl.edu
   # 
   # 7/26/2002 
  
   # get starting time
   # start.time <- proc.time()

   # model parameters
   n <- length(unique(obs))
   k <- length(Y) / n
   p <- dim(X)[2]
   q <- dim(W)[2]

   # check data conformability
   if(dim(obs)[2] != 1) {
      cat("ERROR: obs is not a column vector.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
   }   
   if(length(unique(tabulate(obs))) != 1) {
      cat("ERROR: Panel is not balanced [check obs vector].\n")
      stop("Please respecify and call MCMCpanel() again.\n") 
   }
   if(dim(Y)[2] != 1) {
      cat("ERROR: Y is not a column vector.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
   }
   if(dim(X)[1] != n * k) {
      cat("ERROR: X matrix is not conformable [does not match Y].\n")
      stop("Please respecify and call MCMCpanel() again.\n")  
   }
   if(dim(W)[1] != n * k) {
      cat("ERROR: W matrix is not conformable [does not match Y].\n")
      stop("Please respecify and call MCMCpanel() again.\n")
   }   
 
   # check iteration parameters
   if(mcmc %% thin != 0) {
      cat("ERROR: MCMC interval not evenly divisible by thinning interval.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      } 
   if(mcmc < 0) {
      cat("ERROR: Gibbs interval negative.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      } 
   if(burnin < 0) {
      cat("ERROR: Burnin interval negative.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }       
   if(thin < 0) {
      cat("ERROR: Thinning interval negative.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }   
   totiter <- mcmc + burnin

   # starting values for beta error checking
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
      cat("ERROR: Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }
  
   # sigma2 starting values error checking
   if (is.na(sigma2.start)){
     sigma2.start <- ols.sigma2
   }   
   
   if(sigma2.start <= 0) {
      cat("ERROR: Starting value for sigma2 negative.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }   

   # starting values for D error checking
   if (is.na(D.start)){ # use matrix of ones
     D.start <- .5 * ols.sigma2 * diag(q)
   }   
   if(is.null(dim(D.start))) {
      D.start <- D.start * diag(q)
      }
   if((dim(D.start)[1] != q) || (dim(D.start)[2] != q)) {
      cat("ERROR: Starting value for D not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }

   # set up prior for beta
   if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,p,1)  
      }
   if((dim(b0)[1] != p) || (dim(b0)[2] != 1)) {
      cat("ERROR: N(b0,B0^-1) prior b0 not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }  
   if(is.null(dim(B0))) {
      B0 <- B0 * diag(p)
      }
   if((dim(B0)[1] != p) || (dim(B0)[2] != p)) {
      cat("ERROR: N(b0,B0^-1) prior B0 not conformable.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      } 
   
   # set up prior for sigma2
   if(nu0 <= 0) {
      cat("ERROR: G(nu0,delta0) prior nu0 less than or equal to zero.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }
   if(delta0 <= 0) {
      cat("ERROR: G(nu0,delta0) prior delta0 less than or equal to zero.\n")
      stop("Please respecify and call MCMCpanel() again.\n")      
      }
       
   # set up prior for D
   if(eta0 < q) {
      cat("ERROR: Wishart(eta0,R0) prior eta0 less than or equal to q.\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }   
   if(is.null(dim(R0))) {
      R0 <- R0 * diag(q)
      }
   if((dim(R0)[1] != q) || (dim(R0)[2] != q)) {
      cat("ERROR: Wishart(eta0,R0) prior R0 not comformable [q times q].\n")
      stop("Please respecify and call MCMCpanel() again.\n")
      }
         
   # set up big holder matrix
   sample <- matrix(0, mcmc/thin, p+q*q+1)
   
   # call C++ code to draw sample
   inv.obj <- .C("panelpost",
   samdata = as.double(sample),
   samrow = as.integer(nrow(sample)),
   samcol = as.integer(ncol(sample)),
   obsdata = as.double(obs),
   obsrow = as.integer(nrow(obs)),
   obscol = as.integer(ncol(obs)),   
   ydata = as.double(Y),
   yrow = as.integer(nrow(Y)),
   ycol = as.integer(ncol(Y)),   
   xdata = as.double(X),
   xrow = as.integer(nrow(X)),
   xcol = as.integer(ncol(X)),   
   wdata = as.double(W),
   wrow = as.integer(nrow(W)),
   wcol = as.integer(ncol(W)),   
   burnin = as.integer(burnin),
   gibbs = as.integer(mcmc),
   thin = as.integer(thin),
   seed = as.integer(seed),
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
  
   beta.names <- paste("beta", 1:p, sep = "")
   D.names <- paste("D", 1:(q*q), sep = "")
   sigma2.names <- "sigma2"
   names <- c(beta.names, D.names, sigma2.names)   
   output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
   varnames(output) <- names
   attr(output,"title") <- 
      "MCMCpack Linear Panel Model Posterior Density Sample"
   

   return(output)
}

