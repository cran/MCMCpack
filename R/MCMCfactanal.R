##########################################################################
# sample from the posterior distribution of a factor analysis model
# model in R using linked C++ code in Scythe.
#
# The model is:
#
# x_i = \Lambda \phi_i + \epsilon_i,   \epsilon_i \sim N(0, \Psi)
#
# where \Psi is diagonal and the priors are:
#
# \lambda_{ij} \sim N(l_{ij}, L^{-1}_{ij})
# \phi_i \sim N(0,I)
# \psi^2_{jj} \sim IG(a0_j/2, b0_j/2)
#
#
# Andrew D. Martin
# Washington University
#
# Kevin M. Quinn
# Harvard University
#
# May 7, 2003
#
##########################################################################

"MCMCfactanal" <-
  function(x, factors, lambda.constraints=list(),
           data=list(), burnin = 1000, mcmc = 10000,
           thin=5, verbose = FALSE, seed = 0,
           lambda.start = NA, psi.start = NA,
           l0=0, L0=0, a0=0.001, b0=0.001,
           store.scores = FALSE, std.var=TRUE, ... ) {
    
    # extract X and variable names from the model formula and frame       
    if (is.matrix(x)){
      X <- x
      xvars <- dimnames(X)[[2]]
      xobs <- dimnames(X)[[1]]
      N <- nrow(X)
      K <- ncol(X)
      if (std.var){
        for (i in 1:K){
          X[,i] <- (X[,i]-mean(X[,i]))/sd(X[,i])
        }
      }
      else{
        for (i in 1:K){
          X[,i] <- X[,i]-mean(X[,i])
        }
      }
    }
    else {
      call <- match.call()
      mt <- terms(x, data=data)
      if (attr(mt, "response") > 0) 
        stop("Response not allowed in formula in MCMCfactanal().\n")
      if(missing(data)) data <- sys.frame(sys.parent())
      mf <- match.call(expand.dots = FALSE)
      mf$factors <- mf$lambda.constraints <- mf$burnin <- mf$mcmc <- NULL
      mf$thin <- mf$verbose <- mf$seed <- NULL
      mf$lambda.start <- mf$psi.start <- mf$l0 <- mf$L0 <- NULL
      mf$a0 <- mf$b0 <- mf$store.scores <- mf$std.var <- mf$... <- NULL
      mf$drop.unused.levels <- TRUE
      mf[[1]] <- as.name("model.frame")
      mf <- eval(mf, sys.frame(sys.parent()))
      
      attributes(mt)$intercept <- 0
      X <- model.matrix(mt, mf, contrasts)
      X <- as.matrix(X)         # X matrix
      xvars <- dimnames(X)[[2]] # X variable names
      xobs <- dimnames(X)[[1]]  # observation names
      N <- nrow(X)	      # number of observations      
      K <- ncol(X)              # number of manifest variables

      # standardize X
      if (std.var){
        for (i in 1:K){
          X[,i] <- (X[,i]-mean(X[,i]))/sd(X[,i])
        }
      }
      else{
        for (i in 1:K){
          X[,i] <- X[,i]-mean(X[,i])
        }
      }
    }
    
    if (is.null(xobs)){
      xobs <- 1:N
    }

    check.parameters(burnin, mcmc, thin, "MCMCfactanal")
      
    # give names to the rows of Lambda related matrices
    Lambda.eq.constraints <- matrix(NA, K, factors)
    Lambda.ineq.constraints <- matrix(0, K, factors)
    Lambda.prior.mean <- matrix(0, K, factors)
    Lambda.prior.prec <- matrix(0, K, factors)

  
    if (is.null(colnames(X))){
      rownames(Lambda.eq.constraints) <- paste("V", 1:ncol(X), sep="") 
      rownames(Lambda.ineq.constraints) <- paste("V", 1:ncol(X), sep="")
      rownames(Lambda.prior.mean) <- paste("V", 1:ncol(X), sep="")
      rownames(Lambda.prior.prec) <- paste("V", 1:ncol(X), sep="")
      X.names <- paste("V", 1:ncol(X), sep="")
    }
    if (!is.null(colnames(X))){
      rownames(Lambda.eq.constraints) <- colnames(X)
      rownames(Lambda.ineq.constraints) <- colnames(X)
      rownames(Lambda.prior.mean) <- colnames(X)
      rownames(Lambda.prior.prec) <- colnames(X)
      X.names <- colnames(X)
    }
    
    # setup the equality and inequality contraints on Lambda
    if (length(lambda.constraints) != 0){
      constraint.names <- names(lambda.constraints)  
      for (i in 1:length(constraint.names)){
        name.i <- constraint.names[i]
        lambda.constraints.i <- lambda.constraints[[i]]
        col.index <- lambda.constraints.i[[1]]
        replace.element <- lambda.constraints.i[[2]]
        if (is.numeric(replace.element)){
          Lambda.eq.constraints[rownames(Lambda.eq.constraints)==name.i,
                                col.index] <- replace.element
        }
        if (replace.element=="+"){
          Lambda.ineq.constraints[rownames(Lambda.ineq.constraints)==name.i,
                                  col.index] <- 1
        }
        if (replace.element=="-"){
          Lambda.ineq.constraints[rownames(Lambda.ineq.constraints)==name.i,
                                  col.index] <- -1
        }
      }
    }
    
    testmat <- Lambda.ineq.constraints * Lambda.eq.constraints

    if (min(is.na(testmat))==0){
      if ( min(testmat[!is.na(testmat)]) < 0){
        cat("Constraints on Lambda are logically inconsistent.\n")
        stop("Please respecify and call MCMCfactanal() again.\n")
      }
    }
    Lambda.eq.constraints[is.na(Lambda.eq.constraints)] <- -999
    
    
    # setup prior means and precisions for Lambda
    # prior means
    if (is.matrix(l0)){ # matrix input for l0
      if (nrow(l0)==K && ncol(l0)==factors)
        Lambda.prior.mean <- l0
      else {
        cat("l0 not of correct size for model specification.\n")
        stop("Please respecify and call MCMCfactanl() again.\n")
      }
    }
    else if (is.list(l0)){ # list input for l0
      l0.names <- names(l0)
      for (i in 1:length(l0.names)){
        name.i <- l0.names[i]
        l0.i <- l0[[i]]
        col.index <- l0.i[[1]]
        replace.element <- l0.i[[2]]
        if (is.numeric(replace.element)){
          Lambda.prior.mean[rownames(Lambda.prior.mean)==name.i,
                            col.index] <- replace.element
        }   
      }
    }
    else if (length(l0)==1 && is.numeric(l0)){ # scalar input for l0
      Lambda.prior.mean <- matrix(l0, K, factors)
    }
    else {
      cat("l0 neither matrix, list, nor scalar.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    # prior precisions
    if (is.matrix(L0)){ # matrix input for L0
      if (nrow(L0)==K && ncol(L0)==factors)
        Lambda.prior.prec <- L0
      else {
        cat("L0 not of correct size for model specification.\n")
        stop("Please respecify and call MCMCfactanl() again.\n")
      }
    }
    else if (is.list(L0)){ # list input for L0
      L0.names <- names(L0)
      for (i in 1:length(L0.names)){
        name.i <- L0.names[i]
        L0.i <- L0[[i]]
        col.index <- L0.i[[1]]
        replace.element <- L0.i[[2]]
        if (is.numeric(replace.element)){
          Lambda.prior.prec[rownames(Lambda.prior.prec)==name.i,
                            col.index] <- replace.element
        }   
      }
    }
    else if (length(L0)==1 && is.numeric(L0)){ # scalar input for L0
      Lambda.prior.prec <- matrix(L0, K, factors)
    }
    else {
      cat("L0 neither matrix, list, nor scalar.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    if (min(L0) < 0) {
      cat("L0 contains negative elements.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    
    # Starting values for Lambda
    Lambda <- matrix(0, K, factors)
    if (is.na(lambda.start)){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:factors){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              Lambda[i,j] <- 0
            }
            if(Lambda.ineq.constraints[i,j]>0){
              Lambda[i,j] <- .5
            }
            if(Lambda.ineq.constraints[i,j]<0){
              Lambda[i,j] <- -.5
            }          
          }
          else Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else if (is.matrix(lambda.start)){
      if (nrow(lambda.start)==K && ncol(lambda.start)==factors)
        Lambda  <- lambda.start
      else {
        cat("lambda.start not of correct size for model specification.\n")
        stop("Please respecify and call MCMCfactanal() again.\n")
      }
    }
    else if (length(lambda.start)==1 && is.numeric(lambda.start)){
      Lambda <- matrix(lambda.start, K, factors)
      for (i in 1:K){
        for (j in 1:factors){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else {
      cat("lambda.start neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    
    # starting values for Psi
    if (is.na(psi.start)){
      Psi <- diag(K)*0.5
    }
    else if (is.double(psi.start))
      Psi <- diag(K) * psi.start
    else {
      cat("psi.start neither NA, nor double.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    if (nrow(Psi) != K || ncol(Psi) != K){
      cat("Psi starting value not K by K matrix.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")    
    }
    
    # setup prior for diag(Psi)
    if (length(a0)==1 && is.double(a0))
      a0 <- matrix(a0, K, 1)
    else if (length(a0) == K && is.double(a0))
      a0 <- matrix(a0, K, 1)
    else {
      cat("a0 not properly specified.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    if (length(b0)==1 && is.double(b0))
      b0 <- matrix(b0, K, 1)
    else if (length(b0) == K && is.double(b0))
      b0 <- matrix(b0, K, 1)
    else {
      cat("b0 not properly specified.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    
    # prior for Psi error checking
    if(min(a0) <= 0) {
      cat("IG(a0/2,b0/2) prior parameter a0 less than or equal to zero.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")
    }
    if(min(b0) <= 0) {
      cat("IG(a0/2,b0/2) prior parameter b0 less than or equal to zero.\n")
      stop("Please respecify and call MCMCfactanal() again.\n")      
    }  
    
    # define holder for posterior density sample
    if(store.scores == FALSE) {
      sample <- matrix(data=0, mcmc/thin, K*factors+K)
    }
    else {
      sample <- matrix(data=0, mcmc/thin, K*factors+K+factors*N)
    }

    # Call the C++ code to do the real work
    posterior <- .C("factanalpost",
                    samdata = as.double(sample),
                    samrow = as.integer(nrow(sample)),
                    samcol = as.integer(ncol(sample)),
                    X = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    seed = as.integer(seed),
                    verbose = as.integer(verbose),
                    Lambda = as.double(Lambda),
                    Lambdarow = as.integer(nrow(Lambda)),
                    Lambdacol = as.integer(ncol(Lambda)),
                    Psi = as.double(Psi),
                    Psirow = as.integer(nrow(Psi)),
                    Psicol = as.integer(ncol(Psi)),
                    Lameq = as.double(Lambda.eq.constraints),
                    Lameqrow = as.integer(nrow(Lambda.eq.constraints)),
                    Lameqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lamineq = as.double(Lambda.ineq.constraints),
                    Lamineqrow = as.integer(nrow(Lambda.ineq.constraints)),
                    Lamineqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lampmean = as.double(Lambda.prior.mean),
                    Lampmeanrow = as.integer(nrow(Lambda.prior.mean)),
                    Lampmeancol = as.integer(ncol(Lambda.prior.prec)),
                    Lampprec = as.double(Lambda.prior.prec),
                    Lampprecrow = as.integer(nrow(Lambda.prior.prec)),
                    Lamppreccol = as.integer(ncol(Lambda.prior.prec)),
                    a0 = as.double(a0),
                    a0row = as.integer(nrow(a0)),
                    a0col = as.integer(ncol(a0)),
                    b0 = as.double(b0),
                    b0row = as.integer(nrow(b0)),
                    b0col = as.integer(ncol(b0)),
                    storescores = as.integer(store.scores),
                    PACKAGE="MCMCpack"
                    )
    
    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$samdata, posterior$samrow, posterior$samcol,
                     byrow=TRUE)
    output <- mcmc2(data=sample,start=1, end=mcmc, thin=thin)
    
    Lambda.names <- paste(paste("Lambda",
                                rep(X.names,
                                    each=factors), sep=""),
                          rep(1:factors,K), sep="_")
    
    
    Psi.names <- paste("Psi", X.names, sep="")
    par.names <- c(Lambda.names, Psi.names)
    if (store.scores==TRUE){
      phi.names <- paste(paste("phi",
                               rep(xobs, each=factors), sep="_"),
                         rep(1:factors,factors), sep="_")
      par.names <- c(par.names, phi.names)
    }
    varnames(output) <- par.names

    # get rid of columns for constrained parameters
    output.df <- mcmc2dataframe(output)
    output.var <- diag(var(output.df))
    output.df <- output.df[,output.var != 0]
    
    output <- mcmc2(as.matrix(output.df), start=1, end=mcmc, thin=thin)

    # add constraint info so this isn't lost
    attr(output, "constraints") <- lambda.constraints
    attr(output, "n.manifest") <- K
    attr(output, "n.factors") <- factors
    
    attr(output,"title") <-
      "MCMCpack Factor Analysis Posterior Density Sample"
    return(output)
    
  }
