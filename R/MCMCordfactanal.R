##########################################################################
# sample from the posterior distribution of a factor analysis model
# model in R using linked C++ code in Scythe.
#
# The model is:
#
# x*_i = \Lambda \phi_i + \epsilon_i,   \epsilon_i \sim N(0, I)
#
# \lambda_{ij} \sim N(l0_{ij}, L0^{-1}_{ij})
# \phi_i \sim N(0,I)
#
# and x*_i is the latent variable formed from the observed ordinal
# variable in the usual (Albert and Chib, 1993) way.
#
# Andrew D. Martin
# Washington University
#
# Kevin M. Quinn
# Harvard University
#
# May 12, 2003
#
##########################################################################

"MCMCordfactanal" <-
  function(x, factors, lambda.constraints=list(),
           data=list(), burnin = 1000, mcmc = 10000,
           thin=5, tune=NA, verbose = FALSE, seed = 0,
           lambda.start = NA, l0=0, L0=0,
           store.lambda=TRUE, store.scores=FALSE, ... ) {

    # check for MCMCirtKd special case, this is used to tell the R
    # and C++ code what to echo (1 if regular, 2 if MCMCirtKd)
    # the test is based on the existence of a string
    # "special.case" coming through the ellipses
    args <- list(...)

    if(length(args) > 0) test.string <- args[[1]][1]
    else test.string <- "xxx"
        
    if(test.string=="special.case") {
      case.switch <- 2
      echo.name <- "MCMCirtKd"
    }
    else {
      case.switch <- 1
      echo.name <- "MCMCordfactanal"
    } 
 
    # extract X and variable names from the model formula and frame       
    if (is.matrix(x)){
      X <- as.data.frame(x)
      xvars <- dimnames(X)[[2]]
      xobs <- dimnames(X)[[1]]
      N <- nrow(X)    # number of observations
      K <- ncol(X)    # number of manifest variables
      ncat <- matrix(NA, K, 1) # vector of number of categ. in each man. var. 
      for (i in 1:K){
        X[,i] <- factor(X[,i], ordered=TRUE)
        ncat[i] <- nlevels(X[,i])
        X[,i] <- as.integer(X[,i])
        X[is.na(X[,i]), i] <- -999
      }
      X <- as.matrix(X)
    }
    else {
      call <- match.call()
      mt <- terms(x, data=data)
      if (attr(mt, "response") > 0) 
        stop("Response not allowed in formula in ", echo.name, "().\n")
      if(missing(data)) data <- sys.frame(sys.parent())
      mf <- match.call(expand.dots = FALSE)
      mf$factors <- mf$lambda.constraints <- mf$burnin <- mf$mcmc <- NULL
      mf$thin <- mf$tune <- mf$verbose <- mf$seed <- NULL
      mf$lambda.start <- mf$l0 <- mf$L0 <- NULL
      mf$store.lambda <- mf$store.scores <- mf$... <- NULL
      mf$drop.unused.levels <- TRUE
      mf[[1]] <- as.name("model.frame")
      mf <- eval(mf, sys.frame(sys.parent()))
      attributes(mt)$intercept <- 0
      Xterm.length <- length(attr(mt, "variables"))
      X <- subset(mf,
                  select=as.character(attr(mt, "variables"))[2:Xterm.length])
      N <- nrow(X)	        # number of observations      
      K <- ncol(X)              # number of manifest variables
      ncat <- matrix(NA, K, 1)  # vector of number of categ. in each man. var. 
      for (i in 1:K){
        X[,i] <- factor(X[,i], ordered=TRUE)
        ncat[i] <- nlevels(X[,i])      
        X[,i] <- as.integer(X[,i])
        X[is.na(X[,i]), i] <- -999
      }
      X <- as.matrix(X)
      xvars <- dimnames(X)[[2]] # X variable names
      xobs <- dimnames(X)[[1]]  # observation names
    }
    
    if (is.null(xobs)){
      xobs <- 1:N
    }
    
    check.parameters(burnin, mcmc, thin, echo.name)
    
    # give names to the rows of Lambda related matrices
    Lambda.eq.constraints <- matrix(NA, K, factors+1)
    Lambda.ineq.constraints <- matrix(0, K, factors+1)
    Lambda.prior.mean <- matrix(0, K, factors+1)
    Lambda.prior.prec <- matrix(0, K, factors+1)
    
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
        if(case.switch==1) {
          cat("Constraints on Lambda are logically inconsistent.\n")
        }
        else {
          cat("Constraints on item parameters are logically inconsistent.\n")
        }
        stop("Please respecify and call ", echo.name, "() again\n")
      }
    }
    Lambda.eq.constraints[is.na(Lambda.eq.constraints)] <- -999
    
    # setup prior means and precisions for Lambda
    # prior means
    if (is.matrix(l0)){ # matrix input for l0
      if (nrow(l0)==K && ncol(l0)==(factors+1))
        Lambda.prior.mean <- l0
      else {
        if(case.switch==1) {
          cat("l0 not of correct size for model specification.\n")
        }
        else {
          cat("b0 not of correct size for model specification.\n")
        }
        stop("Please respecify and call ", echo.name, "() again\n")
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
      Lambda.prior.mean <- matrix(l0, K, factors+1)
    }
    else {
      if(case.switch==1) {    
        cat("l0 neither matrix, list, nor scalar.\n")
      }
      else {
        cat("b0 neither matrix, list, nor scalar.\n")      
      }
      stop("Please respecify and call ", echo.name, "() again\n")
    }
    
    # prior precisions
    if (is.matrix(L0)){ # matrix input for L0
      if (nrow(L0)==K && ncol(L0)==(factors+1))
        Lambda.prior.prec <- L0
      else {
        if(case.switch==1) {       
          cat("L0 not of correct size for model specification.\n")
        }
        else {
          cat("B0 not of correct size for model specification.\n")      
        }
        stop("Please respecify and call ", echo.name, "() again\n")
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
      Lambda.prior.prec <- matrix(L0, K, factors+1)
    }
    else {
      if(case.switch==1) {  
        cat("L0 neither matrix, list, nor scalar.\n")
      }
      else {
        cat("B0 neither matrix, list, nor scalar.\n")      
      }
      stop("Please respecify and call ", echo.name, "() again\n")
    }
    if (min(L0) < 0){
      if(case.switch==1) {      
        cat("L0 contains negative elements.\n")
      }
      else {
        cat("B0 contains negative elements.\n")      
      }
      stop("Please respecify and call ", echo.name, "() again\n")
    }
    
    # Starting values for Lambda
    Lambda <- matrix(0, K, factors+1)
    if (is.na(lambda.start)){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              if (j==1){
                if (ncat[i] == 2){
                  probit.out <- glm(as.factor(X[X[,i]!=-999,i])~1,
                                    family=binomial(link=probit))
                  probit.beta <- coef(probit.out)
                  Lambda[i,j] <- probit.beta[1]
                }
                else{
                  polr.out <- polr(ordered(X[X[,i]!=-999,i])~1)
                  Lambda[i,j] <- -polr.out$zeta[1]*.588
                }
              }
              else Lambda[i,j] <- 0
            }
            if(Lambda.ineq.constraints[i,j]>0){
              Lambda[i,j] <- 1.0
            }
            if(Lambda.ineq.constraints[i,j]<0){
              Lambda[i,j] <- -1.0
            }          
          }
          else Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else if (is.matrix(lambda.start)){
      if (nrow(lambda.start)==K && ncol(lambda.start)==(factors+1))
        Lambda  <- lambda.start
      else {
        cat("Starting values not of correct size for model specification.\n")
        stop("Please respecify and call ", echo.name, "() again\n")
      }
    }
    else if (length(lambda.start)==1 && is.numeric(lambda.start)){
      Lambda <- matrix(lambda.start, K, factors+1)
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else {
      cat("Starting values neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call ", echo.name, "() again\n")
    }
    
    # check MH tuning parameter
    if (is.na(tune)){
      tune <- matrix(NA, K, 1)
      for (i in 1:K){
        tune[i] <- 0.05/ncat[i]
      }
    }
    else if (is.double(tune)){
      tune <- matrix(tune, K, 1)
    }
    if(min(tune) < 0) {
      cat("Tuning parameter is negative.\n")
      stop("Please respecify and call ", echo.name, "() again\n")
    }
  
    # starting values for gamma (note: not changeable by user)
    gamma <- matrix(0, max(ncat)+1, K)
    for (i in 1:K){
      if (ncat[i]==2){
        gamma[1,i] <- -300
        gamma[2,i] <- 0
        gamma[3,i] <- 300
      }
      else {
        polr.out <- polr(ordered(X[X[,i]!=-999,i])~1)
        gamma[1,i] <- -300
        gamma[2,i] <- 0
        gamma[3:ncat[i],i] <- (polr.out$zeta[2:(ncat[i]-1)] -
                               polr.out$zeta[1])*.588
        gamma[ncat[i]+1,i] <- 300
      }
    }
    
    # define holder for posterior density sample
    if (store.scores == FALSE && store.lambda == FALSE){
      sample <- matrix(data=0, mcmc/thin, length(gamma))
    }
    else if (store.scores == TRUE && store.lambda == FALSE){
      sample <- matrix(data=0, mcmc/thin, (factors+1)*N + length(gamma))
    }
    else if(store.scores == FALSE && store.lambda == TRUE) {
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+length(gamma))
    }
    else { # store.scores==TRUE && store.lambda==TRUE
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+(factors+1)*N +
                       length(gamma))
    }
    
    # Call the C++ code to do the real work
    posterior <- .C("ordfactanalpost",
                    samdata = as.double(sample),
                    samrow = as.integer(nrow(sample)),
                    samcol = as.integer(ncol(sample)),
                    X = as.integer(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    tune = as.double(tune),
                    seed = as.integer(seed),
                    verbose = as.integer(verbose),
                    Lambda = as.double(Lambda),
                    Lambdarow = as.integer(nrow(Lambda)),
                    Lambdacol = as.integer(ncol(Lambda)),
                    gamma = as.double(gamma),
                    gammarow = as.integer(nrow(gamma)),
                    gammacol = as.integer(ncol(gamma)),
                    ncat = as.integer(ncat),
                    ncatrow = as.integer(nrow(ncat)),
                    ncatcol = as.integer(ncol(ncat)),
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
                    storelambda = as.integer(store.lambda),
                    storescores = as.integer(store.scores),
                    accepts = as.integer(0),
                    outswitch = as.integer(case.switch),
                    PACKAGE="MCMCpack"
                    )
    if(case.switch==1) {
      cat(" overall acceptance rate = ",
          posterior$accepts / ((posterior$burnin+posterior$mcmc)*K), "\n")
    }
    
    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$samdata, posterior$samrow, posterior$samcol,
                     byrow=TRUE)
    output <- mcmc2(data=sample,start=1, end=mcmc, thin=thin)
    
    par.names <- NULL
    if (store.lambda==TRUE){
      if(case.switch==1) {
      Lambda.names <- paste(paste("Lambda",
                                  rep(X.names,
                                      each=(factors+1)), sep=""),
                            rep(1:(factors+1),K), sep="_")
      }
      if(case.switch==2) {
        alpha.hold <- paste("alpha", 1:K, sep = "")
        beta.hold <- paste("beta", 1:K, sep = "")
        beta.hold <- rep(beta.hold, factors, each=factors)
        beta.hold <- paste(beta.hold, 1:factors, sep="_")
  
        Lambda.names <- t(cbind(matrix(alpha.hold, K, 1), 
                                matrix(beta.hold,K,factors,byrow=TRUE)))  
        dim(Lambda.names) <- NULL
      }
      par.names <- c(par.names, Lambda.names)
    }
    
    gamma.names <- paste(paste("gamma",
                               rep(0:(nrow(gamma)-1),
                                   each=K), sep=""),
                         rep(X.names,  nrow(gamma)), sep="_")
    par.names <- c(par.names, gamma.names)
    
    if (store.scores==TRUE){
      if(case.switch==1) {
      phi.names <- paste(paste("phi",
                               rep(xobs, each=(factors+1)), sep="_"),
                         rep(1:(factors+1),(factors+1)), sep="_")
      par.names <- c(par.names, phi.names)
      }
      if(case.switch==2) {
      phi.names <- paste(paste("theta",
                               rep(xobs, each=(factors+1)), sep="_"),
                         rep(0:factors,(factors+1)), sep="_")
      par.names <- c(par.names, phi.names)      
      
      }
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
    attr(output, "accept.rate") <- posterior$accepts /
      (posterior$burnin+posterior$mcmc*K)
    if(case.switch==1) {
      attr(output,"title") <-
        "MCMCpack Ordinal Data Factor Analysis Posterior Density Sample"
    }
    if(case.switch==2) {
      attr(output,"title") <-
        "MCMCpack K-Dimensional Item Response Theory Model Posterior Density Sample"
    }
    return(output)
    
  }

