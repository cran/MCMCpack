#########################################################
## Internvetion Analysis using Changepoint Model
#########################################################

"MCMCintervention"<-
  function(y, data=parent.frame(),  m = 1, intervention = 1,
           prediction.type=c("trend", "ar"),
           change.type = c("fixed", "random", "all"),
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, sigma.mu = NA, sigma.var = NA,
           a = NULL, b = NULL,
           mcmc = 1000, burnin = 1000,  thin = 1, verbose = 0, 
           seed = NA, beta.start = NA, P.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...){
    
    ## form response and model matrices
    y <- as.vector(y)
    n <- length(y)
    ns <- m + 1 # number of states
    if (prediction.type == "trend"){
      X <- matrix(cbind(1, c(1:n)), n, 2)
      xnames <- c("constant", "trend") 
    }
    else if (prediction.type == "ar"){
      y1 <- y[1:(n/2)]
      ar1 <- arima(y1, c(1,0,0))
      muy0 <- as.numeric(ar1$coef[2])
      sigmay0 <- sqrt(as.numeric(as.numeric(ar1[2])/(1 - ar1$coef[1]^2)))
      y0 <- rnorm(1, muy0, sigmay0)
      X <- matrix(cbind(1, c(y0, y[-n])), n, 2)
      xnames <- c("constant", "lag.y") 
    }
    else {
      X <- matrix(1, n, 1)
      xnames <- c("constant") 
    }
    k <- ncol(X)    # number of covariates
    
    ## check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    nstore <- mcmc/thin    
    cl <- match.call()

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    if(!is.na(seed)) set.seed(seed)
  
    ## prior
    mvn.prior <- form.mvn.prior(b0, B0, k)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    if (prediction.type == "ar"){
      if (b0[2] > 1|b0[2] < -1){
        stop("The prior of AR coefficient ",b0[2], " is outside the stationary region! \n") 
      }
    }
    if (is.na(sigma.mu)|is.na(sigma.var)) {
      check.ig.prior(c0, d0)
    }
    else {
      d0 <- 2*(sigma.mu + sigma.mu^3/sigma.var)
      c0 <- 2*(1 + (d0/2)/sigma.mu)
    }
  
    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    
    ## following MCMCregress, set chib as binary
    logmarglike <- loglik <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }
    
    ## initial values
    Y <- matrix(y, n, 1)
    if (m == 0){
      output <- MCMCregress(formula = Y ~ X-1, mcmc=mcmc,
                            burnin=burnin, verbose=verbose, thin=thin,
                            b0 = b0, B0 = solve(B0), c0 = c0, d0 = d0, 
                            marginal.likelihood = marginal.likelihood)
      attr(output, "y")       <- y
      attr(output, "intervention")  <- intervention
      yhatout <- output[, 1:2]%*%t(X)
      attr(output, "yhat") <- matrix(yhatout, nstore, n) ## X_j*beta_j
      attr(output, "yforepred") <- matrix(yhatout, nstore, n) ## yt|Y_t-1, theta
      attr(output, "ybackpred") <- matrix(yhatout, nstore, n) ## yt|Y_t+1, theta

    }
    else{ ## if m > 0
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)  
      Pstart  <-  check.P(P.start, m, a=a, b=b)
      ols <- lm(y~X-1)
      bols <- coef(ols)
      betastart  <- matrix(rep(bols, ns), ns, k, byrow = TRUE)
      Sigmastart <- rep(summary(ols)$sigma^2, ns)
      statestart <- sort(sample(1:ns, n, replace=T))
      AR <- 0
      if (prediction.type == "ar"){
        AR <- 1
      }
      change <- 0
      betaout.N <- nstore*ns*k
      Sigmaout.N <- nstore*ns
      if (change.type == "fixed"){
        change <- 1
        betaout.N <- nstore*ns*k
        Sigmaout.N <- nstore
     }
      if (change.type == "random"){
        change <- 2
        betaout.N <- nstore*k
        Sigmaout.N <- nstore*ns
      }
      
      ## call C++ code to draw sample
      posterior <- .C("MCMCintervention",
                      accept = as.double(0.0),
                      betaout = as.double(rep(0.0, betaout.N)), 
                      Sigmaout = as.double(rep(0.0, Sigmaout.N)), 
                      Pout = as.double(rep(0.0, nstore*ns*ns)), 
                      psout = as.double(rep(0.0, n*ns)),
                      sout = as.double(rep(0.0, nstore*n)),
                      
                      yhatout = as.double(rep(0.0, nstore*n)),
                      yerrorout = as.double(rep(0.0, nstore*n)),
                      yforepredout = as.double(rep(0.0, nstore*n)),
                      ybackpredout = as.double(rep(0.0, nstore*n)),
                      
                      Ydata = as.double(Y),
                      Yrow = as.integer(nrow(Y)),
                      Ycol = as.integer(ncol(Y)),
                      Xdata = as.double(X),
                      Xrow = as.integer(nrow(X)),
                      Xcol = as.integer(ncol(X)),
                      
                      m = as.integer(m),
                      intervention = as.integer(intervention), 
                      burnin = as.integer(burnin),           
                      mcmc = as.integer(mcmc), 
                      thin = as.integer(thin),
                      verbose = as.integer(verbose),
                      
                      lecuyer=as.integer(lecuyer), 
                      seedarray=as.integer(seed.array),
                      lecuyerstream=as.integer(lecuyer.stream),
                      
                      betastart = as.double(betastart),  
                      Sigmastart = as.double(Sigmastart),  
                      Pstart = as.double(Pstart),
                      statestart = as.integer(statestart),    
                      
                      a = as.double(a),
                      b = as.double(b),
                      b0data = as.double(b0),
                      B0data = as.double(B0), 
                      c0 = as.double(c0),
                      d0 = as.double(d0),
                      A0data = as.double(A0), 
                      logmarglikeholder = as.double(0.0),
                      loglikeholder = as.double(0.0),
                      ar = as.integer(AR),
                      change = as.integer(change), 
                      chib = as.integer(chib))                
      
      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
        loglike <- posterior$loglikeholder
      }
      ## pull together matrix and build MCMC object to return
      beta.holder <- matrix(posterior$betaout, nstore, )
      Sigma.holder <- matrix(posterior$Sigmaout, nstore, )
      P.holder    <- matrix(posterior$Pout, nstore, )
      s.holder    <- matrix(posterior$sout, nstore, )
      ps.holder   <- matrix(posterior$psout, n, )
      
      output1 <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      if (change == 0){
        varnames(output1)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c(xnames), "_regime", i, sep = "")
                                     })
        output2 <- mcmc(data=Sigma.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
        varnames(output2)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c("sigma2"), "_regime", i, sep = "")
                                     }) 
      }
      else if (change == 1){
        varnames(output1)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c(xnames), "_regime", i, sep = "")
                                     })
        output2 <- mcmc(data=Sigma.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
        names(output2)  <- c("sigma2") 
      }
      else{
        output2 <- mcmc(data=Sigma.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
        varnames(output2)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c("sigma2"), "_regime", i, sep = "")
                                     }) 
      }
      ## To check the acceptance rate (if AR == 1)
      accept <- posterior$accept

      output <- as.mcmc(cbind(output1, output2))
      attr(output, "title") <- "MCMCintervention Posterior Sample"
      attr(output, "intervention")  <- intervention
      attr(output, "accept")  <- accept
      attr(output, "y")       <- y
      attr(output, "X")       <- X
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglik") <- loglik
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "s.store") <- s.holder
      attr(output, "yhat") <- matrix(posterior$yhatout, nstore, n)## X_j*beta_j
      attr(output, "yerror") <- matrix(posterior$yerrorout, nstore, n)## y_j - X_j*beta_j
      attr(output, "yforepred") <- matrix(posterior$yforepredout, nstore, n)## yt|Y_t-1, theta
      attr(output, "ybackpred") <- matrix(posterior$ybackpredout, nstore, n)## yt|Y_t+1, theta
    }
    return(output)
    
}## end of MCMC function

