################################
## Poisson Changepoint Model  ##
################################

## This version does not allow covariates

## Noninformative prior is not acceptable 
## in case of marginal likelihood computation

## Originally written 12/23/2005 Jong Hee Park
## Modified 02/13/2006 Jong Hee Park 
## Modified 10/21/2006 Jong Hee Park
## Modified 12/13/2006 Jong Hee Park for JStatSoft article

## some utility functions added
## Jong Hee Park 07/14/2006

## two plot functions added
## Jong Hee Park 12/13/2006

##############################################################
## Andrew and Kevin, followings are a list of things to do
##############################################################

## 1. When you write a C code for MCMCpoissonChangepoint, 
##    you need to export "prob.state" in addition to mcmc output. 
##    This "prob.state" contains posterior probabilities of each state
##    and is essential to draw plots and posterior inference. 

## 2. I provide two special plot functions for changepoint models:
##    "plot.post.state" and "plot.post.changepoint." 
##    These plot functions need to be changed based on the changes in 
##    MCMCpoissonChangepoint() output.

## 3. We need to change outputs of MCMCpoissonChangepoint make it accessible 
##    by BayesFactor().

## 4. All helper functions should be straightforward. "rdirichlet.cp" and 
##    "ddirichlet.cp" might be redundant because they're simple a beta 
##    distribution. In the future, I will replace these into beta functions
##    but for the time being, use these functions. Assigning zeros in the right 
##    place is trickier than I thought. (I should have written these using 
##    beta density in the first place...)


 
#########################################################
"MCMCpoissonChangepoint"<-
#########################################################    
    function(data,  m = 1, burnin = 1000, mcmc = 1000, thin = 1, 
            verbose = 0, seed = NA, c0, d0, a = NULL, b = NULL,   
            marginal.likelihood = c("none", "Chib95"), ...) {

    ####################################################
    ## Arguments
    ## m : the number of changepoint
    ## c0 and d0: gamma prior for lambda. NO DEFAULT values 
    ## a and b: beta prior for transition probabilities
    ## By default, the expected duration is computed and 
    ## corresponding a and b values are assigned. The expected
    ## duration is the sample period divided by the number of states.
    ## ex) 300 years and 2 changepoints (3 states)
    ##        expected.duration <- 300/(2+1)  
    ##        b <- 0.1
    ##        a <- b*expected.duration
    ####################################################
    
    # check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin   
    cl <- match.call()
   
    ## ## seeds
    ## seeds <- form.seeds(seed) 
    ## lecuyer <- seeds[[1]]
    ## seed.array <- seeds[[2]]
    ## lecuyer.stream <- seeds[[3]]
    if(!is.na(seed)) set.seed(seed)
    
    ## sample size
    y <- data
    n <- length(y)
    
    ## prior 
    A0 <- trans.mat.prior(m=m, n=n, a=a, b=b) 
    
    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    
    ## storage matrices 
    lambda.store <-  matrix(NA, mcmc/thin, m+1)
    P.store     <-  matrix(NA, mcmc/thin, (m+1)^2)
    ps.store    <-  matrix(0, n, m+1)
    s1.store    <-  matrix(NA, mcmc/thin, n)
    py          <-  rep(0, m+1)
    pdf.P.store <-  matrix(NA, mcmc, m+1)
    lambda1      <-  rep(NA, m+1)
    P1          <-  matrix(NA, m+1, m+1)  
   
    ## starting values
    lambda0      <-  runif(m+1)
    P0          <-  trans.mat.prior(m=m, n=n, a=0.9, b=0.1)
   
    ####################################################################    
    ## MCMC iteration starts!  
    ####################################################################
   
    for (iter in 1:totiter){         
        
    ####################################################################
    ## Step 1: Sampling S
    ####################################################################   
    state.out <- Poisson.state.sampler(m=m, y=y, lambda=lambda0, P=P0)
    s1 <- state.out$s1
    ps1<- state.out$ps1
   
    ####################################################################
    ## Step 2: Sampling lambda 
    ####################################################################
    for (j in 1:(m+1)){
        ej  <-  as.numeric(s1==j)
        yj  <-  y[ej==1]
        nj  <-  length(yj)
        c1  <-  sum(yj) + c0
        d1  <-  nj + d0       
        lambda1[j]   <- rgamma(1, c1, d1)    
    }
    
    ####################################################################    
    ## Step 3: Sampling P
    ####################################################################
    switch  <-  switchg(s1) 
    for (j in 1:(m+1)){
        switch1 <-  A0[j,] + switch[j,]        
        pj      <-  rdirichlet.cp(1, switch1)
        P1[j,]  <-  pj
    }
    
    ## update
    lambda0  <-  lambda1    
    P0      <-  P1
    
    ####################################################################    
    ## end of one iteration
    ####################################################################
    
    ## store
    if (iter > burnin && (iter %% thin == 0)) {
        lambda.store[iter-burnin,]   <-  lambda1    
        P.store[iter-burnin,]       <-  as.vector(t(P1))
        s1.store[iter-burnin,]      <-  s1
        ps.store                    <-  ps.store + ps1
    }                 
    
    ## report 
    if(verbose > 0 && iter%%verbose == 0){
        cat("----------------------------------------------",'\n')
        cat("iteration = ", iter, '\n')
        cat("lambda = ", lambda1, '\n') 
        cat("Transition Matrix", '\n')
        for(i in 1:(m+1)) 
        cat(paste("", P1[i,]), fill=TRUE, labels=paste("{",i,"}:", sep=""), sep=",")
        }
        
   } ## end of MCMC sampling 
    
    ## marginal likelihood calculation if Chib == TRUE
    if (marginal.likelihood == "Chib95"){

    ############################################
    ## Bayes Factor Calculation Starts        ##
    ############################################
    lambda.st    <-  apply(lambda.store, 2, mean)
    P.vec.st    <-  apply(P.store, 2, mean)
    P.st        <-  t(matrix(P.vec.st, m+1, m+1))

    #################################################################
    ## 1. pdf.lambda
    #################################################################
    density.lambda <- matrix(NA, mcmc, m+1)
    for (i in 1:mcmc){
            for (j in 1:(m+1)){   
            # compute proposal density
            ej  <-  as.numeric(s1.store[i,]==j)
            yj  <-  y[ej==1]
            nj  <-  length(yj)
            c1  <-  sum(yj) + c0
            d1  <-  nj + d0          
            density.lambda[i, j]  <-  dgamma(lambda.st[j], c1, d1)     
            }
    }
    pdf.lambda   <-  log(prod(apply(density.lambda, 2, mean)))

    ######################################################################
    ## 2. pdf.P
    ######################################################################
    for (g in 1:mcmc){
        state.out   <-  Poisson.state.sampler(m=m, y=y, lambda=lambda.st, P=P0)
        s1          <-  state.out$s1
        ps1         <-  state.out$ps1       
        switch      <-  switchg(s1)         
    
        for (j in 1:(m+1)){
            switch1 <-  A0[j,] + switch[j,]         
            pj      <-  rdirichlet.cp(1, switch1)
            P1[j,]  <-  pj
            pdf.P.store[g,j]<- ddirichlet.cp(P.st[j,], switch1)
        }
        P0  <- P1
    }
    pdf.P   <-  log(prod(apply(pdf.P.store, 2, mean)))
    
    ####################################################################
    ## likelihood
    ####################################################################
    F   <-  matrix(NA, n, m+1)      
    like<-  rep(NA, n)
    pr1 <-  c(1,rep(0,m))           
    for (t in 1:n){
        py  <-  sapply(c(1:(m+1)), function(i){poisson.pdf(y[t], lambda.st[i])})
        if(t==1) {pstyt1 = pr1}                                           
        else {pstyt1 <- F[t-1,]%*%P.st} 
        unnorm.pstyt <- pstyt1*py       
        pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) 
        F[t,]   <-  pstyt
        like[t] <-  sum(unnorm.pstyt)
    }
    loglik  <-  sum(log(like))
    
    ####################################################################
    ## prior ordinates
    ####################################################################
    nprior.lambda<-  nprior.P <- rep(NA, m+1)
    nprior.lambda<-  sapply(c(1:(m+1)), function(i){dgamma(lambda.st[i], c0, d0, log=TRUE)})
    nprior.P    <-  sapply(c(1:(m+1)), function(i){log(ddirichlet.cp(P.st[i,], A0[i,]))})
    
    prior.lambda <-  sum(nprior.lambda)
    prior.P     <-  sum(nprior.P)
    
    ############################################
    ## Marginal Likelihood                    ##
    ############################################
    numerator   <-  loglik + prior.lambda + prior.P
    denominator <-  pdf.lambda + pdf.P
    logmarglike <-  numerator - denominator
    
    ## print    
    if(verbose > 0){
    cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    cat("Log Marginal Likelihood\n")
    cat("-------------------------------------------------",'\n')
    cat("log(marglike)= ",      logmarglike, '\n')
    cat("log(likelihood)= ",    loglik, '\n')
    cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    }
    } ## end of marginal likelihood calculation
        
    else {marginal <- loglik <- NULL}
    
    ##  return output
    output  <-  as.mcmc(lambda.store)
    varnames(output) <- paste("lambda.", 1:(m+1), sep = "")
    attr(output, "title") <- "MCMCpoissonChangepoint Posterior Sample"
    attr(output, "y") <- data
    attr(output, "call") <- cl
    attr(output, "logmarglike") <- logmarglike
    attr(output, "loglik") <- loglik
    attr(output, "prob.state") <- ps.store/mcmc
    return(output)
    
 }## end of MCMC function

