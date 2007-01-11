########## Utility Functions for Bayesian Times Series Models ##########

## written and maintained by:
##    Jong Hee Park
##    Department of Political Science
##    University of Chicago
##    jhp@uchicago.edu

## NOTE: only the plot functions are documented and exported in the
## NAMESPACE.

##############################################################
## Helper functions for MCMCPoissonChangepoint()
##############################################################
    
    "switchg" <-  function(s1){
    
    ## "switchg" computes the frequency of elements in state vector
    ## for one step ahead Markov process only
    ## example: s1 <- c(1,1,1,2,2,2,3,4,4,4,4,5,5)
    ## switchg(s1)
    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    2    1    0    0    0
    ## [2,]    0    2    1    0    0
    ## [3,]    0    0    0    1    0
    ## [4,]    0    0    0    3    1
    ## [5,]    0    0    0    0    1

            s <- max(s1)
            out <- matrix(0,s,s)

            ## all P(i,i+1) are 1
            for (i in 1:(s-1)){
                out[i,i+1] <- 1}

            ## diagonal elements is (the number of occurrence - 1)
            diag(out) <- table(s1)-1
            return(out)
    }
    
    ## "trans.mat.prior" makes a transition matrix
    
    ## a and b are beta prior for transition matrix
    ## By default, the expected duration is computed and 
    ## corresponding a and b values are assigned. The expected
    ## duration is the sample period divided by the number of states.
    ## ex) 300 years and 2 changepoints (3 states)
    ##        expected.duration <- 300/(2+1)  
    ##        b <- 0.1
    ##        a <- b*expected.duration
 
    "trans.mat.prior" <- function(m, n, a=NULL, b=NULL){
        if (!is.null(a)|!is.null(b)){
            a <- a
            b <- b
        }
        else {expected.duration <- round(n/(m+1))
            b <- 0.1
            a <- b*expected.duration
        }
        trans <- diag(1, m+1)
        # put a as diagonal elements except the last row
        diag(trans)[1:m]<-rep(a, m)
        # put b in trans[i, i+1]
        for (i in 1:m){trans[i, i+1]<-b}
        return(trans)
    } 
    
    "rdirichlet.cp" <- function(n, alpha){
    ## "rdirichlet.cp" picks n random deviates from the Dirichlet function
    ## with shape parameters alpha
    ## Note that alpha can contain zero to deal with transition matrix rowwise.
    ## It returns a matrix.

        col <-  length(alpha)
        out <-  matrix(NA, n, col)
        out[,which(alpha==0)]<- 0
        a   <-  alpha[alpha!=0]
        l   <-  length(a);
        x   <-  matrix(rgamma(l*n,a), ncol=l, byrow=TRUE);
        sm  <-  x%*%rep(1,l);
        dir.out <-  x/as.vector(sm);
        ## combine zero and nonzero parts prob
        out[,which(alpha!=0)] <- dir.out
        return(out)
    }

    "ddirichlet.cp" <- function (x, alpha){
    ## "ddirichlet.cp" is a density function for the Dirichlet distribution
    ## alpha and x can contain zero
    ## but zero should be in the same place
    ## x=c(0.9, 0.1, 0) and alpha=c(6,1,0)
    ## It returns a matrix
        col <- length(alpha)
        out <- rep(NA, col)
        out[which(alpha==0)] <- 0
        a <- alpha[alpha!=0]
        x1 <- x[x!=0]
        dirichlet1 <- function(x1, a) {
            logD <- sum(lgamma(a)) - lgamma(sum(a))
            s <- sum((a - 1) * log(x1))
            exp(sum(s) - logD)
        }
        if (!is.matrix(x1))
            if (is.data.frame(x1))
                x1 <- as.matrix(x1)
            else x1 <- t(x1)
        if (!is.matrix(a))
            a <- matrix(a, ncol = length(a), nrow = nrow(x1),
                byrow = TRUE)
        if (any(dim(x1) != dim(a)))
            stop("Mismatch between dimensions of x and alpha in
            ddirichlet().\n")
        pd <- vector(length = nrow(x1))
        for (i in 1:nrow(x1)) pd[i] <- dirichlet1(x1[i, ], a[i,])
        pd[apply(x1, 1, function(z) any(z < 0 | z > 1))] <- 0
        pd[apply(x1, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
        return(pd)
    }
    
    "poisson.pdf" <- function(y, lambda){
    ## What this function does is same with "prod(dpois(y, lambda))"
        n        <-  length(y)       
        log.like <-  sum(sapply(c(1:n), function(i){-lambda + y[i]*log(lambda)- log(factorial(y[i]))}))
        return(exp(log.like))
    }

    "Poisson.state.sampler" <- function(m, y, lambda, P){
    ## "Poisson.state.sampler" samples state vector (s1) by forward and backward sampling
    ## P is a (m+1 by m+1) transition matrix
    ## F matrix that contains all the information of Pr(st|Yt)
        
        ## Forward sampling: update F matrix
        n   <-  length(y)
        F   <-  matrix(NA, n, m+1)     # storage for the Filtered probabilities       
        pr1 <-  c(1,rep(0, m))         # initial probability Pr(s1=k|Y0, lambda)       
        for (t in 1:n){ 
            py  <-  sapply(c(1:(m+1)), function(i){poisson.pdf(y[t], lambda[i])})
            if(t==1) pstyt1 = pr1           
            else pstyt1     <- F[t-1,]%*%P   
            unnorm.pstyt    <- pstyt1*py      
            pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
            F[t,]   <-  pstyt
        }
    
        ## Backward recursions: sample s1 using F matrix and a transition matrix (P)
        s1      <-  matrix(NA, n, 1)   ## holder for state variables   
        ps1     <-  matrix(NA, n, m+1) ## holder for state probabilities 
        ps1[n,] <-  F[n,]              ## we know last elements of ps1 and s1
        s1[n,1] <-  m+1                 
        t       <-  n-1
    
        while (t>=1){
            st1     <-  s1[t+1]                         
            unnorm.pstyn   <-  F[t,]*P[,st1]      
            ## normalize into a prob. density       
            pstyn   <-  unnorm.pstyn/sum(unnorm.pstyn) 
            if (st1==1) {s1[t]<-1}
            ## If this is not the first period,
            ## draw a state variable from a discrete prob. distribution("pstyn")
            ## using the inverse CDF method.  
            
            ## What inverse CDF does is as follows
            ## In backward recursion, the current state (t) is always at st1, which is one of (1...m+1)
            ## The probability that transition happens at t-1, that is Pr(t-1=st1-1|t=st1) can be found 
            ## by using (dicrete case) inverse CDF method for the probability vector we have (pstyn).
            ## draw u ~ U[0,1] then, if u<pstyn[st1-1], t-1=st1-1, otherwise t-1=st1
            ## This is because Pr(x<t)=Pr(F^-1(u)<t)=Pr(u<F(t))=F(t) 
            ## Thus, when u ~ U[0,1], x=F^-1(u) ~ F

            else {
                pone    <-  pstyn[st1-1]                        
                s1[t]   <-  ifelse (runif(1) < pone, st1-1, st1)} 
                ps1[t,] <-  pstyn  # probabilities pertaining to a certain state                            
                t       <-  t-1
        }
        ## name and report outputs 
        out <-  list(s1, ps1)
        names(out)<-c("s1","ps1")
    return(out)
    }
    
 
##############################################################
##  Plot functions for MCMCPoissonChangepoint()
##############################################################

    ## "plotPostState" draws a plot of posterior distribution of states
 
    "plotPostState" <- 
    function(mcmcout, y, m, legend.control=NULL, cex=0.8, lwd=1.2){
    
    ## To use legend.control, a user has to specify xy coordinate for legend.control.
    ## ex) legend.control = c(x, y)
    
    ## extract posterior state probability vector        
    out <- attr(mcmcout, "prob.state")
       
    ## y must be a ts object. Otherwise, it is coerced as a ts object.
    if (!is.ts(y)) y <- ts(y)
    time.frame <- as.vector(time(y))     
    
    ## draw a plot for each state variable
        plot(1, 0, xlim=range(time.frame), ylim=c(0,1), type="n", main="", xlab="Time", cex=cex, lwd=lwd,
            ylab=expression(paste("Pr(", S[t], "= k |", Y[t], ")")))
            for (i in 1:(m+1))
            points(time.frame, out[,i], type="o", lty=i, lwd=lwd, col=i, cex=cex) 
            
    ## put a legend to indicate each line for the state 
        if (!is.null(legend.control)){ 
            if (length(legend.control)!=2)
            stop("You should specify x and y coordinate for a legend.")
            else
            legend(legend.control[1], legend.control[2], legend = paste("State", 1:(m+1), sep=""), 
                col=1:(m+1), lty=1:(m+1), lwd = rep(lwd, m+1), pch=rep(1, m+1), bty="n")
        }  
    }
    
    ## "plotPostChangepoint" draws a plot of posterior changepoint probability
    "plotPostChangepoint" <-
    function(mcmcout, y, m, xlab="Time", ylab=""){

    ## extract posterior state probability vector
    out <- attr(mcmcout, "prob.state")

    ## y must be a ts object. Otherwise, it is coerced as a ts object.
    if (!is.ts(y)) y <- ts(y)
    time.frame <- as.vector(time(y))   

        ## divide space
        if (m == 1){
        pr.st <- c(0, diff(out[,(m+1)]))
        plot(time.frame, pr.st, type="h", main="", xlab=xlab, ylab=ylab)
        cp <- which(cumsum(pr.st) > 0.5)[1]-1
        abline(v=cp+ time.frame[1], lty=3, col="red")
        }
        else {
        par(mfrow=c(m,1));par(mar=c(2,4,1,1))
        cp <- rep(NA, m) # holder for expected changepoints
            for (i in 2:m){
                pr.st <- c(0, diff(out[,i]))
                pr.st <- ifelse(pr.st <0, 0, pr.st)
                plot(time.frame, pr.st, type="h", main="", xlab=xlab, ylab=ylab)
                cp[i-1] <- which(cumsum(pr.st)>0.5)[1]-1
                abline(v=cp[i-1]+ time.frame[1], lty=3, col="red")

            }
            pr.st <- c(0, diff(out[,(m+1)]))
            plot(time.frame, pr.st, type="h", main="", xlab=xlab, ylab=ylab)
            cp[m] <- which(cumsum(pr.st)>0.5)[1]-1
            abline(v=cp[m]+ time.frame[1], lty=3, col="red")
        }

    cp.means <- rep(NA, m+1)
    cp.start <- c(1, cp + 1)
    cp.end <- c(cp, length(y))
    cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    cat("Expected changepoint(s) ", cp + time.frame[1], '\n')
    for (i in 1:(m+1)) cp.means[i] <- mean(y[cp.start[i]:cp.end[i]])
    cat("Local means for each regime are ", cp.means, '\n')
    cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    }
