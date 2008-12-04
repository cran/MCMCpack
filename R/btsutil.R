##########################################################################
## Utility Functions for Bayesian Times Series Models
##
## written and maintained by:
##    Jong Hee Park
##    Department of Political Science
##    University of Chicago
##    jhp@uchicago.edu
##
## Revised on 09/12/2007 JHP	  
##
## NOTE: only the plot functions are documented and exported in the
## NAMESPACE.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################



##############################################################
## Helper functions for MCMCPoissonChangepoint()
##############################################################
    
	## switch a state vector into a matrix containing the number of states
    "switchg" <-  function(s1){    
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
    
   ## "plotState" draws a plot of posterior distribution of states 
    "plotState" <-
    function (mcmcout, legend.control = NULL, cex = 0.8, lwd = 1.2)
    {
    out <- attr(mcmcout, "prob.state")
    y <- attr(mcmcout, "y")
    m <- attr(mcmcout, "m")
    if (!is.ts(y))
        y <- ts(y)
    time.frame <- as.vector(time(y))
    plot(1, 0, xlim = range(time.frame), ylim = c(0, 1), type = "n",
        main = "", xlab = "Time", cex = cex, lwd = lwd, ylab = expression(paste("Pr(",
            S[t], "= k |", Y[t], ")")))
    for (i in 1:(m + 1)) points(time.frame, out[, i], type = "o",
        lty = i, lwd = lwd, col = i, cex = cex)
    if (!is.null(legend.control)) {
        if (length(legend.control) != 2)
            stop("You should specify x and y coordinate for a legend.")
        else legend(legend.control[1], legend.control[2], legend = paste("State",
            1:(m + 1), sep = ""), col = 1:(m + 1), lty = 1:(m +
            1), lwd = rep(lwd, m + 1), pch = rep(1, m + 1), bty = "n")
    }
}

    ## "plotChangepoint" draws a plot of posterior changepoint probability
    "plotChangepoint" <-
    function (mcmcout, xlab = "Time", ylab = "", verbose = FALSE)
    {
    out <- attr(mcmcout, "prob.state")
    y <- attr(mcmcout, "y")
    m <- attr(mcmcout, "m")
    if (!is.ts(y))
        y <- ts(y)
    time.frame <- as.vector(time(y))
    if (m == 1) {
        pr.st <- c(0, diff(out[, (m + 1)]))
        plot(time.frame, pr.st, type = "h", main = "", xlab = xlab,
            ylab = ylab)
        cp <- which(cumsum(pr.st) > 0.5)[1] - 1
        abline(v = cp + time.frame[1], lty = 3, col = "red")
    }
    else {
        par(mfrow = c(m, 1))
        par(mar = c(2, 4, 1, 1))
        cp <- rep(NA, m)
        for (i in 2:m) {
            pr.st <- c(0, diff(out[, i]))
            pr.st <- ifelse(pr.st < 0, 0, pr.st)
            plot(time.frame, pr.st, type = "h", main = "", xlab = xlab,
                ylab = ylab)
            cp[i - 1] <- which(cumsum(pr.st) > 0.5)[1] - 1
            abline(v = cp[i - 1] + time.frame[1], lty = 3, col = "red")
        }
        pr.st <- c(0, diff(out[, (m + 1)]))
        plot(time.frame, pr.st, type = "h", main = "", xlab = xlab,
            ylab = ylab)
        cp[m] <- which(cumsum(pr.st) > 0.5)[1] - 1
        abline(v = cp[m] + time.frame[1], lty = 3, col = "red")
    }
    cp.means <- rep(NA, m + 1)
    cp.start <- c(1, cp + 1)
    cp.end <- c(cp, length(y))
	if (verbose == TRUE){
		cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
		cat("Expected changepoint(s) ", cp + time.frame[1], "\n")
		for (i in 1:(m + 1)) cp.means[i] <- mean(y[cp.start[i]:cp.end[i]])
			cat("Local means for each regime are ", cp.means, "\n")
			cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
	}
}
    ## prior check for transition matrix
    "check.P" <-
    function(P.start = NA, m=m, n=n, a=a, b=b){
        if (is.na(P.start)[1]){
            P <- trans.mat.prior(m=m, n=n, a=0.9, b=0.1)}
        else if ((dim(P.start)[1]==m+1)&(dim(P.start)[2]==m+1)){
            if ((max(P.start)>1)||(min(P.start)<0)){
                stop("Error: P starting values are out of unit range.\n")
            }
            else
                P <- P.start
            }
        else {
            stop("Error: P starting values are not conformable.\n")
            }
        return(P)
    }

    ## priro check for mean
    "check.theta" <-
    function(theta.start = NA, ns = ns, y = y, min, max){
        if (is.na(theta.start)[1]){ # draw from uniform with range(y)
            theta <- runif(ns, min=min, max=max)}
        else if (length(theta.start)==ns)
            theta <- theta.start
        else if (length(theta.start)!=ns) {
            stop("Error: theta starting values are not conformable.\n")
            }
        return(theta)
    }
