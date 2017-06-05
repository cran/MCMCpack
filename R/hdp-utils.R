#' Posterior Changepoint Probabilities from HDP-HMM
#'
#' Plot the posterior density of regime change.
#'
#' @param mcmcout The \code{mcmc} object containing the posterior density
#' sample from a changepoint model. Note that this must be from a
#'   HDP-HMM sampler.
#'
#' @param main Title of the plot
#'
#' @param xlab Label for the x-axis.
#'
#' @param ylab Label for the y-axis.
#'
#' @param start The time of the first observation to be shown in the time
#' series plot.
#'
#' @export
#'
#' @seealso \code{\link{HDPHMMpoisson}}, \code{\link{HDPHMMnegbin}}, \code{\link{HDPHSMMnegbin}}
#'
#' @keywords hplot
"plotHDPChangepoint" <-
  function (mcmcout, main="Posterior Changepoint Probabilities", xlab = "Time", ylab = "",
            start = 1){
      pr.st <- attr(mcmcout, "pr.chpt")
      y <- attr(mcmcout, "y")
      K <- attr(mcmcout, "K")
      if (!is.ts(y))
          y <- ts(y, start)
      time.frame <- as.vector(time(y))
      stopifnot(K > 1)
      plot(time.frame, pr.st, type = "h", lwd=2, main = main, ylim=c(0, max(pr.st)), xlab = xlab, ylab = ylab, axes=FALSE)
      axis(1, tick = FALSE, col="darkgrey", lwd=0.5)
      axis(2, tick = FALSE, col="darkgrey", lwd=0.5)
      for (i in 1:length(axTicks(1))) lines(c(axTicks(1)[i], axTicks(1)[i]), c(0,1), col="darkgrey", lwd=0.5)
      for (i in 1:length(axTicks(2))) lines(c(axTicks(2)[1], max(time.frame)),
                                            c(axTicks(2)[i], axTicks(2)[i]), col="darkgrey", lwd=0.5)
      invisible()
  }

## prior check for transition matrix
"check.hdp.P" <- function(P.start = NA, K = K, n = n, theta = theta){
  if (is.na(P.start)[1]){
    P <- matrix((1 - theta) / (K - 1), nrow = K, ncol = K)
    diag(P) <- theta
  } else if ((dim(P.start)[1] == K)&(dim(P.start)[2] == K)){
    if ((max(P.start) > 1) || (min(P.start) < 0)) {
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
