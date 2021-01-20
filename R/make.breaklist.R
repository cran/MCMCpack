#' Vector of break numbers
#'
#' This function generates a vector of break numbers using the output of
#' \code{testpanelSubjectBreak}.  The function performs a pairwise comparison
#' of models using Bayes Factors.
#'
#' @param BF output of \code{testpanelSubjectBreak}.
#'
#' @param threshold The Bayes Factor threshold to pick the best model.  If a
#' Bayes factor of two models is smaller than \code{threshold}, the model with
#' a smaller number of break is chosen to avoid the over-identification
#' problem.  Users can change threshold into any positive number.  The default
#' value of 3 is chosen as it indicates the existence of "substantial evidence"
#' in favor of the model in the numerator according to Jeffreys' scale.
#'
#' @return Vector fo break numbers.
#'
#' @export
#'
#' @seealso \code{\link{testpanelSubjectBreak}}
#'
#' @references Jong Hee Park, 2012. ``Unified Method for Dynamic and
#'   Cross-Sectional Heterogeneity: Introducing Hidden Markov Panel
#'   Models.''  \emph{American Journal of Political Science}.56:
#'   1040-1054. <doi: 10.1111/j.1540-5907.2012.00590.x>
#'
#' Harold Jeffreys, 1961. The Theory of Probability. Oxford University Press.
"make.breaklist" <- function(BF, threshold=3){
  eBF <- exp(BF)
  N <- nrow(BF)
  out <- rep(NA, N)
  for (i in 1:N){
    order.i <- order(eBF[i,], decreasing=TRUE)
    if(sum(is.na(eBF[i,]))>0){
      out[i] <- 1
    }
    else{
      if(eBF[i, order.i[1]] / eBF[i, order.i[2]] > threshold){
        out[i] <- order.i[1]
      }
      else if (eBF[i, order.i[1]] / eBF[i, order.i[2]] < threshold & order.i[1]<=order.i[2]){
        out[i] <- order.i[1]
      }
      else if (eBF[i, order.i[1]] / eBF[i, order.i[2]] < threshold & order.i[1]>order.i[2]){
        out[i] <- order.i[2]
      }
      else{
        cat("\n Error occurs at i = ", i)
      }
    }
  }
  return(out-1)
}
