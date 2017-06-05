#' @export
"print.qrssvs"<-function(x, ...){
  x.orig<-x
  cat("Quantile regression stochastic search \nvariable selection (QR-SSVS) output:\nStart = ",
      attr(x,"mcpar")[1], "\nEnd = ",
      attr(x,"mcpar")[2], "\nThinning interval = ",
      attr(x,"mcpar")[3], "\n")

  attr(x, "mcpar") <- NULL
  attr(x, "class") <- NULL
  NextMethod("print", ...)
  invisible(x.orig)
}

#' Calculate the marginal posterior probabilities of predictors being included
#' in a quantile regression model.
#'
#' This function extracts the marginal probability table produced by
#' \code{summary.qrssvs}.
#'
#' @param qrssvs An object of class \code{qrssvs}. Typically this will be the
#' \code{gamma} component of the list returned by \code{SSVSquantreg}.
#'
#' @return A table with the predictors listed together with their posterior
#' marginal posterior probability of inclusion.
#'
#' @export
#'
#' @author Craig Reed
#'
#' @seealso \code{\link[MCMCpack]{SSVSquantreg}}
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' set.seed(1)
#' epsilon<-rnorm(100)
#' set.seed(2)
#' x<-matrix(rnorm(1000),100,10)
#' y<-x[,1]+x[,10]+epsilon
#' qrssvs<-SSVSquantreg(y~x)
#' mptable(qrssvs$gamma)
#' }
#'
"mptable" <- function(qrssvs){
  if (!is(qrssvs, "qrssvs")){
    stop("Can only be used on objects of class qrssvs.\n")
  }
  ssvs.start <- attr(qrssvs, "mcpar")[1]
  ssvs.end <- attr(qrssvs, "mcpar")[2]
  ssvs.thin <- attr(qrssvs, "mcpar")[3]
  nstore <- (ssvs.end-ssvs.start)/ssvs.thin + 1
  probs<-apply(qrssvs,2,function(z){length(which(z==1))})/nstore
  return(data.frame(Probability=probs))
}

#' Shows an ordered list of the most frequently visited models sampled during
#' quantile regression stochastic search variable selection (QR-SSVS).
#'
#' Given output from quantile regression stochastic search variable selection,
#' this function returns a table of the 'best' models together with their
#' associated empirical posterior probability.
#'
#' @param qrssvs An object of class \code{qrssvs}. Typically this will be the
#' \code{gamma} component of the list returned by \code{SSVSquantreg}.
#'
#' @param nmodels The number of models to tabulate.
#'
#' @param abbreviate Logical: should the names of the predictors be
#' abbreviated?
#'
#' @param minlength If \code{abbreviate} is set to \code{TRUE}, the minimum
#' length of the abbreviations.
#'
#' @return A table with the models and their associated posterior probability.
#' The models are arranged in descending order of probability.
#'
#' @export
#'
#' @author Craig Reed
#'
#' @seealso \code{\link[MCMCpack]{SSVSquantreg}}
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' set.seed(1)
#' epsilon<-rnorm(100)
#' set.seed(2)
#' x<-matrix(rnorm(1000),100,10)
#' y<-x[,1]+x[,10]+epsilon
#' qrssvs<-SSVSquantreg(y~x)
#' topmodels(qrssvs$gamma)
#' }
#'
"topmodels"<-function(qrssvs, nmodels=5, abbreviate=FALSE, minlength=3){
  if (!is(qrssvs, "qrssvs")){
    stop("Can only be used on objects of class qrssvs.\n")
  }
  ssvs.start <- attr(qrssvs, "mcpar")[1]
  ssvs.end <- attr(qrssvs, "mcpar")[2]
  ssvs.thin <- attr(qrssvs, "mcpar")[3]
  nstore <- (ssvs.end-ssvs.start)/ssvs.thin + 1
  xnames <- attr(qrssvs, "xnames")
  if (abbreviate){
    xnames <- abbreviate(xnames, minlength)
  }
  model.list<-apply(qrssvs,1,function(z)xnames[which(z==1)])
  model.vector<-sapply(model.list, function(z)paste(z, collapse=","))
  model.count<-sort(table(model.vector), decreasing=T)/nstore
  if (nmodels>length(model.count)){
    warning("Number of models requested exceeds total number of models visited.\n")
  }
  if (rownames(model.count)[1]==""){
    rownames(model.count)[1]<-"Null model"
  }
  return(data.frame(Probability=model.count[1:(min(nmodels, length(model.count)))]))
}

#' Plot output from quantile regression stochastic search variable selection
#' (QR-SSVS).
#'
#' This function produces a Trellis plot of the predictors on the y-axis versus
#' the marginal posterior probability of inclusion on the x-axis.
#'
#' @param x An object of class \code{qrssvs}. Typically this will be the
#' \code{gamma} component of the list returned by \code{SSVSquantreg}.
#'
#' @param ... Further arguments
#'
#' @return An object with class \code{"trellis"}. The associated
#' \code{\link[lattice:update.trellis]{update}} and
#' \code{\link[lattice:print.trellis]{print}} methods are documented in the
#' "Lattice" package.
#'
#' @export
#'
#' @author Craig Reed
#'
#' @seealso \code{\link[MCMCpack]{SSVSquantreg}},
#' \code{\link[MCMCpack]{mptable}}, \code{\link[lattice:Lattice]{Lattice}} for
#' a brief introduction to lattice displays and links to further documentation.
#'
#' @references Deepayan Sarkar. 2008. \emph{lattice: Lattice Graphics.} R package version
#' 0.17-17
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' set.seed(1)
#' epsilon<-rnorm(100)
#' set.seed(2)
#' x<-matrix(rnorm(1000),100,10)
#' y<-x[,1]+x[,10]+epsilon
#' qrssvs<-SSVSquantreg(y~x)
#' plot(qrssvs$gamma)
#' ## Modify the graph by increasing the fontsize on the axes
#' qrssvsplot<-plot(qrssvs$gamma)
#' update(qrssvsplot, scales=list(cex=3))
#' }
#'
"plot.qrssvs"<-function(x, ...){
  probs<-mptable(x)
  dotplot(as.matrix(probs),
          panel=function(x, y, ...){
            panel.abline(v=0.5, lty=3)
            panel.dotplot(x, y, ...)
          },
          origin=0, type=c("p","h"), pch=16,
          xlim=c(-0.05,1.05),
          scales = list(x = list(at = c(0,0.2,0.4,0.5,0.6,0.8,1))),
          xlab="Marginal inclusion probability", ...)
}

#' Summarising the results of quantile regression stochastic search variable
#' selection (QR-SSVS).
#'
#' This function produces a table of predictors and their associated marginal
#' posterior probability of inclusion. It also returns the median probability
#' model (see the details section).
#'
#' The median probability model is defined to be the model that contains any
#' predictor with marginal posterior probability greater than or equal to 0.5.
#' If the goal is to select a single model e.g. for prediction, Barbieri and
#' Berger (2004) recommend the median probability model. In some cases, this
#' will coincide with the maximum probability model.
#'
#' @aliases summary.qrssvs print.summary.qrssvs
#' @name summaryqrssvs
#'
#' @param object An object of class \code{qrssvs}. Typically this will be the
#' \code{gamma} component of the list returned by \code{SSVSquantreg}.
#'
#' @param ... Further arguments.
#'
#' @export
#'
#' @author Craig Reed
#'
#' @seealso \code{\link[MCMCpack]{SSVSquantreg}},
#' \code{\link[MCMCpack]{mptable}}, \code{\link[MCMCpack]{topmodels}}
#'
#' @references Maria M. Barbieri, and James O. Berger (2004). "Optimal predictive model
#' selection". \emph{Annals of Statistics}, 32, 870-897.
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#' set.seed(1)
#' epsilon<-rnorm(100)
#' set.seed(2)
#' x<-matrix(rnorm(1000),100,10)
#' y<-x[,1]+x[,10]+epsilon
#' qrssvs<-SSVSquantreg(y~x)
#' summary(qrssvs$gamma)
#' }
#'
"summary.qrssvs"<-function(object, ...){
  covnames <- attr(object, "xnames")
  probs<-mptable(object)
  median.model<-covnames[probs>=0.5]
  results<-probs
  attr(results, "median.model")<-median.model
  attr(results, "tau")<-attr(object, "tau")
  class(results)<-"summary.qrssvs"
  return(results)
}

#' @export
"print.summary.qrssvs"<-function(x, digits=max(3, .Options$digits-3), ...){
  attr(x, "class")<-"data.frame"
  cat("\nMarginal inclusion probability of each predictor:\n\n")
  print(x, digits=digits, ...)
  cat("\nFor tau = ", attr(x,"tau"), ", the median probability model \nincludes the following predictors:\n\n",
      paste(attr(x, "median.model"), collapse=", "), ".\n\n", sep="")
  invisible(x)
}

