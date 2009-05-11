\name{summary.qrssvs}
\alias{summary.qrssvs}
\alias{print.summary.qrssvs}
\title{Summarising the results of quantile regression stochastic search variable selection (QR-SSVS).}

\description{This function produces a table of predictors and their associated marginal posterior probability of inclusion. It also returns the median probability model (see the details section).}

\usage{\method{summary}{qrssvs}(object, \dots)}

\arguments{
  \item{object}{An object of class \code{qrssvs}. Typically this will be the \code{gamma} component of the list returned by \code{SSVSquantreg}.}

  \item{\dots}{Further arguments.}
}

\details{The median probability model is defined to be the model that contains any predictor with marginal posterior probability greater than or equal to 0.5. If the goal is to select a single model e.g. for prediction, Barbieri and Berger (2004) recommend the median probability model. In some cases, this will coincide with the maximum probability model.}

\author{Craig Reed}

\references{
   
  Maria M. Barbieri, and James O. Berger (2004). "Optimal predictive model selection". \emph{Annals of Statistics}, 32, 870-897.
}

\examples{
\dontrun{
set.seed(1)
epsilon<-rnorm(100)
set.seed(2)
x<-matrix(rnorm(1000),100,10)
y<-x[,1]+x[,10]+epsilon
qrssvs<-SSVSquantreg(y~x)
summary(qrssvs$gamma)
}
}

\keyword{models}

\seealso{
\code{\link[MCMCpack]{SSVSquantreg}},
\code{\link[MCMCpack]{mptable}},
\code{\link[MCMCpack]{topmodels}}}
