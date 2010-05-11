\name{mptable}
\alias{mptable}

\title{Calculate the marginal posterior probabilities of predictors being included in a quantile regression model.}

\description{This function extracts the marginal probability table produced by \code{summary.qrssvs}.}

\usage{mptable(qrssvs)}

\arguments{
  \item{qrssvs}{An object of class \code{qrssvs}. Typically this will be the \code{gamma} component of the list returned by \code{SSVSquantreg}.}
}

\value{A table with the predictors listed together with their posterior marginal posterior probability of inclusion.}

\author{Craig Reed}

\examples{
\dontrun{
set.seed(1)
epsilon<-rnorm(100)
set.seed(2)
x<-matrix(rnorm(1000),100,10)
y<-x[,1]+x[,10]+epsilon
qrssvs<-SSVSquantreg(y~x)
mptable(qrssvs$gamma)
}
}

\keyword{models}

\seealso{
\code{\link[MCMCpack]{SSVSquantreg}}}
