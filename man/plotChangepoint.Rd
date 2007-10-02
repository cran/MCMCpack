\name{plotChangepoint}
\alias{plotChangepoint}
\title{Changepoint Location Plots}
\description{Plot the posterior distribution over the location of the changepoints.}

\usage{
   plotChangepoint(mcmcout, xlab="Time", ylab="")
}

\arguments{

\item{mcmcout}{The \code{mcmc} object containing the posterior density sample from a changepoint model.  Note that this must have a \code{prob.state} attribute.}

\item{xlab}{Label for the x-axis.}

\item{ylab}{Label for the y-axis.}

}

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\keyword{hplot}

\seealso{\code{\link{MCMCpoissonChangepoint}}}
