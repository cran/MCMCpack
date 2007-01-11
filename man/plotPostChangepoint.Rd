\name{plotPostChangepoint}
\alias{plotPostChangepoint}
\title{Changepoint Location Plots}
\description{Plot the posterior distribution over the location of the changepoints.}

\usage{
   plotPostChangepoint(mcmcout, y, m, xlab="Time", ylab="")
}

\arguments{

\item{mcmcout}{The \code{mcmc} object containing the posterior density sample from a changepoint model.  Note that this must have a \code{prob.state} attribute.}

\item{y}{The data.}

\item{m}{The number of changepoints.}

\item{xlab}{Label for the x-axis.}

\item{ylab}{Label for the y-axis.}

}

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\keyword{hplot}

\seealso{\code{\link{MCMCpoissonChangepoint}}}
