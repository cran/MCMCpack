% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distn.R
\name{Wishart}
\alias{Wishart}
\alias{dwish}
\alias{rwish}
\title{The Wishart Distribution}
\usage{
rwish(v, S)

dwish(W, v, S)
}
\arguments{
\item{v}{Degrees of freedom (scalar).}

\item{S}{Inverse scale matrix \eqn{(p \times p)}.}

\item{W}{Positive definite matrix W \eqn{(p \times p)}.}
}
\value{
\code{dwish} evaluates the density at positive definite matrix W.
\code{rwish} generates one random draw from the distribution.
}
\description{
Density function and random generation from the Wishart distribution.
}
\details{
The mean of a Wishart random variable with \code{v} degrees of freedom and
inverse scale matrix \code{S} is \eqn{vS}.
}
\examples{

density <- dwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
draw <- rwish(3, matrix(c(1,.3,.3,1),2,2))

}
\keyword{distribution}
