\name{InvWishart}
\alias{diwish}
\alias{riwish}
\alias{InvWishart}
\title{The Inverse Wishart Distribution}
\description{
  Density function and random generation from the Inverse Wishart distribution.
}
  
\usage{
   diwish(W, v, S)
   riwish(v, S)
}

\arguments{
    \item{W}{Positive definite matrix W \eqn{(p \times p)}{(p x p)}.}
    \item{v}{Degrees of freedom (scalar).}
    \item{S}{Scale matrix \eqn{(p \times p)}{(p x p)}.}}

\value{
  \code{diwish} evaluates the density at positive definite matrix W.
  \code{riwish} generates one random draw from the distribution.
}

\details{
  The mean of an inverse Wishart random variable with \code{v} degrees
  of freedom and scale matrix \code{S} is \eqn{(v-p-1)^{-1}S}{1/(v-p-1)
    S}. 
}

\examples{
density <- diwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
draw <- riwish(3, matrix(c(1,.3,.3,1),2,2))
}

\keyword{distribution}



