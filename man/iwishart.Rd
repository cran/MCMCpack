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
    \item{v}{Inverse Wishart degrees of freedom (scalar).}
    \item{S}{Inverse Wishart scale matrix \eqn{(p \times p)}{(p x p)}.}}

\value{
  \code{diwish} evaluates the density at positive definite matrix W.
  \code{riwish} generates one random draw from the distribution.
}

\examples{
density <- diwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
draw <- riwish(3, matrix(c(1,.3,.3,1),2,2))
}

\keyword{distribution}



