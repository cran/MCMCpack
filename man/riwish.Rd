\name{riwish}
\alias{riwish}
\title{Generate Random Draw from inverse Wishart Distribution}
\description{
  This function generates a random draw from the inverse Wishart distribution.
  }
  
\usage{
   riwish(v, S)
   }

\arguments{
    \item{v}{inverse Wishart degrees of freedom (scalar).}
    \item{S}{inverse Wishart scale matrix \eqn{(p \times p)}{(p x p)}.}}

\value{
   A \eqn{(p \times p)}{(p x p)} matrix.
}

\examples{
A <- riwish(3, matrix(c(1,.3,.3,1),2,2))
}

\keyword{distribution}
