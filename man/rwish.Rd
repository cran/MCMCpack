\name{rwish}
\alias{rwish}
\title{Generate Random Draw from Wishart Distribution}
\description{
  This function generates a random draw from the Wishart distribution.
  }
  
\usage{
   rwish(v, S)
   }

\arguments{
    \item{v}{Wishart degrees of freedom (scalar).}
    \item{S}{Wishart scale matrix \eqn{(p \times p)}{(p x p)}.}}

\value{
   A \eqn{(p \times p)}{(p x p)} matrix.
}

\examples{
A <- rwish(3, matrix(c(1,.3,.3,1),2,2))
}

\keyword{distribution}

\seealso{\code{\link{dwish}}}
