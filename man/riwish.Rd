\name{riwish}
\alias{riwish}
\title{Generate Random Draw from Inverse Wishart Distribution}
\description{
  This function generates a random draw from the Inverse Wishart
  distribution.
}
  
\usage{
   riwish(v, S)
}

\arguments{
  \item{v}{Inverse Wishart degrees of freedom (scalar).}
  \item{S}{Inverse Wishart scale matrix \eqn{(p \times p)}{(p x p)}.}}

\value{
  A \eqn{(p \times p)}{(p x p)} matrix.
}

\examples{
A <- riwish(3, matrix(c(1,.3,.3,1),2,2))
}

\keyword{distribution}

\seealso{\code{\link{diwish}}}
