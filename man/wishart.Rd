\name{Wishart}
\alias{dwish}
\alias{rwish}
\alias{Wishart}
\title{The Wishart Distribution}
\description{
  Density function and random generation from the Wishart distribution.
}
  
\usage{
   dwish(W, v, S)
   rwish(v, S)
}

\arguments{
    \item{W}{Positive definite matrix W \eqn{(p \times p)}{(p x p)}.}
    \item{v}{Wishart degrees of freedom (scalar).}
    \item{S}{Wishart scale matrix \eqn{(p \times p)}{(p x p)}.}
}

\value{
  \code{dwish} evaluates the density at positive definite matrix W.
  \code{rwish} generates one random draw from the distribution.
}


\examples{
density <- dwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
draw <- rwish(3, matrix(c(1,.3,.3,1),2,2))
}

\keyword{distribution}



