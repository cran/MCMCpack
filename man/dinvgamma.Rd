\name{dinvgamma}
\alias{dinvgamma}
\title{Evaluate the Density of the Inverse Gamma Distribution}
\description{
  This evaluates the density of the Inverse Gamma Distribution at the
  value x.
}

\usage{
   dinvgamma(x, shape, rate=1)
}

\arguments{
  \item{x}{Scalar location to evaluate density.}
  \item{shape}{Scalar shape parameter.}
  \item{rate}{Scalar rate parameter (default value one).}
}

\value{
  The density evaluated at \code{x}.  The parameterization is
  consistent with the Gamma Distribution in the base package.
}

\examples{
density <- dinvgamma(4.2, 1.1)
}

\keyword{distribution}

\seealso{
  \code{\link{rinvgamma}},
  \code{\link[base]{GammaDist}}
  }
