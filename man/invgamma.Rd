\name{InvGamma}
\alias{dinvgamma}
\alias{rinvgamma}
\alias{InvGamma}
\title{The Inverse Gamma Distribution}
\description{
  Density function and random generation from the inverse Gamma distribution.
}

\usage{
rinvgamma(n, shape, rate = 1)
dinvgamma(x, shape, rate = 1)
}

\arguments{
  \item{x}{Scalar location to evaluate density.}
  \item{n}{Number of draws from the distribution.}
  \item{shape}{Scalar shape parameter.}
  \item{rate}{Scalar rate parameter (default value one).}
}

\value{
  \code{dinvgamma} evaluates the density at \code{x}.  \code{rinvgamma} takes
  \code{n} draws from the inverse Gamma distribution.  The parameterization is
  consistent with the Gamma Distribution in the stats package.
}

\examples{
density <- dinvgamma(4.2, 1.1)
draws <- rinvgamma(10, 3.2)
}

\keyword{distribution}

\seealso{\code{\link[stats]{GammaDist}}}
