\name{InvGamma}
\alias{dinvgamma}
\alias{rinvgamma}
\alias{InvGamma}
\title{The Inverse Gamma Distribution}
\description{
  Density function and random generation from the inverse gamma distribution.
}

\usage{
rinvgamma(n, shape, scale = 1)
dinvgamma(x, shape, scale = 1)
}

\arguments{
  \item{x}{Scalar location to evaluate density.}
  \item{n}{Number of draws from the distribution.}
  \item{shape}{Scalar shape parameter.}
  \item{scale}{Scalar scale parameter (default value one).}
}

\value{
  \code{dinvgamma} evaluates the density at \code{x}.  \code{rinvgamma} takes
  \code{n} draws from the inverse Gamma distribution.  The parameterization is
  consistent with the Gamma Distribution in the stats package.
}

\details{
An inverse gamma random variable with shape \eqn{a}{a} and scale
\eqn{b}{b} has mean \eqn{\frac{b}{a-1}}{b/(a-1)} (assuming \eqn{a>1}{a>1}) and
variance \eqn{\frac{b^2}{(a-1)^2(a-2)}}{(b^2)/((a-1)^2 (a-2))} (assuming
\eqn{a>2}{a>2}).
}

\examples{
density <- dinvgamma(4.2, 1.1)
draws <- rinvgamma(10, 3.2)
}

\keyword{distribution}

\seealso{\code{\link[stats]{GammaDist}}}
