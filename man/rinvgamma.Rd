\name{rinvgamma}
\alias{rinvgamma}
\title{Generate Random Draw from Inverse Gamma Distribution}
\description{
  This function generates a random draw from the Inverse Gamma distribution.
}
  
\usage{
   rinvgamma(n, shape, rate = 1)
}

\arguments{
  \item{n}{Number of draws from the distribution.}
  \item{shape}{Scalar shape parameter.}
  \item{rate}{Scalar rate parameter (default value one).}
}

\value{
  \code{n} draws from the Inverse Gamma distribution.  The
  parameterization is consistent with the Gamma Distribution in the
  base package.
}

\examples{
  draws <- rinvgamma(10, 3.2)
}

\keyword{distribution}

\seealso{
  \code{\link{dinvgamma}},
  \code{\link[base]{GammaDist}}
}


