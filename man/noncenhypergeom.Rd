\name{NoncenHypergeom}
\alias{NoncenHypergeom}
\alias{rnoncenhypergeom}
\alias{dnoncenhypergeom}
\title{The Noncentral Hypergeometric Distribution}
\description{
  Evaluates the density at a single point or all points, and generate random 
  draws from the Noncentral Hypergeometric distribution.
}

\usage{
dnoncenhypergeom(x=NA, n1, n2, m1, psi)
rnoncenhypergeom(n, n1, n2, m1, psi)
}

\arguments{
  \item{x}{The location to evaluate the density.  If \code{x} is NA,
    then a matrix is returned with the density evaluated at all possible
    points.}
  \item{n}{The number of draws to make from the distribution.}
  \item{n1}{The size of group one.}
  \item{n2}{The size of group two.}
  \item{m1}{The observed number of positive outcomes (in both groups).}
  \item{psi}{Odds ratio.}
}

\value{
  \code{dnoncenhypergeom} evaluates the density at point \code{x},
  or a matrix with
  the first column containing the possible values of the random
  variable, and the second column containing the probabilities.
  \code{rnoncenhypergeom} returns a list of \code{n} random draws from
  the distribution.
}

\details{
  The Noncentral Hypergeometric is particularly useful for conditional
  inference for \eqn{(2 \times 2)}{(2 x 2)} tables.  We use the
  parameterization and algorithms of Liao and Rosen (2001).  The
  underlying R code is based on their published code.  See their
  article for details of the parameterization.
}

\source{
  J. G. Liao and Ori Rosen. 2001. ``Fast and Stable Algorithms for Computing
  and Sampling From the Noncentral Hypergeometric Distribution." \emph{The
  American Statistician.} 55: 366-369.
}

\examples{
  density <- dnoncenhypergeom(NA, 500, 500, 500, 6.0)
  draws <- rnoncenhypergeom(10, 500, 500, 500, 6.0)
}

\keyword{distribution}


