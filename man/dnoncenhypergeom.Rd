\name{dnoncenhypergeom}
\alias{dnoncenhypergeom}
\title{Evaluate Density of Noncentral Hypergeometric Distribution}
\description{
  This function evaluates the density of the the Noncentral  
  Hypergeometric distribution at a single point or all points.
}

\usage{
dnoncenhypergeom(x=NA, n1, n2, m1, psi)
}

\arguments{
  \item{x}{The location to evaluate the density.  If \code{x} is NA,
    then a matrix is returned with the density evaluated at all possible
    points.}
  \item{n1}{The size of group one.}
  \item{n2}{The size of group two.}
  \item{m1}{The observed number of positive outcomes (in both groups).}
  \item{psi}{Odds ratio.}
}

\value{
  The density evaluated at point \code{x}, or a matrix with
  the first column containing the possible values of the random
  variable, and the second column containing the probabilities.
}

\details{
  The Noncentral Hypergeometric is particularly useful for conditional
  inference for \eqn{(2 \times 2)}{(2 x 2)} tables.  We use the
  parameterization and algorithms of Liao and Rosen (2001).  The
  underlying R code is based on their published code.  See their
  article for details of the parameterization.
}

\source{
  L. G. Liao and Ori Rosen. 2001. ``Fast and Stable Algorithms for Computing
  and Sampling From the Noncentral Hypergeometric Distribution." \emph{The
  American Statistician.} 55: 366-369.
}

\examples{
  dnoncenhypergeom(NA, 500, 500, 500, 6.0)
}

\keyword{distribution}

\seealso{
   \code{\link{rnoncenhypergeom}}
}
