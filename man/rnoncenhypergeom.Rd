\name{rnoncenhypergeom}
\alias{rnoncenhypergeom}
\title{Generate Random Draw from Noncentral Hypergeometric Distribution}
\description{
  This function generates random draws from the Noncentral  
  Hypergeometric distribution.
}

\usage{
  rnoncenhypergeom(n, n1, n2, m1, psi)
}

\arguments{
  \item{n}{The number of draws to make from the distribution.}
  \item{n1}{The size of group one.}
  \item{n2}{The size of group two.}
  \item{m1}{The observed number of positive outcomes (in both groups).}
  \item{psi}{Odds ratio.}
}

\value{
  An list of \code{n} random draws from the distribution.
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
  rnoncenhypergeom(1, 500, 500, 500, 6.0)
}

\keyword{distribution}

\seealso{
   \code{\link{dnoncenhypergeom}}
}
