\name{ddirichlet}
\alias{ddirichlet}
\title{Evaluate Density of Dirichlet Distribution}
\description{
  Functions to compute the density of the Dirichlet distribution.
}
\usage{
  ddirichlet(x, alpha)
}
\arguments{
  \item{x}{A vector containing a single deviate or matrix
    containing one random deviate per row.} 
  \item{alpha}{Vector of shape parameters, or matrix of shape
    parameters corresponding to each row above. }
}
\details{
  The Dirichlet distribution is the multidimensional generalization of
  the beta distribution. 
}
\value{
  \code{ddirichlet} returns a vector containing the evaluated Dirichlet
  density for the corresponding rows of \code{x}.
}

\author{
  Code is taken from Greg's Miscellaneous Functions (gregmisc).  His
  code was based on code posted by Ben Bolker to R-News on Fri Dec 15
  2000.
}

\seealso{
  \code{\link[base]{Beta}},
  \code{\link{rdirichlet}}
}

\examples{
  ddirichlet(c(.1,.2,.7), c(1,1,1))
}
\keyword{distribution}
