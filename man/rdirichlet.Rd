\name{rdirichlet}
\alias{rdirichlet}
\title{Generate Random Draws from the Dirichlet Distribution}
\description{
  Generate random deviates from the Dirichlet distribution.
}
\usage{
rdirichlet(n, alpha)
}
\arguments{
  \item{n}{Number of random vectors to generate. }
  \item{alpha}{Vector of shape parameters, or matrix of shape
    parameters corresponding to draw.}
}
\details{
  The Dirichlet distribution is the multidimensional generalization of
  the beta distribution. 
}

\value{
  \code{rdirichlet} returns a matrix with \code{n} rows, each containing
  a single Dirichlet random deviate.
}
\author{
  Code is taken from Greg's Miscellaneous Functions (gregmisc).  His
  code was based on code posted by Ben Bolker to R-News on Fri Dec 15
  2000.
}
  
\seealso{ \code{\link[base]{Beta}}, \code{\link{ddirichlet}} }

\examples{
  x <- rdirichlet(20, c(1,1,1) )
}

\keyword{distribution}
