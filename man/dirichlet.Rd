\name{Dirichlet}
\alias{Dirichlet}
\alias{ddirichlet}
\alias{rdirichlet}
\title{The Dirichlet Distribution}
\description{
  Density function and random generation from the Dirichlet distribution.
}
\usage{
ddirichlet(x, alpha)
rdirichlet(n, alpha)
}
\arguments{
  \item{x}{A vector containing a single deviate or matrix containing
    one random deviate per row.}
  \item{n}{Number of random vectors to generate. }
  \item{alpha}{Vector of shape parameters, or matrix of shape
    parameters corresponding to the number of draw.}
}
\details{
  The Dirichlet distribution is the multidimensional generalization of
  the beta distribution. 
}

\value{
  \code{ddirichlet} gives the density. \code{rdirichlet} returns a
  matrix with \code{n} rows, each containing a single Dirichlet random
  deviate.
}
\author{
  Code is taken from Greg's Miscellaneous Functions (gregmisc).  His
  code was based on code posted by Ben Bolker to R-News on 15 Dec
  2000.
}
  
\seealso{ \code{\link[stats]{Beta}} }

\examples{
  density <- ddirichlet(c(.1,.2,.7), c(1,1,1))
  draws <- rdirichlet(20, c(1,1,1) )
}

\keyword{distribution}
