\name{xpnd}
\alias{xpnd}
\title{Expand a Vector into a Symmetric Matrix}

\description{
  This function takes a vector of appropriate length (typically
  created using \code{vech}) and creates a symmetric matrix.
}
  
\usage{
   xpnd(x, nrow)
}

\arguments{
  \item{x}{A list of elements to expand into symmetric matrix.}
  \item{nrow}{The number of rows (and columns) in the returned matrix.}
}

\value{
  An \eqn{(nrows \times nrows)}{(nrows * nrows)} symmetric matrix.
}

\details{
  This function is particularly useful when dealing with variance
  covariance matrices. Note that R stores matrices in column major
  order, and that the items in \code{x} will be recycled to fill the
  matrix if need be.
}


\examples{
  xpnd(c(1,2,3,4,4,5,6,7,8,9),4)
}

\keyword{manip}

\seealso{\code{\link{vech}}}
