\name{dwish}
\alias{dwish}
\title{Evaluate the PDF of the Wishart Distribution}
\description{
  This evaluates the PDF of the Wishart Distribution at the 
  positive definite matrix W.
  }
  
\usage{
   dwish(W, v, S)
   }

\arguments{
    \item{W}{Positive definite matrix W \eqn{(p \times p)}{(p x p)}.}
    \item{v}{inverse Wishart degrees of freedom (scalar).}
    \item{S}{inverse Wishart scale matrix \eqn{(p \times p)}{(p x p)}.}}

\value{
  The density evaluated at W.
}

\examples{
x <- dwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
}

\keyword{distribution}


