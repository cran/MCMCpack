% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scythe.R
\name{read.Scythe}
\alias{read.Scythe}
\title{Read a Matrix from a File written by Scythe}
\usage{
read.Scythe(infile = NA)
}
\arguments{
\item{infile}{The file to be read. This can include path information.}
}
\value{
A matrix containing the data stored in the read file.
}
\description{
This function reads a matrix from an ASCII file in the form produced by the
Scythe Statistical Library.  Scythe output files contain the number of rows
and columns in the first row, followed by the data.
}
\examples{

  \dontrun{
  mymatrix <- read.Scythe("myfile.txt")
  }

}
\references{
Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.
\emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu.s3-website-us-east-1.amazonaws.com/}.
}
\seealso{
\code{\link{write.Scythe}}
}
\keyword{file}
