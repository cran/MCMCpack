\name{read.Scythe}
\alias{read.Scythe}
\title{Read a Matrix from a File written by Scythe}
\description{
  This function reads a matrix from an ASCII file in the form produced
  by the Scythe Statistical Library.  Scythe
  output files contain the number of rows and columns in the first row,
  followed by the data.
}

\usage{
   read.Scythe(infile=NA)
}

\arguments{
    \item{infile}{The file to be read. This can include path information.}
}

\value{
  A matrix containing the data stored in the read file.
}

\examples{
  \dontrun{mymatrix <- read.Scythe("myfile.txt")}
}

\references{
  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.
  \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
}

\keyword{file}

\seealso{\code{\link{write.Scythe}}}
