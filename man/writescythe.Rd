\name{write.Scythe}
\alias{write.Scythe}
\title{Write a Matrix to a File to be Read by Scythe}
\description{
  This function writes a matrix to an ASCII file that can be 
  read by the Sycthe Statistical Library.  Scythe requires that
  input files contain the number of rows and columns in the first
  row, followed by the data.
}
  
\usage{
  write.Scythe(outmatrix, outfile=NA, overwrite=FALSE)
}

\arguments{
  \item{outmatrix}{The matrix to be written to a file.}
  \item{outfile}{The file to be written. This can include path information.}
  \item{overwrite}{A logical that determines whether an existing file
    should be over-written.  By default, it protects the user from
    over-writing existing files.}
}

\value{
  A zero if the file is properly written.
}

\examples{
  \dontrun{write.Scythe(mymatrix, "myfile.txt")}
}

\references{
  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  
  \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
}

\keyword{file}

\seealso{\code{\link{write.Scythe}}}

  
