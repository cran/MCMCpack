\name{plot.qrssvs}
\alias{plot.qrssvs}
\title{Plot output from quantile regression stochastic search variable selection (QR-SSVS).}

\description{This function produces a Trellis plot of the predictors on the y-axis versus the marginal posterior probability of inclusion on the x-axis.}

\usage{
\method{plot}{qrssvs}(x, \dots)
}

\arguments{
  \item{x}{An object of class \code{qrssvs}. Typically this will be the \code{gamma} component of the list returned by \code{SSVSquantreg}.}

  \item{\dots}{Further arguments}
}

\value{An object with class \code{"trellis"}. The associated \code{\link[lattice:update.trellis]{update}} and \code{\link[lattice:print.trellis]{print}} methods are documented in the "Lattice" package.}

\author{Craig Reed}

\references{
   
   Deepayan Sarkar. 2008. \emph{lattice: Lattice Graphics.} R package version 0.17-17  
}

\examples{
\dontrun{
set.seed(1)
epsilon<-rnorm(100)
set.seed(2)
x<-matrix(rnorm(1000),100,10)
y<-x[,1]+x[,10]+epsilon
qrssvs<-SSVSquantreg(y~x)
plot(qrssvs$gamma)
## Modify the graph by increasing the fontsize on the axes
qrssvsplot<-plot(qrssvs$gamma)
update(qrssvsplot, scales=list(cex=3))
}
}

\keyword{models}

\seealso{
\code{\link[MCMCpack]{SSVSquantreg}},
\code{\link[MCMCpack]{mptable}}, 
\code{\link[lattice:Lattice]{Lattice}} for a brief introduction to
  lattice displays and links to further documentation.
}


