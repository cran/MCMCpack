\name{tomogplot}
\alias{tomogplot}
\title{Tomography Plot}
\description{
  tomogplot is used to produce a tomography plot (see King, 1997) for a
  series of partially observed 2 x 2 contingency tables.
}
  
\usage{
tomogplot(r0, r1, c0, c1, xlab="fraction of r0 in c0 (p0)",
          ylab="fraction of r1 in c0 (p1)", bgcol="white", ...)
   }

\arguments{
  \item{r0}{(ntables * 1) vector of row sums from row 0}

  \item{r1}{(ntables * 1) vector of row sums from row 1}

  \item{c0}{(ntables * 1) vector of column sums from column 0}

  \item{c1}{(ntables * 1) vector of column sums from column 1}

  \item{xlab}{the x axis label for the plot}

  \item{ylab}{the y axis label for the plot}

  \item{bgcol}{the background color for the plot}
  
  \item{...}{further arguments to be passed}     
}


\details{
  Consider the following partially observed 2 by 2 contingency table:\cr
  \cr
  \tabular{llll}{
               \tab | \eqn{Y=0} \tab | \eqn{Y=1} \tab |   \cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=0} \tab | \eqn{Y_0}{Y0} \tab |  \tab | \eqn{r_0}{r0}\cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=1} \tab | \eqn{Y_1}{Y1} \tab |  \tab | \eqn{r_1}{r1}\cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
              \tab | \eqn{c_0}{c0} \tab | \eqn{c_1}{c1} \tab | \eqn{N}\cr    
  }

  where \eqn{r_0}{r0}, \eqn{r_1}{r1}, \eqn{c_0}{c0}, \eqn{c_1}{c1}, and
  \eqn{N}  are non-negative integers that are
  observed. The interior cell entries are not observed. It is
  assumed that \eqn{Y_0|r_0 \sim \mathcal{B}inomial(r_0,
    p_0)}{Y0|r0 ~ Binomial(r0, p0)} and
  \eqn{Y_1|r_1 \sim \mathcal{B}inomial(r_1, p_1)}{Y1|r1 ~
    Binomial(r1,p1)}.

  tomogplot plots the bounds on the MLEs for (p0, p1). 
  }

  \keyword{hplot}
  
  \references{
    Gary King, 1997. A Solution to the Ecological Inference Problem.
    Princeton: Princeton University Press.
    
    Jonathan Wakefield. 2001. ``Ecological Inference for 2 x 2 Tables,''
    Center for Statistics and the Social Sciences Working Paper 
    no. 12. University of Washington.  
}

\examples{
r0 <- rpois(100, 500) 
r1 <- rpois(100, 200) 
c0 <- rpois(100, 100) 
c1 <- (r0 + r1) - c0 
tomogplot(r0, r1, c0, c1) 
}

\seealso{\code{\link{MCMCbaselineEI}}, \code{\link{MCMChierEI}},
  \code{\link{MCMCdynamicEI}}, 
  \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}

