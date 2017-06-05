##########################################################################
## Tomography Plots for Ecological Inference
##
## produces tomography plots (see King, 1997, A Solution to the
## Ecological Inference Problem, Princeton University Press)
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 11/9/2002
##
## Modification added suggested by David Hugh-Jones 6/10/2006
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Tomography Plot
#'
#' tomogplot is used to produce a tomography plot (see King, 1997) for a series
#' of partially observed 2 x 2 contingency tables.
#'
#' Consider the following partially observed 2 by 2 contingency table:
#'
#' \tabular{llll}{
#'            \tab | \eqn{Y=0} \tab | \eqn{Y=1} \tab |           \cr
#'  --------- \tab ---------   \tab ---------   \tab ---------   \cr
#'  \eqn{X=0} \tab | \eqn{Y_0} \tab |           \tab | \eqn{r_0} \cr
#'  --------- \tab ---------   \tab ---------   \tab ---------   \cr
#'  \eqn{X=1} \tab | \eqn{Y_1} \tab |           \tab | \eqn{r_1} \cr
#'  --------- \tab ---------   \tab ---------   \tab ---------   \cr
#'            \tab | \eqn{c_0} \tab | \eqn{c_1} \tab | \eqn{N}
#' }
#'
#' where \eqn{r_0}, \eqn{r_1}, \eqn{c_0}, \eqn{c_1}, and \eqn{N} are
#' non-negative integers that are observed. The interior cell entries
#' are not observed. It is assumed that \eqn{Y_0|r_0 \sim
#' \mathcal{B}inomial(r_0, p_0)} and \eqn{Y_1|r_1 \sim
#' \mathcal{B}inomial(r_1, p_1)}.
#'
#' This function plots the bounds on the maximum likelihood estimatess for (p0,
#' p1).
#'
#' @param r0 An \eqn{(ntables \times 1)} vector of row sums from
#' row 0.
#'
#' @param r1 An \eqn{(ntables \times 1)} vector of row sums from
#' row 1.
#'
#' @param c0 An \eqn{(ntables \times 1)} vector of column sums
#' from column 0.
#'
#' @param c1 An \eqn{(ntables \times 1)} vector of column sums
#' from column 1.
#'
#' @param xlab The x axis label for the plot.
#'
#' @param ylab The y axis label for the plot.
#'
#' @param bgcol The background color for the plot.
#'
#' @param ... further arguments to be passed
#'
#' @export
#'
#' @seealso \code{\link{MCMChierEI}}, \code{\link{MCMCdynamicEI}},
#' \code{\link{dtomogplot}}
#'
#' @references Gary King, 1997. \emph{A Solution to the Ecological Inference
#' Problem}.  Princeton: Princeton University Press.
#'
#' Jonathan C. Wakefield. 2004. ``Ecological Inference for 2 x 2 Tables.''
#' \emph{Journal of the Royal Statistical Society, Series A}. 167(3): 385445.
#'
#' @keywords hplot
#'
#' @examples
#'
#' r0 <- rpois(100, 500)
#' r1 <- rpois(100, 200)
#' c0 <- rpois(100, 100)
#' c1 <- (r0 + r1) - c0
#' tomogplot(r0, r1, c0, c1)
#'
"tomogplot" <-
  function(r0, r1, c0, c1,
           xlab="fraction of r0 in c0 (p0)",
           ylab="fraction of r1 in c0 (p1)",
           bgcol="white", ...) {
    if (length(r0) != length(r1)) {
      stop("r0 and r1 different lengths in tomogplot().\n")
    }
    if (length(r0) != length(c0)) {
      stop("r0 and c0 different lengths in tomogplot().\n")
    }
    if (length(r0) != length(c1)) {
      stop("r0 and c1 different lengths in tomogplot().\n")
    }

    intercept <-  c0/r1
    slope <- -1 * r0/r1
    N <- length(r0)

    par(pty="s")
    plot(0:1, 0:1, type="n", main="", xlab=xlab, ylab=ylab)
    rect(0, 0, 1, 1, col=bgcol, lty=0)

    for (year in 1:N) {
       if (is.finite(intercept[year]) & is.finite(slope[year]))
			abline(intercept[year], slope[year])
		else abline(v=c0[year]/(c0[year]+c1[year]))
    }

    rect(-0.05, -0.05, 1.05, 0, col="white", lty=0)
    rect(-0.05, -0.05, 0, 1.05, col="white", lty=0)
    rect(-0.05, 1, 1.05, 1.05, col="white", lty=0)
    rect(1, -0.05, 1.05, 1.05, col="white", lty=0)
    box()
    par(pty="m")
    return(0)
  }

## produces temporally organized tomography plots
## (see King, 1997, A Solution to the Ecological Inference
## Problem, Princeton University Press)
##
## KQ 11/9/2002
## Modification added suggested by David Hugh-Jones 6/10/2006

#' Dynamic Tomography Plot
#'
#' dtomogplot is used to produce a tomography plot (see King, 1997) for a
#' series of temporally ordered, partially observed 2 x 2 contingency tables.
#'
#' Consider the following partially observed 2 by 2 contingency table:
#'
#' \tabular{llll}{
#'            \tab | \eqn{Y=0} \tab | \eqn{Y=1} \tab |           \cr
#'  --------- \tab ---------   \tab ---------   \tab ---------   \cr
#'  \eqn{X=0} \tab | \eqn{Y_0} \tab |           \tab | \eqn{r_0} \cr
#'  --------- \tab ---------   \tab ---------   \tab ---------   \cr
#'  \eqn{X=1} \tab | \eqn{Y_1} \tab |           \tab | \eqn{r_1} \cr
#'  --------- \tab ---------   \tab ---------   \tab ---------   \cr
#'            \tab | \eqn{c_0} \tab | \eqn{c_1} \tab | \eqn{N}
#' }
#' where \eqn{r_0}, \eqn{r_1}, \eqn{c_0}, \eqn{c_1}, and
#' \eqn{N} are non-negative integers that are observed. The interior cell
#' entries are not observed. It is assumed that \eqn{Y_0|r_0 \sim
#' \mathcal{B}inomial(r_0, p_0)} and \eqn{Y_1|r_1 \sim \mathcal{B}inomial(r_1, p_1)}.
#'
#' This function plots the bounds on the maximum likelihood estimates for (p0,
#' p1) and color codes them by the elements of time.vec.
#'
#' @param r0 An \eqn{(ntables \times 1)} vector of row sums from
#' row 0.
#'
#' @param r1 An \eqn{(ntables \times 1)} vector of row sums from
#' row 1.
#'
#' @param c0 An \eqn{(ntables \times 1)} vector of column sums
#' from column 0.
#'
#' @param c1 An \eqn{(ntables \times 1)} vector of column sums
#' from column 1.
#'
#' @param time.vec Vector of time periods that correspond to the elements of
#' \eqn{r_0}, \eqn{r_1}, \eqn{c_0}, and \eqn{c_1}.
#'
#' @param delay Time delay in seconds between the plotting of the tomography
#' lines. Setting a positive delay is useful for visualizing temporal
#' dependence.
#'
#' @param xlab The x axis label for the plot.
#'
#' @param ylab The y axis label for the plot.
#'
#' @param color.palette Color palette to be used to encode temporal patterns.
#'
#' @param bgcol The background color for the plot.
#'
#' @param ... further arguments to be passed
#'
#' @export
#'
#' @seealso \code{\link{MCMChierEI}},
#' \code{\link{MCMCdynamicEI}},\code{\link{tomogplot}}
#'
#' @references Gary King, 1997. \emph{A Solution to the Ecological Inference
#' Problem}.  Princeton: Princeton University Press.
#'
#' Jonathan C. Wakefield. 2004. ``Ecological Inference for 2 x 2 Tables.''
#' \emph{Journal of the Royal Statistical Society, Series A}. 167(3): 385445.
#'
#' Kevin Quinn. 2004. ``Ecological Inference in the Presence of Temporal
#' Dependence." In \emph{Ecological Inference: New Methodological Strategies}.
#' Gary King, Ori Rosen, and Martin A. Tanner (eds.). New York: Cambridge
#' University Press.
#'
#' @keywords hplot
#'
#' @examples
#'
#' \dontrun{
#' ## simulated data example 1
#' set.seed(3920)
#' n <- 100
#' r0 <- rpois(n, 2000)
#' r1 <- round(runif(n, 100, 4000))
#' p0.true <- pnorm(-1.5 + 1:n/(n/2))
#' p1.true <- pnorm(1.0 - 1:n/(n/4))
#' y0 <- rbinom(n, r0, p0.true)
#' y1 <- rbinom(n, r1, p1.true)
#' c0 <- y0 + y1
#' c1 <- (r0+r1) - c0
#'
#' ## plot data
#' dtomogplot(r0, r1, c0, c1, delay=0.1)
#'
#' ## simulated data example 2
#' set.seed(8722)
#' n <- 100
#' r0 <- rpois(n, 2000)
#' r1 <- round(runif(n, 100, 4000))
#' p0.true <- pnorm(-1.0 + sin(1:n/(n/4)))
#' p1.true <- pnorm(0.0 - 2*cos(1:n/(n/9)))
#' y0 <- rbinom(n, r0, p0.true)
#' y1 <- rbinom(n, r1, p1.true)
#' c0 <- y0 + y1
#' c1 <- (r0+r1) - c0
#'
#' ## plot data
#' dtomogplot(r0, r1, c0, c1, delay=0.1)
#' }
#'
"dtomogplot" <-
  function(r0, r1, c0, c1, time.vec=NA, delay=0,
           xlab="fraction of r0 in c0 (p0)",
           ylab="fraction of r1 in c0 (p1)",
           color.palette=heat.colors,
           bgcol="black", ...) {
    if (length(r0) != length(r1)){
      stop("r0 and r1 different lengths in dtomogplot().\n")
    }
    if (length(r0) != length(c0)){
      stop("r0 and c0 different lengths in dtomogplot().\n")
    }
    if (length(r0) != length(c1)){
      stop("r0 and c1 different lengths in dtomogplot().\n")
    }
    if (length(r0) != length(time.vec) & !is.na(time.vec)[1]){
      stop("r0 and time.vec different lengths in dtomogplot().\n")
    }

    intercept <-  c0/r1
    slope     <- -1 * r0/r1
    N <- length(r0)
    if (is.na(time.vec)[1])
      time.vec <- 1:N
    col.vec <- color.palette(N)

    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2,1), ncol=2), widths=c(1,lcm(w)))
    par(las=1)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar=mar)
    par(pty="m")
    plot.new()
    plot.window(xlim=c(0,1), ylim=range(time.vec), xaxs="i",
                yaxs="i")
    rect(0, time.vec[-length(time.vec)], 1, time.vec[-1], col=col.vec)
    axis(4)
    box()
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
    par(pty="s")
    plot(0:1, 0:1, type="n", main="", xlab=xlab, ylab=ylab)
    rect(0, 0, 1, 1, col=bgcol, lty=0)

    for (year in 1:N) {
      time.last <- proc.time()[3]
      time.next <- proc.time()[3]
      while ( (time.next - time.last) < delay ){
        time.next <- proc.time()[3]
      }
       if (is.finite(intercept[year]) & is.finite(slope[year]))
         abline(intercept[year], slope[year], col=col.vec[year])
       else abline(v=c0[year]/(c0[year]+c1[year]), col=col.vec[year])
    }

    rect(-0.05, -0.05, 1.05, 0, col="white", lty=0)
    rect(-0.05, -0.05, 0, 1.05, col="white", lty=0)
    rect(-0.05, 1, 1.05, 1.05, col="white", lty=0)
    rect(1, -0.05, 1.05, 1.05, col="white", lty=0)
    box()
    par(pty="m")
    return(0)
  }
