############################################################################
# produces tomography plots (see King 1997)
# 
# Kevin M. Quinn
# University of Washington
#
# Andrew D. Martin
# Washington University
#
# November 9, 2002
#
##########################################################################
tomogplot <- function(r0, r1, c0, c1, 
                       xlab="fraction of r0 in c0 (p0)",
                       ylab="fraction of r1 in c0 (p1)",
                       bgcol="white", ...){
  if (length(r0) != length(r1)){
    stop("r0 and r1 different lengths")
  }
  if (length(r0) != length(c0)){
    stop("r0 and c0 different lengths")
  }
  if (length(r0) != length(c1)){
    stop("r0 and c1 different lengths")
  }


  
  intercept <-  c0/r1
  slope     <- -1 * r0/r1
  N <- length(r0)

  par(pty="s")
  plot(0:1, 0:1, type="n", main="", xlab=xlab, ylab=ylab)
  rect(0, 0, 1, 1, col=bgcol, lty=0)
  
  for (year in 1:N){
    abline(intercept[year], slope[year])
  }

  rect(-0.05, -0.05, 1.05, 0, col="white", lty=0)
  rect(-0.05, -0.05, 0, 1.05, col="white", lty=0)
  rect(-0.05, 1, 1.05, 1.05, col="white", lty=0)
  rect(1, -0.05, 1.05, 1.05, col="white", lty=0)
  box()
  
}














  






