##########################################################################
# mcmc2dataframe converts an mcmc object to a dataframe
#
# Depends: coda
#
# Kevin M. Quinn
# University of Washington
#
# Andrew D. Martin
# Washington University
#
# November 10, 2002
#
##########################################################################


mcmc2dataframe <- function(obj){
  if (!is.mcmc(obj))
    stop("input object not of type mcmc")

  objdf <- as.data.frame(matrix(obj, nrow(obj), ncol(obj)))
  colnames(objdf) <- varnames(obj)
  rownames(objdf) <- seq(from=start(obj), to=end(obj), by=thin(obj))
  return(objdf)
}

