########## Functions to Manipulate mcmc Objects ##########

# mcmc2 creates mcmc objects.  note: taken almost verbatim from
# the coda package of Plummer, Best, Cowles, and Vines.
#
# KQ 11/10/2002

"mcmc2" <-
  function (data = NA, start = 1, end = numeric(0), thin = 1) {
    if (is.matrix(data)) {
      niter <- nrow(data)
      nvar <- ncol(data)
    }
    else {
      niter <- length(data)
      nvar <- 1
    }
    thin <- round(thin)
    if (length(start) > 1) 
      stop("Invalid start in mcmc2().\n")
    if (length(end) > 1) 
      stop("Invalid end in mcmc2().\n")
    if (length(thin) != 1) 
      stop("Invalid thin in mcmc2().\n")
    if (missing(end)) 
      end <- start + (niter - 1) * thin
    else if (missing(start)) 
      start <- end - (niter - 1) * thin
    nobs <- floor((end - start)/thin + 1.0) # only change from coda mcmc()
    if (niter < nobs) 
      stop("Start, end and thin incompatible with data in mcmc2().\n")
    else {
      end <- start + thin * (nobs - 1)
      if (nobs < niter) 
        data <- data[1:nobs, , drop = FALSE]
    }
    attr(data, "mcpar") <- c(start, end, thin)
    attr(data, "class") <- "mcmc"
    data
  }



# mcmc2dataframe converts an mcmc object to a dataframe (requires
# coda package)
#
# KQ 11/10/2002

"mcmc2dataframe" <-
  function(obj){
    if (!is.mcmc(obj))
      stop("Input object not of type mcmc in mcmc2dataframe().\n")

    objdf <- as.data.frame(matrix(obj, nrow(obj), ncol(obj)))
    colnames(objdf) <- varnames(obj)
    rownames(objdf) <- seq(from=start(obj), to=end(obj), by=thin(obj))
    return(objdf)
  }
