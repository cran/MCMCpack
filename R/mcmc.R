# Taken almost verbatim from the coda package of Plummer, Best,
# Cowles, and Vines
mcmc2 <- function (data = NA, start = 1, end = numeric(0), thin = 1) 
{
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
        stop("Invalid start")
    if (length(end) > 1) 
        stop("Invalid end")
    if (length(thin) != 1) 
        stop("Invalid thin")
    if (missing(end)) 
        end <- start + (niter - 1) * thin
    else if (missing(start)) 
        start <- end - (niter - 1) * thin
    nobs <- floor((end - start)/thin + 1.0) # only change from coda mcmc()
    if (niter < nobs) 
        stop("Start, end and thin incompatible with data")
    else {
        end <- start + thin * (nobs - 1)
        if (nobs < niter) 
            data <- data[1:nobs, , drop = FALSE]
    }
    attr(data, "mcpar") <- c(start, end, thin)
    attr(data, "class") <- "mcmc"
    data
}
