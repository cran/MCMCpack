#########################################################################
# sample from teh posterior distribution of Wakefield's baseline model
# for ecological inference in R using linked C++ code in Scythe
# 
# Kevin M. Quinn
# University of Washington
#
# Andrew D. Martin
# Washington University
#
# October 22, 2002
#
##########################################################################

MCMCbaselineEI <- function(r0, r1, c0, c1, burnin=1000, mcmc=50000, thin=10,
                           tune=2.65316, verbose=FALSE, seed=0,
                           alpha0=1, beta0=1, alpha1=1, beta1=1, method="NA",
                           ...){


# Error checking
if (length(r0) != length(r1)){
  cat("ERROR: length(r0) != length(r1).\n")
  stop("Please check data and try again.\n")
}

if (length(r0) != length(c0)){
  cat("ERROR: length(r0) != length(c0).\n")
  stop("Please check data and try again.\n")
}

if (length(r0) != length(c1)){
  cat("ERROR: length(r0) != length(c1).\n")
  stop("Please check data and try again.\n")
}

if (length(r1) != length(c0)){
  cat("ERROR: length(r1) != length(c0).\n")
  stop("Please check data and try again.\n")
}

if (length(r1) != length(c1)){
  cat("ERROR: length(r1) != length(c1).\n")
  stop("Please check data and try again.\n")
}

if (length(c0) != length(c1)){
  cat("ERROR: length(c0) != length(c1).\n")
  stop("Please check data and try again.\n")
 }

if (min((r0+r1) == (c0+c1))==0){
  cat("ERROR: rows and columns do not sum to same thing.\n")
  stop("Please check data and try again.\n")
}

if (burnin < 0 ){
  cat("ERROR: parameter burnin < 0.\n")
  stop("Please respecify and try again.\n")
}

if (mcmc < 0 ){
  cat("ERROR: parameter mcmc < 0.\n")
  stop("Please respecify and try again.\n")
}

if (mcmc%%thin != 0 ){
  cat("ERROR: parameter mcmc not evenly divisible by parameter thin.\n")
  stop("Please respecify and try again.\n")
 }

if (alpha0 <= 0 ){
  cat("ERROR: parameter alpha0 <= 0.\n")
  stop("Please respecify and try again.\n")
 }

if (beta0 <= 0 ){
  cat("ERROR: parameter beta0 <= 0.\n")
  stop("Please respecify and try again.\n")
 }

if (alpha1 <= 0 ){
  cat("ERROR: parameter alpha1 <= 0.\n")
  stop("Please respecify and try again.\n")
 }

if (beta1 <= 0 ){
  cat("ERROR: parameter beta1 <= 0.\n")
  stop("Please respecify and try again.\n")
 }

if (!(method %in% c("DA", "NA")) ){
  cat("ERROR: parameter method not DA or NA.\n")
  stop("Please respecify and try again.\n")
}



# setup matrix to hold output from sampling
ntables = length(r0)
sample <- matrix(0, mcmc/thin, ntables*4)

# call C++ code to do the sampling
if (method=="NA"){
  C.sample <- .C("baselineNA",
                 samdata = as.double(sample),
                 samrow = as.integer(nrow(sample)),
                 samcol = as.integer(ncol(sample)),
                 r0 = as.double(r0),
                 r1 = as.double(r1),
                 c0 = as.double(c0),
                 c1 = as.double(c1),
                 ntables = as.integer(ntables),
                 burnin = as.integer(burnin),
                 mcmc = as.integer(mcmc),
                 thin = as.integer(thin),
                 alpha0 = as.double(alpha0),
                 beta0 = as.double(beta0),
                 alpha1 = as.double(alpha1),
                 beta1 = as.double(beta1),
                 verbose = as.integer(verbose),
                 tune = as.double(tune),
                 seed = as.integer(seed),
                 accepts = as.integer(0),
                 PACKAGE="MCMCpack"
                 )
}
else if (method=="DA"){
  C.sample <- .C("baselineDA",
                 samdata = as.double(sample),
                 samrow = as.integer(nrow(sample)),
                 samcol = as.integer(ncol(sample)),
                 r0 = as.double(r0),
                 r1 = as.double(r1),
                 c0 = as.double(c0),
                 c1 = as.double(c1),
                 ntables = as.integer(ntables),
                 burnin = as.integer(burnin),
                 mcmc = as.integer(mcmc),
                 thin = as.integer(thin),
                 alpha0 = as.double(alpha0),
                 beta0 = as.double(beta0),
                 alpha1 = as.double(alpha1),
                 beta1 = as.double(beta1),
                 verbose = as.integer(verbose),
                 tune = as.double(tune),
                 seed = as.integer(seed),
                 PACKAGE="MCMCpack"
                 )
}

sample <- matrix(C.sample$samdata, C.sample$samrow, C.sample$samcol,
                 byrow=TRUE)
if (method=="NA"){
  sample <- sample[,1:(ntables*2)]
  output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
  p0names <- paste("p0table", 1:ntables, sep="")
  p1names <- paste("p1table", 1:ntables, sep="")
  varnames(output) <- c(p0names, p1names)
  cat(" overall acceptance rate = ",
      C.sample$accepts / (C.sample$burnin+C.sample$mcmc) / ntables, "\n")

}
else{
  output <- mcmc2(data=sample, start=1, end=mcmc, thin=thin)
  p0names <- paste("p0table", 1:ntables, sep="")
  p1names <- paste("p1table", 1:ntables, sep="")
  y0names <- paste("y0table", 1:ntables, sep="")
  y1names <- paste("y1table", 1:ntables, sep="")
  varnames(output) <- c(p0names, p1names, y0names, y1names)
}

attr(output, "title") <- paste("MCMCpack Wakefield's Baseline EI Model Posterior Density Sample, method =", method)

return(output)
  
}




















