########## Helper Functions for MCMC Simulations ##########

# check MCMC parameters for conformability
#
# Kevin Rompala 5/1/2002

"check.parameters" <-
  function(burnin, mcmc, thin, fcn, tune=NA) {

    if(mcmc %% thin != 0) {
      cat("Gibbs interval not evenly divisible by thinning interval.\n")
      stop("Please respecify and call ", fcn, "() again.\n")
    }
    if(mcmc < 0) {
      cat("Gibbs interval negative.\n")
      stop("Please respecify and call ", fcn, "() again.\n")
    }
    if(burnin < 0) {
      cat("Burnin interval negative.\n")
      stop("Please respecify and call ", fcn, "() again.\n")
    }
    if(thin < 0) {
      cat("Thinning interval negative.\n")
      stop("Please respecify and call ", fcn, "() again.\n")
    }
    if(!is.na(tune) & tune <= 0) {
      cat("Tuning parameter negative.\n")
      stop("Please respecify and call ", fcn, "() again.\n")      
    }
    
    return(0)
  }
