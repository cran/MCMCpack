##########################################################################
# sample from a K-dimensional two-parameter item response model with
# probit link. This is just a wrapper function that calls
# MCMCordfactanal.
#
# Andrew D. Martin
# Washington University
#
# Kevin M. Quinn
# Harvard University
#
# June 8, 2003
#
##########################################################################

"MCMCirtKd" <-
  function(datamatrix, dimensions, item.constraints=list(),
           burnin = 1000, mcmc = 10000,
           thin=5, verbose = FALSE, seed = 0,
           alphabeta.start = NA, b0=0, B0=0,
           store.item=FALSE, store.ability=TRUE, ... ) {

    datamatrix <- t(as.matrix(datamatrix))   
    
    post <- MCMCordfactanal(x=datamatrix, factors=dimensions,
                            lambda.constraints=item.constraints,
                            burnin=burnin, mcmc=mcmc, thin=thin,
                            tune=NA, verbose=verbose, seed=seed,
                            lambda.start=alphabeta.start,
                            l0=b0, L0=B0, store.lambda=store.item,
                            store.scores=store.ability,
                            special.case="special.case")
    return(post)
  }

