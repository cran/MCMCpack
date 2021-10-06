##########################################################################
## sample from a K-dimensional two-parameter item response model with
## probit link. This is just a wrapper function that calls
## MCMCordfactanal.
##
## Andrew D. Martin
## Washington University
##
## Kevin M. Quinn
## Harvard University
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## June 8, 2003
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for K-Dimensional Item Response Theory Model
#'
#' This function generates a sample from the posterior distribution of a
#' K-dimensional item response theory (IRT) model, with standard normal priors
#' on the subject abilities (ideal points), and normal priors on the item
#' parameters.  The user supplies data and priors, and a sample from the
#' posterior distribution is returned as an mcmc object, which can be
#' subsequently analyzed with functions provided in the coda package.
#'
#' \code{MCMCirtKd} simulates from the posterior distribution using standard
#' Gibbs sampling using data augmentation (a normal draw for the subject
#' abilities, a multivariate normal draw for the item parameters, and a
#' truncated normal draw for the latent utilities). The simulation proper is
#' done in compiled C++ code to maximize efficiency.  Please consult the coda
#' documentation for a comprehensive list of functions that can be used to
#' analyze the posterior sample.
#'
#' The default number of burnin and mcmc iterations is much smaller than the
#' typical default values in MCMCpack.  This is because fitting this model is
#' extremely computationally expensive.  It does not mean that this small of a
#' number of scans will yield good estimates.  The priors of this model need to
#' be proper for identification purposes.  The user is asked to provide prior
#' means and precisions \emph{(not variances)} for the item parameters and the
#' subject parameters.
#'
#' The model takes the following form.  We assume that each subject
#' has an ability (ideal point) denoted \eqn{\theta_j} \eqn{(K \times
#' 1)}, and that each item has a difficulty parameter \eqn{\alpha_i}
#' and discrimination parameter \eqn{\beta_i} \eqn{(K \times 1)}.  The
#' observed choice by subject \eqn{j} on item \eqn{i} is the observed
#' data matrix which is \eqn{(I \times J)}.  We assume that the choice
#' is dictated by an unobserved utility:
#'
#' \deqn{z_{i,j} = - \alpha_i + \beta_i'\theta_j + \varepsilon_{i,j}}
#'
#' Where the \eqn{\varepsilon_{i,j}}s are assumed to be distributed
#' standard normal.  The parameters of interest are the subject
#' abilities (ideal points) and the item parameters.
#'
#' We assume the following priors.  For the subject abilities (ideal points) we
#' assume independent standard normal priors:
#'
#' \deqn{\theta_{j,k} \sim \mathcal{N}(0,1)}
#'
#' These cannot be changed by the user.
#' For each item parameter, we assume independent normal priors:
#'
#' \deqn{\left[\alpha_i, \beta_i \right]' \sim \mathcal{N}_{(K+1)} (b_{0,i},B_{0,i})}
#'
#' Where \eqn{B_{0,i}} is a diagonal matrix.  One can specify a
#' separate prior mean and precision for each item parameter.
#'
#' The model is identified by the constraints on the item parameters (see
#' Jackman 2001).  The user cannot place constraints on the subject abilities.
#' This identification scheme differs from that in \code{MCMCirt1d}, which uses
#' constraints on the subject abilities to identify the model.  In our
#' experience, using subject ability constraints for models in greater than one
#' dimension does not work particularly well.
#'
#' As is the case with all measurement models, make sure that you have plenty
#' of free memory, especially when storing the item parameters.
#'
#' @param datamatrix The matrix of data.  Must be 0, 1, or missing values.  It
#' is of dimensionality subjects by items.
#'
#' @param dimensions The number of dimensions in the latent space.
#'
#' @param item.constraints List of lists specifying possible equality
#'   or simple inequality constraints on the item parameters. A
#'   typical entry in the list has one of three forms:
#'   \code{rowname=list(d,c)} which will constrain the dth item
#'   parameter for the item named rowname to be equal to c,
#'   \code{rowname=list(d,"+")} which will constrain the dth item
#'   parameter for the item named rowname to be positive,
#'   and\code{rowname=list(d, "-")} which will constrain the dth item
#'   parameter for the item named rowname to be negative. If x is a
#'   matrix without row names defaults names of ``V1", ``V2", ... ,
#'   etc will be used. In a K dimensional model, the first item
#'   parameter for item \eqn{i} is the difficulty parameter
#'   (\eqn{\alpha_i}), the second item parameter is the discrimation
#'   parameter on dimension 1 (\eqn{\beta_{i,1}}), the third item
#'   parameter is the discrimation parameter on dimension 2
#'   (\eqn{\beta_{i,2}}), ..., and the (K+1)th item parameter is the
#'   discrimation parameter on dimension K (\eqn{\beta_{i,1}}).  The
#'   item difficulty parameters (\eqn{\alpha}) should generally not be
#'   constrained.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 then
#' every \code{verbose}th iteration will be printed to the screen.
#'
#' @param seed The seed for the random number generator.  If NA, the Mersenne
#' Twister generator is used with default seed 12345; if an integer is passed
#' it is used to seed the Mersenne twister.  The user can also pass a list of
#' length two to use the L'Ecuyer random number generator, which is suitable
#' for parallel computation.  The first element of the list is the L'Ecuyer
#' seed, which is a vector of length six or NA (if NA a default seed of
#' \code{rep(12345,6)} is used).  The second element of list is a positive
#' substream number. See the MCMCpack specification for more details.
#'
#' @param alphabeta.start The starting values for the \eqn{\alpha} and
#' \eqn{\beta} difficulty and discrimination parameters. If
#' \code{alphabeta.start} is set to a scalar the starting value for all
#' unconstrained item parameters will be set to that scalar. If
#' \code{alphabeta.start} is a matrix of dimension \eqn{(K+1) \times
#' items} then the \code{alphabeta.start} matrix is used as the
#' starting values (except for equality-constrained elements). If
#' \code{alphabeta.start} is set to \code{NA} (the default) then starting
#' values for unconstrained elements are set to values generated from a series
#' of proportional odds logistic regression fits, and starting values for
#' inequality constrained elements are set to either 1.0 or -1.0 depending on
#' the nature of the constraints.
#'
#' @param b0 The prior means of the \eqn{\alpha} and \eqn{\beta}
#' difficulty and discrimination parameters, stacked for all items.  If a
#' scalar is passed, it is used as the prior mean for all items.
#'
#' @param B0 The prior precisions (inverse variances) of the independent normal
#' prior on the item parameters.  Can be either a scalar or a matrix of
#' dimension \eqn{(K+1) \times items}.
#'
#' @param store.item A switch that determines whether or not to store the item
#' parameters for posterior analysis.  \emph{NOTE: In applications with many
#' items this takes an enormous amount of memory. If you have many items and
#' want to want to store the item parameters you may want to thin the chain
#' heavily}.  By default, the item parameters are not stored.
#'
#' @param store.ability A switch that determines whether or not to store the
#' subject abilities for posterior analysis. \emph{NOTE: In applications with
#' many subjects this takes an enormous amount of memory. If you have many
#' subjects and want to want to store the ability parameters you may want to
#' thin the chain heavily}. By default, the ability parameters are all stored.
#' @param drop.constant.items A switch that determines whether or not items
#' that have no variation should be deleted before fitting the model. Default =
#' TRUE.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the posterior sample.  This object can
#' be summarized by functions provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MCMCpack]{MCMCirt1d}}, \code{\link[MCMCpack]{MCMCordfactanal}}
#'
#' @references James H. Albert. 1992. ``Bayesian Estimation of Normal Ogive
#' Item Response Curves Using Gibbs Sampling." \emph{Journal of Educational
#' Statistics}.  17: 251-269.
#'
#' Joshua Clinton, Simon Jackman, and Douglas Rivers. 2004. ``The Statistical
#' Analysis of Roll Call Data."  \emph{American Political Science Review}.  98:
#' 355-370.
#'
#' Simon Jackman. 2001. ``Multidimensional Analysis of Roll Call Data via
#' Bayesian Simulation.'' \emph{Political Analysis.} 9: 227-241.
#'
#' Valen E. Johnson and James H. Albert. 1999. \emph{Ordinal Data Modeling}.
#' Springer: New York.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' Douglas Rivers.  2004.  ``Identification of Multidimensional Item-Response
#' Models."  Stanford University, typescript.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
#'    data(SupremeCourt)
#'    # note that the rownames (the item names) are "1", "2", etc
#'    posterior1 <- MCMCirtKd(t(SupremeCourt), dimensions=1,
#'                    burnin=5000, mcmc=50000, thin=10,
#'                    B0=.25, store.item=TRUE,
#'                    item.constraints=list("1"=list(2,"-")))
#'    plot(posterior1)
#'    summary(posterior1)
#'
#'
#'    data(Senate)
#'    Sen.rollcalls <- Senate[,6:677]
#'    posterior2 <- MCMCirtKd(Sen.rollcalls, dimensions=2,
#'                    burnin=5000, mcmc=50000, thin=10,
#'                    item.constraints=list(rc2=list(2,"-"), rc2=c(3,0),
#'                                          rc3=list(3,"-")),
#'                    B0=.25)
#'    plot(posterior2)
#'    summary(posterior2)
#'    }
#'
"MCMCirtKd" <-
  function(datamatrix, dimensions, item.constraints=list(),
           burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA,
           alphabeta.start = NA, b0=0, B0=0,
           store.item=FALSE, store.ability=TRUE,
           drop.constant.items=TRUE, ... ) {

    datamatrix <- as.matrix(datamatrix)

    post <- MCMCordfactanal(x=datamatrix, factors=dimensions,
                            lambda.constraints=item.constraints,
                            burnin=burnin, mcmc=mcmc, thin=thin,
                            tune=NA, verbose=verbose, seed=seed,
                            lambda.start=alphabeta.start,
                            l0=b0, L0=B0, store.lambda=store.item,
                            store.scores=store.ability,
                            drop.constantvars=drop.constant.items,
                            model="MCMCirtKd")
    return(post)
  }
