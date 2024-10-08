##########################################################################
## sample from the posterior distribution of Wakefield's hierarchical model
## for ecological inference in R using linked C++ code in Scythe
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 10/22/2002
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for Wakefield's Hierarchial Ecological Inference
#' Model
#'
#' `MCMChierEI' is used to fit Wakefield's hierarchical ecological
#' inference model for partially observed 2 x 2 contingency tables.
#'
#' Consider the following partially observed 2 by 2 contingency table
#' for unit \eqn{t} where \eqn{t=1,\ldots,ntables}:
#'
#' \tabular{llll}{
#'             \tab | \eqn{Y=0}    \tab | \eqn{Y=1}    \tab |              \cr
#'   --------- \tab ------------   \tab ------------   \tab ------------   \cr
#'   \eqn{X=0} \tab | \eqn{Y_{0t}} \tab |              \tab | \eqn{r_{0t}} \cr
#'   --------- \tab ------------   \tab ------------   \tab ------------   \cr
#'   \eqn{X=1} \tab | \eqn{Y_{1t}} \tab |              \tab | \eqn{r_{1t}} \cr
#'   --------- \tab ------------   \tab ------------   \tab ------------   \cr
#'             \tab | \eqn{c_{0t}} \tab | \eqn{c_{1t}} \tab | \eqn{N_t}
#' }
#'
#'  Where \eqn{r_{0t}}, \eqn{r_{1t}}, \eqn{c_{0t}}, \eqn{c_{1t}}, and
#' \eqn{N_t} are non-negative integers that are observed. The interior
#' cell entries are not observed. It is assumed that
#' \eqn{Y_{0t}|r_{0t} \sim \mathcal{B}inomial(r_{0t}, p_{0t})} and
#' \eqn{Y_{1t}|r_{1t} \sim \mathcal{B}inomial(r_{1t}, p_{1t})}.  Let
#' \eqn{\theta_{0t} = log(p_{0t}/(1-p_{0t}))}, and \eqn{\theta_{1t} =
#' log(p_{1t}/(1-p_{1t}))}.
#'
#' The following prior distributions are assumed: \eqn{\theta_{0t}
#' \sim \mathcal{N}(\mu_0, \sigma^2_0)}, \eqn{\theta_{1t} \sim
#' \mathcal{N}(\mu_1, \sigma^2_1)}.  \eqn{\theta_{0t}} is assumed to
#' be a priori independent of \eqn{\theta_{1t}} for all t.  In
#' addition, we assume the following hyperpriors: \eqn{\mu_0 \sim
#' \mathcal{N}(m_0, M_0)}, \eqn{\mu_1 \sim \mathcal{N}(m_1, M_1)},
#' \eqn{\sigma^2_0 \sim \mathcal{IG}(a_0/2, b_0/2)}, and
#' \eqn{\sigma^2_1 \sim \mathcal{IG}(a_1/2, b_1/2)}.
#'
#' The default priors have been chosen to make the implied prior
#' distribution for \eqn{p_{0}} and \eqn{p_{1}} \emph{approximately}
#' uniform on (0,1).
#'
#' Inference centers on \eqn{p_0}, \eqn{p_1}, \eqn{\mu_0},
#' \eqn{\mu_1}, \eqn{\sigma^2_0}, and \eqn{\sigma^2_1}.  Univariate
#' slice sampling (Neal, 2003) along with Gibbs sampling is used to
#' sample from the posterior distribution.
#'
#' See Section 5.4 of Wakefield (2003) for discussion of the priors
#' used here.  \code{MCMChierEI} departs from the Wakefield model in
#' that the \code{mu0} and \code{mu1} are here assumed to be drawn
#' from independent normal distributions whereas Wakefield assumes
#' they are drawn from logistic distributions.
#'
#' @param r0 \eqn{(ntables \times 1)} vector of row sums from row 0.
#'
#' @param r1 \eqn{(ntables \times 1)} vector of row sums from row 1.
#'
#' @param c0 \eqn{(ntables \times 1)} vector of column sums from
#'   column 0.
#'
#' @param c1 \eqn{(ntables \times 1)} vector of column sums from
#'   column 1.
#'
#' @param burnin The number of burn-in scans for the sampler.
#'
#' @param mcmc The number of mcmc scans to be saved.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' mcmc iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen. If \code{verbose} is greater than 0 then
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
#' @param m0 Prior mean of the \eqn{\mu_0} parameter.
#'
#' @param M0 Prior variance of the \eqn{\mu_0} parameter.
#'
#' @param m1 Prior mean of the \eqn{\mu_1} parameter.
#'
#' @param M1 Prior variance of the \eqn{\mu_1} parameter.
#'
#' @param a0 \code{a0/2} is the shape parameter for the inverse-gamma
#'   prior on the \eqn{\sigma^2_0} parameter.
#'
#' @param b0 \code{b0/2} is the scale parameter for the inverse-gamma
#'   prior on the \eqn{\sigma^2_0} parameter.
#'
#' @param a1 \code{a1/2} is the shape parameter for the inverse-gamma
#'   prior on the \eqn{\sigma^2_1} parameter.
#'
#' @param b1 \code{b1/2} is the scale parameter for the inverse-gamma
#'   prior on the \eqn{\sigma^2_1} parameter.
#'
#' @param ... further arguments to be passed
#'
#' @return An mcmc object that contains the sample from the posterior
#'   distribution.  This object can be summarized by functions
#'   provided by the coda package.
#'
#' @export
#'
#' @seealso \code{\link{MCMCdynamicEI}},
#' \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}
#'
#' @references Jonathan C. Wakefield. 2004. ``Ecological Inference for 2 x 2
#' Tables.'' \emph{Journal of the Royal Statistical Society, Series A}. 167(3):
#' 385445.
#'
#' Radford Neal. 2003. ``Slice Sampling" (with discussion). \emph{Annals of
#' Statistics}, 31: 705-767.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.wustl.edu.s3-website-us-east-1.amazonaws.com/}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
#' ## simulated data example
#' set.seed(3920)
#' n <- 100
#' r0 <- round(runif(n, 400, 1500))
#' r1 <- round(runif(n, 100, 4000))
#' p0.true <- pnorm(rnorm(n, m=0.5, s=0.25))
#' p1.true <- pnorm(rnorm(n, m=0.0, s=0.10))
#' y0 <- rbinom(n, r0, p0.true)
#' y1 <- rbinom(n, r1, p1.true)
#' c0 <- y0 + y1
#' c1 <- (r0+r1) - c0
#'
#' ## plot data
#' tomogplot(r0, r1, c0, c1)
#'
#' ## fit exchangeable hierarchical model
#' post <- MCMChierEI(r0,r1,c0,c1, mcmc=40000, thin=5, verbose=100,
#'                     seed=list(NA, 1))
#'
#' p0meanHier <- colMeans(post)[1:n]
#' p1meanHier <- colMeans(post)[(n+1):(2*n)]
#'
#' ## plot truth and posterior means
#' pairs(cbind(p0.true, p0meanHier, p1.true, p1meanHier))
#'    }
#'
"MCMChierEI" <-
  function(r0, r1, c0, c1, burnin=5000, mcmc=50000, thin=1,
           verbose=0, seed=NA,
           m0=0, M0=2.287656,
           m1=0, M1=2.287656,
           a0=0.825, b0=0.0105,
           a1=0.825, b1=0.0105, ...){

    # Error checking
    if (length(r0) != length(r1)){
      cat("length(r0) != length(r1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    if (length(r0) != length(c0)){
      cat("length(r0) != length(c0).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    if (length(r0) != length(c1)){
      cat("length(r0) != length(c1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    if (length(r1) != length(c0)){
      cat("length(r1) != length(c0).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    if (length(r1) != length(c1)){
      cat("length(r1) != length(c1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    if (length(c0) != length(c1)){
      cat("length(c0) != length(c1).\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    if (min((r0+r1) == (c0+c1))==0){
      cat("Rows and columns do not sum to same thing.\n")
      stop("Please check data and try MCMChierEI() again.\n")
    }

    check.mcmc.parameters(burnin, mcmc, thin)

    # seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]


    if (M0 <= 0 ){
      cat("Parameter M0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }

    if (M1 <= 0 ){
      cat("Parameter M1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }

    if (a0 <= 0 ){
      cat("Parameter a0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }

    if (a1 <= 0 ){
      cat("Parameter a1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }

    if (b0 <= 0 ){
      cat("Parameter b0 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }

    if (b1 <= 0 ){
      cat("Parameter b1 <= 0.\n")
      stop("Please respecify and try MCMChierEI() again.\n")
    }

    # setup matrix to hold output from sampling
    ntables = length(r0)
    sample <- matrix(0, mcmc/thin, ntables*2+4)

    # call C++ code to do the sampling
    C.sample <- .C("hierEI",
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
                   mu0priormean = as.double(m0),
                   mu0priorvar = as.double(M0),
                   mu1priormean = as.double(m1),
                   mu1priorvar = as.double(M1),
                   a0 = as.double(a0),
                   b0 = as.double(b0),
                   a1 = as.double(a1),
                   b1 = as.double(b1),
                   verbose = as.integer(verbose),
                   lecuyer = as.integer(lecuyer),
                   seedarray = as.integer(seed.array),
                   lecuyerstream = as.integer(lecuyer.stream),
                   PACKAGE="MCMCpack"
                   )

    sample <- matrix(C.sample$samdata, C.sample$samrow, C.sample$samcol,
                     byrow=TRUE)

    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    p0names <- paste("p0table", 1:ntables, sep="")
    p1names <- paste("p1table", 1:ntables, sep="")
    varnames(output) <- c(p0names, p1names, "mu0", "mu1", "sigma^2.0",
                          "sigma^2.1")

    attr(output, "title") <- "MCMCpack Wakefield's Hierarchical EI Model Posterior Sample"

    return(output)

  }
