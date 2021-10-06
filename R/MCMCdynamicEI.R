##########################################################################
## sample from the posterior of Quinn's dynamic ecological inference model
## in R using linked C++ code in Scythe
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 10/25/2002
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

#' Markov Chain Monte Carlo for Quinn's Dynamic Ecological Inference Model
#'
#' MCMCdynamicEI is used to fit Quinn's dynamic ecological inference
#' model for partially observed 2 x 2 contingency tables.
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
#' Where \eqn{r_{0t}}, \eqn{r_{1t}}, \eqn{c_{0t}}, \eqn{c_{1t}}, and
#' \eqn{N_t} are non-negative integers that are observed. The interior
#' cell entries are not observed. It is assumed that
#' \eqn{Y_{0t}|r_{0t} \sim \mathcal{B}inomial(r_{0t}, p_{0t})} and
#' \eqn{Y_{1t}|r_{1t} \sim \mathcal{B}inomial(r_{1t}, p_{1t})}.  Let
#' \eqn{\theta_{0t} = log(p_{0t}/(1-p_{0t}))}, and \eqn{\theta_{1t} =
#' log(p_{1t}/(1-p_{1t}))}.
#'
#' The following prior distributions are assumed:
#'
#' \deqn{p(\theta_0|\sigma^2_0) \propto \sigma_0^{-ntables} \exp \left(-\frac{1}{2\sigma^2_0} \theta'_{0} P \theta_{0}\right)}
#'
#' and
#'
#' \deqn{p(\theta_1|\sigma^2_1) \propto \sigma_1^{-ntables} \exp \left(-\frac{1}{2\sigma^2_1} \theta'_{1} P \theta_{1}\right)}
#'
#' where \eqn{P_{ts}} = \eqn{-W_{ts}} for \eqn{t} not equal to \eqn{s}
#' and \eqn{P_{tt}} = \eqn{\sum_{s \ne t}W_{ts}}.  The
#' \eqn{\theta_{0t}} is assumed to be a priori independent of
#' \eqn{\theta_{1t}} for all t.  In addition, the following
#' hyperpriors are assumed: \eqn{\sigma^2_0 \sim \mathcal{IG}(a_0/2,
#' b_0/2)}, and \eqn{\sigma^2_1 \sim \mathcal{IG}(a_1/2, b_1/2)}.
#'
#' Inference centers on \eqn{p_0}, \eqn{p_1}, \eqn{\sigma^2_0}, and
#' \eqn{\sigma^2_1}.  Univariate slice sampling (Neal, 2003) together
#' with Gibbs sampling is used to sample from the posterior
#' distribution.
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
#' @param thin The thinning interval used in the simulation.  The
#'   number of mcmc iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the
#'   progress of the sampler is printed to the screen.  If
#'   \code{verbose} is greater than 0 then every \code{verbose}th
#'   iteration will be printed to the screen.
#'
#' @param seed The seed for the random number generator.  If NA, the
#'   Mersenne Twister generator is used with default seed 12345; if an
#'   integer is passed it is used to seed the Mersenne twister.  The
#'   user can also pass a list of length two to use the L'Ecuyer
#'   random number generator, which is suitable for parallel
#'   computation.  The first element of the list is the L'Ecuyer seed,
#'   which is a vector of length six or NA (if NA a default seed of
#'   \code{rep(12345,6)} is used).  The second element of list is a
#'   positive substream number. See the MCMCpack specification for
#'   more details.
#'
#' @param W Weight (\emph{not precision}) matrix structuring the
#'   temporal dependence among elements of \eqn{\theta_{0}} and
#'   \eqn{\theta_{1}}. The default value of 0 will construct a weight
#'   matrix that corresponds to random walk priors for
#'   \eqn{\theta_{0}} and \eqn{\theta_{1}}. The default assumes that
#'   the tables are equally spaced throughout time and that the
#'   elements of \eqn{r0}, \eqn{r1}, \eqn{c0}, and \eqn{c1} are
#'   temporally ordered.
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
#' @seealso \code{\link{MCMChierEI}},
#'   \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}
#'
#' @references
#'
#' Kevin Quinn. 2004. ``Ecological Inference in the Presence of Temporal
#' Dependence." In \emph{Ecological Inference: New Methodological Strategies}.
#' Gary King, Ori Rosen, and Martin A. Tanner (eds.). New York: Cambridge
#' University Press.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Radford Neal. 2003. ``Slice Sampling" (with discussion). \emph{Annals of
#' Statistics}, 31: 705-767.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' Jonathan C. Wakefield. 2004. ``Ecological Inference for 2 x 2 Tables.''
#' \emph{Journal of the Royal Statistical Society, Series A}. 167(3): 385445.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
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
#' ## fit dynamic model
#' post1 <- MCMCdynamicEI(r0,r1,c0,c1, mcmc=40000, thin=5, verbose=100,
#'                     seed=list(NA, 1))
#'
#' ## fit exchangeable hierarchical model
#' post2 <- MCMChierEI(r0,r1,c0,c1, mcmc=40000, thin=5, verbose=100,
#'                     seed=list(NA, 2))
#'
#' p0meanDyn <- colMeans(post1)[1:n]
#' p1meanDyn <- colMeans(post1)[(n+1):(2*n)]
#' p0meanHier <- colMeans(post2)[1:n]
#' p1meanHier <- colMeans(post2)[(n+1):(2*n)]
#'
#' ## plot truth and posterior means
#' pairs(cbind(p0.true, p0meanDyn, p0meanHier, p1.true, p1meanDyn, p1meanHier))
#'
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
#'
#' ## fit dynamic model
#' post1 <- MCMCdynamicEI(r0,r1,c0,c1, mcmc=40000, thin=5, verbose=100,
#'                     seed=list(NA, 1))
#'
#' ## fit exchangeable hierarchical model
#' post2 <- MCMChierEI(r0,r1,c0,c1, mcmc=40000, thin=5, verbose=100,
#'                     seed=list(NA, 2))
#'
#' p0meanDyn <- colMeans(post1)[1:n]
#' p1meanDyn <- colMeans(post1)[(n+1):(2*n)]
#' p0meanHier <- colMeans(post2)[1:n]
#' p1meanHier <- colMeans(post2)[(n+1):(2*n)]
#'
#' ## plot truth and posterior means
#' pairs(cbind(p0.true, p0meanDyn, p0meanHier, p1.true, p1meanDyn, p1meanHier))
#'    }
#'
"MCMCdynamicEI" <-
  function(r0, r1, c0, c1, burnin=5000, mcmc=50000,
           thin=1, verbose=0, seed=NA,
           W=0, a0=0.825, b0=0.0105, a1=0.825,
           b1=0.0105, ...){

    # Error checking
    if (length(r0) != length(r1)){
      cat("length(r0) != length(r1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r0) != length(c0)){
      cat("length(r0) != length(c0).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r0) != length(c1)){
      cat("length(r0) != length(c1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r1) != length(c0)){
      cat("length(r1) != length(c0).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(r1) != length(c1)){
      cat("length(r1) != length(c1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (length(c0) != length(c1)){
      cat("length(c0) != length(c1).\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    if (min((r0+r1) == (c0+c1))==0){
      cat("Rows and columns do not sum to same thing.\n")
      stop("Please check data and try MCMCdynamicEI() again.\n")
    }

    check.mcmc.parameters(burnin, mcmc, thin)

    # seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    if (a0 <= 0 ){
      cat("Parameter a0 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }

    if (b0 <= 0 ){
      cat("Parameter b0 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }

    if (a1 <= 0 ){
      cat("Parameter a1 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }

    if (b1 <= 0 ){
      cat("Parameter b1 <= 0.\n")
      stop("Please respecify and try MCMCdynamicEI() again.\n")
    }

    ntables = length(r0)

    if (W==0){ # construct weight matrix for a simple random walk assuming
               # tables are temporally ordered and 1 time unit apart
      W <- matrix(0, ntables, ntables)
      for (i in 2:(ntables)){
        W[i,i-1] <- 1
        W[i-1,i] <- 1
      }
    }

    # setup matrix to hold output from sampling
    sample <- matrix(0, mcmc/thin, ntables*2+2)

    # call C++ code to do the sampling
    C.sample <- .C("dynamicEI",
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
                   W = as.double(W),
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
    output <- mcmc(data=sample, start=(burnin+1), end=burnin+mcmc, thin=thin)
    p0names <- paste("p0table", 1:ntables, sep="")
    p1names <- paste("p1table", 1:ntables, sep="")
    varnames(output) <- c(p0names, p1names, "sigma^2_0", "sigma^2_1")

    attr(output, "title") <- "MCMCpack Quinn's Dynamic EI Model Posterior Sample"


    return(output)

  }
