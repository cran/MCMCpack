\name{MCMChierEI}
\alias{MCMChierEI}
\title{Markov Chain Monte Carlo for Wakefield's Hierarchial Ecological
  Inference Model}
\description{
  `MCMChierEI' is used to fit Wakefield's hierarchical ecological inference
  model for partially observed 2 x 2 contingency tables.
  }
  
\usage{
MCMChierEI(r0, r1, c0, c1, burnin=5000, mcmc=50000, thin=1,
           verbose=0, seed=NA,
           m0=0, M0=2.287656, m1=0, M1=2.287656, a0=0.825, b0=0.0105,
           a1=0.825, b1=0.0105, ...)
   }

\arguments{
  \item{r0}{\eqn{(ntables \times 1)}{(ntables * 1)} vector of row
  sums from row 0.}

  \item{r1}{\eqn{(ntables \times 1)}{(ntables * 1)} vector of row
  sums from row 1.}

  \item{c0}{\eqn{(ntables \times 1)}{(ntables * 1)} vector of
  column sums from column 0.}

  \item{c1}{\eqn{(ntables \times 1)}{(ntables * 1)} vector of
  column sums from column 1.}

  \item{burnin}{The number of burn-in scans for the sampler.}

  \item{mcmc}{The number of mcmc scans to be saved.}

  \item{thin}{The thinning interval used in the simulation.  The number of
    mcmc iterations must be divisible by this value.}

  \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen. If \code{verbose} is greater
    than 0 then every \code{verbose}th iteration will be printed to the
    screen. }

    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details.}

  \item{m0}{Prior mean of the \eqn{\mu_0}{mu0} parameter.}

  \item{M0}{Prior variance of the \eqn{\mu_0}{mu0} parameter.}

  \item{m1}{Prior mean of the \eqn{\mu_1}{mu1} parameter.}

  \item{M1}{Prior variance of the \eqn{\mu_1}{mu1} parameter.}

  \item{a0}{\code{a0/2} is the shape parameter for the inverse-gamma
    prior on the \eqn{\sigma^2_0}{sigma^2_0} parameter.}
  
  \item{b0}{\code{b0/2} is the scale parameter for the inverse-gamma
    prior on the \eqn{\sigma^2_0}{sigma^2_0} parameter.}
             
  \item{a1}{\code{a1/2} is the shape parameter for the inverse-gamma
    prior on the \eqn{\sigma^2_1}{sigma^2_1} parameter.}
  
  \item{b1}{\code{b1/2} is the scale parameter for the inverse-gamma
    prior on the \eqn{\sigma^2_1}{sigma^2_1} parameter.}
  
  \item{...}{further arguments to be passed}     
}

\value{
  An mcmc object that contains the sample from the posterior distribution.
  This object can  be summarized by functions provided by the coda package.
}

\details{
  Consider the following partially observed 2 by 2 contingency table for
  unit \eqn{t} where \eqn{t=1,\ldots,ntables}:\cr
  \cr
  \tabular{llll}{
               \tab | \eqn{Y=0} \tab | \eqn{Y=1} \tab |   \cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=0} \tab | \eqn{Y_{0t}}{Y0[t]} \tab |  \tab |\eqn{r_{0t}}{r0[t]}\cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=1} \tab | \eqn{Y_{1t}}{Y1[t]} \tab |  \tab | \eqn{r_{1t}}{r1[t]}\cr
   - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
   \tab | \eqn{c_{0t}}{c0[t]} \tab | \eqn{c_{1t}}{c1[t]} \tab | \eqn{N_t}{N[t]}\cr    
  }
  Where \eqn{r_{0t}}{r0[t]}, \eqn{r_{1t}}{r1[t]},
  \eqn{c_{0t}}{c0[t]}, \eqn{c_{1t}}{c1[t]}, and
  \eqn{N_t}{N[t]}  are non-negative integers that are
  observed. The interior cell entries are not observed. It is
  assumed that \eqn{Y_{0t}|r_{0t} \sim \mathcal{B}inomial(r_{0t},
    p_{0t})}{Y0[t]|r0[t] ~ Binomial(r0[t], p0[t])} and 
  \eqn{Y_{1t}|r_{1t} \sim \mathcal{B}inomial(r_{1t}, p_{1t})}{Y1[t]|r1[t] ~
    Binomial(r1[t],p1[t])}. Let \eqn{\theta_{0t} =
    log(p_{0t}/(1-p_{0t}))}{theta0[t] = log(p0[t]/(1-p0[t]))},
  and  \eqn{\theta_{1t} = log(p_{1t}/(1-p_{1t}))}{theta1[t] =
    log(p1[t]/(1-p1[t]))}.

   The following prior distributions are
  assumed: \eqn{\theta_{0t} \sim \mathcal{N}(\mu_0,
    \sigma^2_0)}{\theta0[t] ~ Normal(mu0, sigma^2_0)},
  \eqn{\theta_{1t} \sim \mathcal{N}(\mu_1,
    \sigma^2_1)}{\theta1[t] ~ Normal(mu1, sigma^2_1)}.
  \eqn{\theta_{0t}}{theta0[t]} is assumed to be a priori independent of
  \eqn{\theta_{1t}}{theta1[t]} for all t.
  In addition, we assume the
  following hyperpriors:
  \eqn{\mu_0 \sim \mathcal{N}(m_0,
    M_0)}{mu0 ~ Normal(m0, M0)},
  \eqn{\mu_1 \sim \mathcal{N}(m_1,
    M_1)}{mu1 ~ Normal(m1,
    M1)},
  \eqn{\sigma^2_0 \sim \mathcal{IG}(a_0/2, b_0/2)}{\sigma^2_0 ~
    InvGamma(a0/2, b0/2)}, and
  \eqn{\sigma^2_1 \sim \mathcal{IG}(a_1/2, b_1/2)}{\sigma^2_1 ~
    InvGamma(a1/2, b1/2)}.

   The default priors have been chosen to make the implied prior
   distribution for \eqn{p_{0}}{p0} and \eqn{p_{1}}{p1}
   \emph{approximately} uniform on (0,1). 
  
   Inference centers on \eqn{p_0}{p0}, \eqn{p_1}{p1}, \eqn{\mu_0}{mu0},
   \eqn{\mu_1}{mu1}, \eqn{\sigma^2_0}{sigma^2_0}, and
   \eqn{\sigma^2_1}{sigma^2_1}.
   Univariate slice sampling (Neal, 2003) along with Gibbs sampling is
   used to sample from the posterior distribution.

   See Section 5.4 of Wakefield (2003) for discussion of the priors used
   here. \code{MCMChierEI} departs from the Wakefield model in that the
   \code{mu0} and \code{mu1} are here assumed to be drawn from
   independent normal distributions whereas Wakefield assumes they are
   drawn from logistic distributions. 
}
  
\references{
Jonathan C. Wakefield.  2003. ``Ecological inference for 2x2 tables."  Read 
before the Royal Statistical Society, on November 12th, 2003.

   Radford Neal. 2003. ``Slice Sampling" (with discussion). \emph{Annals of
   Statistics}, 31: 705-767. 
   
   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004. 
    \emph{Scythe Statistical 
   Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}

\examples{
   \dontrun{
## simulated data example 
set.seed(3920)
n <- 100
r0 <- round(runif(n, 400, 1500))
r1 <- round(runif(n, 100, 4000))
p0.true <- pnorm(rnorm(n, m=0.5, s=0.25))
p1.true <- pnorm(rnorm(n, m=0.0, s=0.10))
y0 <- rbinom(n, r0, p0.true)
y1 <- rbinom(n, r1, p1.true)
c0 <- y0 + y1
c1 <- (r0+r1) - c0

## plot data
tomogplot(r0, r1, c0, c1)

## fit exchangeable hierarchical model
post <- MCMChierEI(r0,r1,c0,c1, mcmc=40000, thin=5, verbose=100,
                    seed=list(NA, 1))

p0meanHier <- colMeans(post)[1:n]
p1meanHier <- colMeans(post)[(n+1):(2*n)]

## plot truth and posterior means
pairs(cbind(p0.true, p0meanHier, p1.true, p1meanHier))
   }
}

\keyword{models}

\seealso{\code{\link{MCMCdynamicEI}},
  \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}

