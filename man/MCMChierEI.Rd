\name{MCMChierEI}
\alias{MCMChierEI}
\title{Markov chain Monte Carlo for Wakefield's Hierarchial Ecological
  Inference Model}
\description{
  `MCMChierEI' is used to fit Wakefield's hierarchical ecological inference
  model for partially observed 2 x 2 contingency tables.
  }
  
\usage{
MCMChierEI(r0, r1, c0, c1, burnin=1000, mcmc=50000, thin=1,
           verbose=FALSE, tune=2.65316, seed=NA,
           m0=0, M0=10, m1=0, M1=10, a0=1, b0=0.5, a1=1,
           b1=0.5, ...)
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

  \item{tune}{Tuning parameter for the Metropolis-Hasting sampling. The
    default value has been chosen to work well across most datasets. The
  user should only change the tuning parameter if s/he has carefully
  inspected the C++ code where the sampling takes place and understands
  what s/he is doing. Note that the acceptance rate will typically be
  much lower (say around 0.01) than would ordinarily be advisable. The
  reason for this low acceptance rate is because the proposal
  distribution is evenly spread around the entire tomography line while
  the posterior is concentrated in one area of the line. While many
  candidate draws are rejected, this algorithm mixes better than
  algorithms that use more local proposal distributions.}

  \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  Information is printed
    if TRUE.}

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
  An mcmc object that contains the posterior density sample.
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

   Inference centers on \eqn{p_0}{p0}, \eqn{p_1}{p1}, \eqn{\mu_0}{mu0},
   \eqn{\mu_1}{mu1}, \eqn{\sigma^2_0}{sigma^2_0}, and
   \eqn{\sigma^2_1}{sigma^2_1}.
   The Metropolis-Hastings algorithm is used to sample from the posterior
   density. 
}
  
\references{
  Jonathan Wakefield. 2001. ``Ecological Inference for 2 x 2 Tables."
   Center for Statistics and the Social Sciences Working Paper
   no. 12. University of Washington. 

   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004. 
    \emph{Scythe Statistical 
   Library 1.0.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}

\examples{
   \dontrun{
   c0 <- rpois(5, 500)
   c1 <- c(200, 140, 250, 190, 75)
   r0 <- rpois(5, 400)
   r1 <- (c0 + c1) - r0
   posterior <- MCMChierEI(r0,r1,c0,c1, mcmc=200000, thin=50)
   plot(posterior)
   summary(posterior) 
   }
}

\keyword{models}

\seealso{\code{\link{MCMCbaselineEI}}, \code{\link{MCMCdynamicEI}},
  \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}

