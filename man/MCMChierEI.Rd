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
           m0=0, M0=10, m1=0, M1=10, nu0=1, delta0=0.5, nu1=1,
           delta1=0.5, verbose=FALSE, tune=2.65316, seed=0, ...)
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

  \item{tune}{Tuning parameter for the Metropolis-Hasting sampling.}


  \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  Information is printed
    if TRUE.}

  \item{seed}{The seed for the random number generator.  The code uses the
  Mersenne Twister, which requires an integer as an input.  If nothing
  is provided, the Scythe default seed is used.}

  \item{m0}{Prior mean of the \eqn{\mu_0}{mu0} parameter.}

  \item{M0}{Prior variance of the \eqn{\mu_0}{mu0} parameter.}

  \item{m1}{Prior mean of the \eqn{\mu_1}{mu1} parameter.}

  \item{M1}{Prior variance of the \eqn{\mu_1}{mu1} parameter.}

  \item{nu0}{Shape parameter for the inverse-gamma prior on the
    \eqn{\sigma^2_0}{sigma^2_0} parameter.}
  
  \item{delta0}{Scale parameter for the inverse-gamma prior on the
    \eqn{\sigma^2_0}{sigma^2_0} parameter.}
             
  \item{nu1}{Shape parameter for the inverse-gamma prior on the
    \eqn{\sigma^2_1}{sigma^2_1} parameter.}
  
  \item{delta1}{Scale parameter for the inverse-gamma prior on the
    \eqn{\sigma^2_1}{sigma^2_1} parameter.}


  \item{...}{further arguments to be passed}     
}

\value{
  An \code{mcmc} object that contains the posterior density sample.
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
  \eqn{\sigma^2_0 \sim \mathcal{IG}(\nu_0/2, \delta_0/2)}{\sigma^2_0 ~
    InvGamma(nu0/2, delta0/2)}, and
  \eqn{\sigma^2_1 \sim \mathcal{IG}(\nu_1/2, \delta_1/2)}{\sigma^2_1 ~
    InvGamma(nu1/2, delta1/2)}.

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

   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2002.  \emph{Scythe Statistical 
   Library 0.3.} \url{http://scythe.wustl.edu}.
   
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
posterior <- MCMChierEI(r0,r1,c0,c1)}
}

\keyword{models}

\seealso{\code{\link{MCMCbaselineEI}}, \code{\link{MCMCdynamicEI}},
  \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}

