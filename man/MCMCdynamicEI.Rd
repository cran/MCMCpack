\name{MCMCdynamicEI}
\alias{MCMCdynamicEI}
\title{Markov chain Monte Carlo for Quinn's Dynamic Ecological
  Inference Model}
\description{
  MCMCdynamicEI is used to fit Quinn's dynamic ecological inference
  model for partially observed 2 x 2 contingency tables.
  }
  
\usage{
MCMCdynamicEI(r0, r1, c0, c1, burnin=5000, mcmc=200000, thin=200,
              tune=2.65316, verbose=FALSE, seed=0, W=0, nu0=1,
              delta0=0.03, nu1=1, delta1=0.03, ...)
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
    
  \item{W}{Weight (\emph{not precision}) matrix structuring the temporal
    dependence among elements of  \eqn{\theta_{0}}{theta0} and
    \eqn{\theta_{1}}{theta1}. The default value of 0 will construct a
    weight matrix that corresponds to random walk priors for
    \eqn{\theta_{0}}{theta0} and \eqn{\theta_{1}}{theta1}. The default
    assumes that the tables are equally spaced throughout time and that
    the elements of \eqn{r0}, \eqn{r1}, \eqn{c0}, and \eqn{c1} are
    temporally ordered.}

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
  An mcmc object that contains the posterior density sample.
  This object can  be summarized by functions provided by the coda package.
}

\details{
  Consider the following partially observed 2 by 2 contingency table for
  unit \eqn{t} where \eqn{t=1,\ldots,ntables}{t=1,...,ntables}:\cr
  \cr
  \tabular{llll}{
               \tab | \eqn{Y=0} \tab | \eqn{Y=1} \tab |   \cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=0} \tab | \eqn{Y_{0t}}{Y0[t]} \tab |  \tab | \eqn{r_{0t}}{r0[t]}\cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=1} \tab | \eqn{Y_{1t}}{Y1[t]} \tab |  \tab | \eqn{r_{1t}}{r1[t]}\cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \tab | \eqn{c_{0t}}{c0[t]} \tab | \eqn{c_{1t}}{c1[t]} \tab | \eqn{N_t}{N[t]}\cr    
  }
  Where \eqn{r_{0t}}{r0-t}, \eqn{r_{1t}}{r1[t]},
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
  assumed:
  \deqn{p(\theta_0|\sigma^2_0) \propto \sigma_0^{-ntables}
    \exp \left(-\frac{1}{2\sigma^2_0}
    \theta'_{0} P \theta_{0}\right)}{p(theta0|sigma^2_0) propto
    sigma^(-ntables)_0 exp(-1/(2*sigma^2_0) theta0' * P * theta0)}
  and
  \deqn{p(\theta_1|\sigma^2_1) \propto \sigma_1^{-ntables}
    \exp \left(-\frac{1}{2\sigma^2_1}
    \theta'_{1} P \theta_{1}\right)}{p(theta1|sigma^2_1) propto
    sigma^(-ntables)_1 exp(-1/(2*sigma^2_1) theta1' * P * theta1)}
  where \eqn{P_{ts}}{P[t,s]} = \eqn{-W_{ts}}{-W[t,s]} for \eqn{t} not
  equal to \eqn{s} and \eqn{P_{tt}}{P[t,t]} =
  \eqn{\sum_{s \ne t}W_{ts}}{sum(W[t,])}.
  The \eqn{\theta_{0t}}{theta0[t]} is assumed to be a priori independent of
  \eqn{\theta_{1t}}{theta1[t]} for all t.
  In addition, the
  following hyperpriors are assumed:
  \eqn{\sigma^2_0 \sim \mathcal{IG}(\nu_0/2, \delta_0/2)}{\sigma^2_0 ~
    InvGamma(nu0/2, delta0/2)}, and
  \eqn{\sigma^2_1 \sim \mathcal{IG}(\nu_1/2, \delta_1/2)}{\sigma^2_1 ~
    InvGamma(nu1/2, delta1/2)}.

  Inference centers on \eqn{p_0}{p0}, \eqn{p_1}{p1},
  \eqn{\sigma^2_0}{sigma^2_0}, and \eqn{\sigma^2_1}{sigma^2_1}.
 The Metropolis-Hastings algorithm is used to sample from the posterior
  density. 
  }
  
\references{
  Jonathan Wakefield. 2001. ``Ecological Inference for 2 x 2 Tables."
  Center for Statistics and the Social Sciences Working Paper
  no. 12. University of Washington. 

  Kevin M. Quinn. 2002. ``Ecological Inference in the Presence of
  Temporal Dependence." Paper prepared for Ecological Inference
  Conference, Harvard University, June 17-18, 2002.

  Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2003. 
  \emph{Scythe Statistical 
    Library 0.4.} \url{http://scythe.wustl.edu}.
  
  Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
  \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
  \url{http://www-fis.iarc.fr/coda/}.
}

\examples{
   \dontrun{
   r0 <- rpois(10, 300)
   r1 <- rpois(10, 200)
   c0 <- 100 + 1:10*7 + rpois(10, 30)
   c1 <- (r0+r1) - c0
   posterior <- MCMCdynamicEI(r0, r1, c0, c1,
      burnin=10000, mcmc=5000000, thin=1000)
   plot(posterior)
   summary(posterior) 
   }
}

\keyword{models}

\seealso{\code{\link{MCMCbaselineEI}}, \code{\link{MCMChierEI}},
  \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}
