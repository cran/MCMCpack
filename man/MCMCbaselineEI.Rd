\name{MCMCbaselineEI}
\alias{MCMCbaselineEI}
\title{Markov chain Monte Carlo for Wakefield's Baseline Ecological
  Inference Model}
\description{
  MCMCbaselineEI is used to fit Wakefield's baseline ecological
  inference model for partially observed 2 x 2 contingency tables.
}
  
\usage{
MCMCbaselineEI(r0, r1, c0, c1, burnin=1000, mcmc=50000, thin=10,
               tune=2.65316, verbose=FALSE, seed=0, alpha0=1, beta0=1,
               alpha1=1, beta1=1, method="NA", ...)
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

  \item{alpha0}{alpha parameter for the beta prior on \eqn{p_0}{p0}.}

  \item{beta0}{beta parameter for the beta prior on \eqn{p_0}{p0}.}

  \item{alpha1}{alpha parameter for the beta prior on \eqn{p_1}{p1}.}

  \item{beta1}{beta parameter for the beta prior on \eqn{p_1}{p1}.}

  \item{method}{Parameter determining whether a data augmentation
    algorithm should be used on the exact posterior (``DA"), or a
    Metropolis-Hastings algorithm should be used on
    Wakefield's normal approximation to the posterior (``NA"). For tables
    with large row and column sums, the preferred method is ``NA."}

  \item{...}{further arguments to be passed}     
}

\value{
  An \code{mcmc} object that contains the posterior density sample.
  This object can  be summarized by functions provided by the coda package.
}

\details{
  Consider the following partially observed 2 by 2 contingency table:\cr
  \cr
  \tabular{llll}{
               \tab | \eqn{Y=0} \tab | \eqn{Y=1} \tab |   \cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=0} \tab | \eqn{Y_0}{Y0} \tab |  \tab | \eqn{r_0}{r0}\cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
    \eqn{X=1} \tab | \eqn{Y_1}{Y1} \tab |  \tab | \eqn{r_1}{r1}\cr
    - - - - - \tab - - - - - \tab - - - - - \tab - - - - - \cr
              \tab | \eqn{c_0}{c0} \tab | \eqn{c_1}{c1} \tab | \eqn{N}\cr    
  }

  where \eqn{r_0}{r0}, \eqn{r_1}{r1}, \eqn{c_0}{c0}, \eqn{c_1}{c1}, and
  \eqn{N}  are non-negative integers that are
  observed. The interior cell entries are not observed. It is
  assumed that \eqn{Y_0|r_0 \sim \mathcal{B}inomial(r_0,
    p_0)}{Y0|r0 ~ Binomial(r0, p0)} and
  \eqn{Y_1|r_1 \sim \mathcal{B}inomial(r_1, p_1)}{Y1|r1 ~
    Binomial(r1,p1)}. Inference centers on
  \eqn{p_0}{p0} and \eqn{p_1}{p1}.
  Wakefield's baseline model starts with the
  assumption that a priori
  \eqn{p_0 \sim \mathcal{B}eta(\alpha_0, \beta_0)}{p0 ~
    Beta(alpha0, beta0)} and
  \eqn{p_1 \sim \mathcal{B}eta(\alpha_1, \beta_1)}{p1 ~
    Beta(alpha1, beta1)}.        
  }
  
\references{
   Jonathan Wakefield. 2001. ``Ecological Inference for 2 x 2 Tables,"
   Center for Statistics and the Social Sciences Working Paper
   no. 12. University of Washington. 

   Andrew D. Martin, Kevin M. Quinn, and Daniel B. Pemstein.  2002.   
   \emph{Scythe Statistical 
   Library 0.3.} \url{http://scythe.wustl.edu}.
   
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}

\examples{
   \dontrun{posterior <- MCMCbaselineEI(300, 200, 100, 400)}
}

\keyword{models}

\seealso{\code{\link{MCMChierEI}}, \code{\link{MCMCdynamicEI}},
  \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}

