\name{mcmc2}
\alias{mcmc2}
\title{Create Markov Chain Monte Carlo Objects}

\usage{
mcmc2(data= NA, start = 1, end = numeric(0), thin = 1)
}

\arguments{
  \item{data}{A vector or  matrix of MCMC output.}

  \item{start}{The iteration number of the first observation.}

  \item{end}{The iteration number of the last observation.}

  \item{thin}{The thinning interval between consecutive observations.}
}

\description{

   The function mcmc2 is used to create a Markov chain Monte Carlo object.
   The data are taken to be a vector, or a matrix with one column per
   variable.  This works slightly differently than the \code{coda} function
   to make more efficient use of memory.

   An mcmc object may be summarized by the \code{summary} function
   and visualized with the \code{plot} function.

   MCMC objects resemble time series (\code{ts}) objects and have
   methods for the generic functions \code{time}, \code{start},
   \code{end}, \code{frequency}, and \code{window}.
}

\author{Martyn Plummer (slightly modified by Kevin Quinn)}

\note{
   The format of the mcmc class has changed between coda version 0.3
   and 0.4.  Older mcmc objects will now cause \code{is.mcmc} to
   fail with an appropriate warning message.  Obsolete mcmc objects can
   be upgraded with the \code{mcmcUpgrade} function, which is part of the
   \code{coda} package.
}

\keyword{models}


\seealso{
   \code{\link{mcmc2dataframe}}
}
