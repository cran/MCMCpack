//////////////////////////////////////////////////////////////////////////
// cMCMCpoisson.cc is C++ code to estimate a Poisson regression model with
//   a multivariate normal prior
//
// Andrew D. Martin
// Dept. of Political Science
// Washington University in St. Louis
// admartin@wustl.edu
//
// Kevin M. Quinn
// Dept. of Government
// Harvard University
// kevin_quinn@harvard.edu
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// updated to the new version of Scythe 7/26/2004 KQ
// updated to Scythe 1.0.X 7/7/2007 ADM
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef CMCMCPOISSON_CC
#define CMCMCPOISSON_CC

#include "MCMCrng.h"
#include "MCMCfcds.h"
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "rng.h"
#include "mersenne.h"
#include "lecuyer.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;

static double poisson_logpost(const Matrix<>& Y,
			      const Matrix<>& X,
			      const Matrix<>& beta,
			      const Matrix<>& beta_prior_mean,
			      const Matrix<>& beta_prior_prec){

  // likelihood
  const Matrix<> eta = X * beta;
  const Matrix<> mu = exp(eta);
  double loglike = 0.0;
  for (unsigned int i=0; i<Y.rows(); ++i)
    loglike += -mu[i] + Y[i] * eta[i];

  // prior
  double logprior = 0.0;
  if (beta_prior_prec(0,0) != 0){
    logprior = lndmvn(beta, beta_prior_mean, invpd(beta_prior_prec));
  }

  return (loglike + logprior);
}


/* MCMCpoisson implementation.  Takes Matrix<> reference and fills with the
 * posterior.
 */
template <typename RNGTYPE>
void MCMCpoisson_impl (rng<RNGTYPE>& stream, const Matrix<>& Y,
		       const Matrix<>& X, const Matrix<>& tune, Matrix<>& beta, const Matrix<>& b0,
		       const Matrix<>& B0,  const Matrix<>& V, unsigned int burnin,
		       unsigned int mcmc, unsigned int thin, unsigned int verbose,  Matrix<>& result) {

  // define constants
  const unsigned int tot_iter = burnin + mcmc;  // total number iterations
  const unsigned int nstore = mcmc / thin;      // number of draws to store
  const unsigned int k = X.cols();

  // storage matrix or matrices
  Matrix<> storemat(nstore, k);

  // proposal parameters
  const Matrix<> propV = tune * invpd(B0 + invpd(V)) * tune;
  const Matrix<> propC = cholesky(propV) ;

  double logpost_cur = poisson_logpost(Y, X, beta, b0, B0);

  // MCMC loop
  int count = 0;
  int accepts = 0;
  for (unsigned int iter = 0; iter < tot_iter; ++iter){

    // sample beta
    const Matrix<> beta_can = gaxpy(propC, stream.rnorm(k,1,0,1), beta);
    const double logpost_can = poisson_logpost(Y,X,beta_can, b0, B0);
    const double ratio = ::exp(logpost_can - logpost_cur);

    if (stream.runif() < ratio){
      beta = beta_can;
      logpost_cur = logpost_can;
      ++accepts;
    }

    // store values in matrices
    if (iter >= burnin && (iter % thin==0)){
      storemat(count,_) = beta;
      ++count;
    }

    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n\nMCMCpoisson iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("beta = \n");
      for (unsigned int j=0; j<k; ++j)
	Rprintf("%10.5f\n", beta[j]);
      Rprintf("Metropolis acceptance rate for beta = %3.5f\n\n",
	      static_cast<double>(accepts) /
	      static_cast<double>(iter+1));
    }

    R_CheckUserInterrupt(); // allow user interrupts
  }// end MCMC loop

  result = storemat;
  if (verbose > 0){
    Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    Rprintf("The Metropolis acceptance rate for beta was %3.5f",
	    static_cast<double>(accepts) / static_cast<double>(tot_iter));
    Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  }
}
extern "C"{
  void cMCMCpoisson(double *sampledata, const int *samplerow,
		   const int *samplecol, const double *Ydata,
		   const int *Yrow, const int *Ycol, const double *Xdata,
		   const int *Xrow, const int *Xcol, const int *burnin,
		   const int *mcmc, const int *thin, const double *tunedata,
		   const int *tunerow, const int *tunecol, const int *uselecuyer,
		   const int *seedarray, const int *lecuyerstream,
		   const int *verbose, const double *betastartdata,
		   const int *betastartrow, const int *betastartcol,
		   const double *b0data, const int *b0row, const int *b0col,
		   const double *B0data, const int *B0row, const int *B0col,
		   const double *Vdata, const int *Vrow, const int *Vcol) {

   // pull together Matrix objects
    const Matrix <> Y(*Yrow, *Ycol, Ydata);
    const Matrix <> X(*Xrow, *Xcol, Xdata);
    const Matrix <> tune(*tunerow, *tunecol, tunedata);
    Matrix <> beta(*betastartrow, *betastartcol,
      betastartdata);
    const Matrix <> b0(*b0row, *b0col, b0data);
    const Matrix <> B0(*B0row, *B0col, B0data);
    const Matrix <> V(*Vrow, *Vcol, Vdata);

    Matrix<> storagematrix;
    MCMCPACK_PASSRNG2MODEL(MCMCpoisson_impl, Y, X, tune, beta, b0, B0,
			   V, *burnin, *mcmc, *thin, *verbose, storagematrix);

   const unsigned int size = *samplerow * *samplecol;
   for (unsigned int i=0; i<size; ++i)
     sampledata[i] = storagematrix(i);
  }
}

#endif
