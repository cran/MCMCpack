// MCMCpoisson.cc is C++ code to estimate a Poisson regression model with
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
// Copyright (C) 2004 Andrew D. Martin and Kevin M. Quinn
// 
// updated to the new version of Scythe 7/26/2004 KQ

#include <iostream>
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace SCYTHE;
using namespace std;

static double poisson_logpost(const Matrix<double>& Y, 
			      const Matrix<double>& X, 
			      const Matrix<double>& beta,
			      const Matrix<double>& beta_prior_mean,
			      const Matrix<double>& beta_prior_prec){
  
  // likelihood
  const Matrix<double> eta = X * beta;
  const Matrix<double> mu = exp(eta);
  double loglike = 0.0;
  for (int i=0; i<Y.rows(); ++i) 
    loglike += -mu[i] + Y[i] * eta[i];
  
  // prior
  double logprior = 0.0;
  if (beta_prior_prec(0,0) != 0){
    logprior = lndmvn(beta, beta_prior_mean, invpd(beta_prior_prec));
  }
  
  return (loglike + logprior);
}

extern "C"{

  
  void MCMCpoisson(double *sampledata, const int *samplerow, 
		   const int *samplecol, const double *Ydata, 
		   const int *Yrow, const int *Ycol, const double *Xdata, 
		   const int *Xrow, const int *Xcol, const int *burnin, 
		   const int *mcmc, const int *thin, const double *tunedata, 
		   const int *tunerow, const int *tunecol, const int *lecuyer, 
		   const int *seedarray, const int *lecuyerstream, 
		   const int *verbose, const double *betastartdata, 
		   const int *betastartrow, const int *betastartcol, 
		   const double *b0data, const int *b0row, const int *b0col, 
		   const double *B0data, const int *B0row, const int *B0col, 
		   const double *Vdata, const int *Vrow, const int *Vcol) {  
    
    // pull together Matrix objects
    const Matrix <double> Y = r2scythe(*Yrow, *Ycol, Ydata);
    const Matrix <double> X = r2scythe(*Xrow, *Xcol, Xdata);
    const Matrix <double> tune = r2scythe(*tunerow, *tunecol, tunedata);
    Matrix <double> beta = r2scythe(*betastartrow, *betastartcol, 
				    betastartdata);
    const Matrix <double> b0 = r2scythe(*b0row, *b0col, b0data);
    const Matrix <double> B0 = r2scythe(*B0row, *B0col, B0data);
    const Matrix <double> V = r2scythe(*Vrow, *Vcol, Vdata);
    
    // define constants
    const int tot_iter = *burnin + *mcmc;  // total number of mcmc iterations
    const int nstore = *mcmc / *thin;      // number of draws to store
    const int k = X.cols();
    
    // storage matrix or matrices
    Matrix<double> storemat(nstore, k);
    
    // initialize rng stream
    rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
    
    // proposal parameters
    const Matrix<double> propV = tune * invpd(B0 + invpd(V)) * tune;
    const Matrix<double> propC = cholesky(propV) ;

    double logpost_cur = poisson_logpost(Y, X, beta, b0, B0);
    
    // MCMC loop
    int count = 0;
    int accepts = 0;
    for (int iter = 0; iter < tot_iter; ++iter){

      // sample beta
      const Matrix<double> beta_can = gaxpy(propC, stream->rnorm(k,1), beta);
      
      const double logpost_can = poisson_logpost(Y,X,beta_can, b0, B0);
      const double ratio = ::exp(logpost_can - logpost_cur); 
      
      if (stream->runif() < ratio){
	beta = beta_can;
	logpost_cur = logpost_can;
	++accepts;
      }
      
      // store values in matrices
      if (iter >= *burnin && ((iter % *thin)==0)){ 
	for (int j = 0; j < k; j++)
	  storemat(count, j) = beta[j];
	++count;
      }
      
      // print output to stdout
      if(*verbose == 1 && iter % 500 == 0){
	Rprintf("\n\nMCMCpoisson iteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("beta = \n");
	for (int j=0; j<k; ++j)
	  Rprintf("%10.5f\n", beta[j]);
	Rprintf("Metropolis acceptance rate for beta = %3.5f\n\n", 
		static_cast<double>(accepts) / 
		static_cast<double>(iter+1));	
      }
      
      void R_CheckUserInterrupt(void); // allow user interrupts    
    }// end MCMC loop

     delete stream; // clean up random number stream
    
    // return output
    const int size = *samplerow * *samplecol;
    for (int i=0; i<size; ++i)
      sampledata[i] = storemat[i];
    
    Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    Rprintf("The Metropolis acceptance rate for beta was %3.5f", 
	    static_cast<double>(accepts) / static_cast<double>(tot_iter));
    Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    
  }
  
}

