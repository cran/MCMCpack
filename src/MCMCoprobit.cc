// MCMCoprobit.cc is C++ code to estimate a ordinalprobit regression 
//   model with a multivariate normal prior
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
// fixed a bug pointed out by Alexander Raach 1/16/2005 KQ
//

#include<iostream>

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

extern "C"{

  
  void MCMCoprobit(double *sampledata, const int *samplerow, 
		   const int *samplecol, const int *Y, 
		   const double *Xdata, 
		   const int *Xrow, const int *Xcol, const int *burnin, 
		   const int *mcmc, const int *thin, const double* tune,
		   const int *lecuyer, const int *seedarray, 
		   const int *lecuyerstream, const int *verbose, 
		   const double *betadata, const int *betarow, 
		   const int *betacol, const double* gammadata, 
		   const int* gammarow, const int* gammacol,
		   const double *b0data, const int *b0row, const int *b0col, 
		   const double *B0data, const int *B0row, const int *B0col) {  
    
    // pull together Matrix objects
    const Matrix <double> X = r2scythe(*Xrow, *Xcol, Xdata);
    Matrix <double> beta = r2scythe(*betarow, *betacol, 
				    betadata);
    Matrix <double> gamma = r2scythe(*gammarow, *gammacol, 
				     gammadata);
    const Matrix <double> b0 = r2scythe(*b0row, *b0col, b0data);
    const Matrix <double> B0 = r2scythe(*B0row, *B0col, B0data);

    // define constants and from cross-product matrices
    const int tot_iter = *burnin + *mcmc;  // total number of mcmc iterations
    const int nstore = *mcmc / *thin;      // number of draws to store
    const int k = X.cols();
    const int N = X.rows();
    const int ncat = gamma.rows() - 1;
    const Matrix<double> XpX = crossprod(X);

    // storage matrix or matrices
    Matrix<double> storemat(nstore, k+ncat+1);
    
    // initialize rng stream
    rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
    
    // initialize Z
    Matrix<double> gamma_p = gamma;
    Matrix<double> Z(N,1);
    Matrix<double> Xbeta = X * beta;


    // Gibbs loop
    int count = 0;
    int accepts = 0;
    for (int iter = 0; iter < tot_iter; ++iter){
      
      // [gamma | Z, beta]
      for (int i=2; i<(ncat); ++i){
	if (i==(ncat-1)){
	  gamma_p[i] = stream->rtbnorm_combo(gamma[i], ::pow(tune[0], 2.0), 
					     gamma_p[i-1]);
	}
	else {
	  gamma_p[i] = stream->rtnorm_combo(gamma[i], ::pow(tune[0], 2.0), 
					    gamma_p[i-1], 
					    gamma[i+1]);
	}
      }
      double loglikerat = 0.0;
      double loggendenrat = 0.0;
      
      // loop over observations and construct the acceptance ratio
      for (int i=0; i<N; ++i){
	if (Y[i] == ncat){
	  loglikerat = loglikerat 
	    + log(1.0  - 
		  pnorm(gamma_p[Y[i]-1] - Xbeta[i]) ) 
	    - log(1.0 - 
		  pnorm(gamma[Y[i]-1] - Xbeta[i]) );
	}
	else if (Y[i] == 1){
	  loglikerat = loglikerat 
	    + log(pnorm(gamma_p[Y[i]] - Xbeta[i])  ) 
	    - log(pnorm(gamma[Y[i]] - Xbeta[i]) );
	}
	else{
	  loglikerat = loglikerat 
	    + log(pnorm(gamma_p[Y[i]] - Xbeta[i]) - 
		  pnorm(gamma_p[Y[i]-1] - Xbeta[i]) ) 
	    - log(pnorm(gamma[Y[i]] - Xbeta[i]) - 
		  pnorm(gamma[Y[i]-1] - Xbeta[i]) );
	}
      }
      for (int j=2; j<ncat; ++j){	   
	loggendenrat = loggendenrat
	  + log(pnorm(gamma[j+1], gamma[j], tune[0]) - 
		pnorm(gamma_p[j-1], gamma[j], tune[0]) )  
	  - log(pnorm(gamma_p[j+1], gamma_p[j], tune[0]) - 
		pnorm(gamma[j-1], gamma_p[j], tune[0]) );
	  
      }
      double logacceptrat = loglikerat + loggendenrat;
      if (stream->runif() <= exp(logacceptrat)){
	gamma = gamma_p;
	++accepts;
      }
      
      
      // [Z| gamma, beta, y] 
      //Matrix<double> Z_mean = X * beta;
      for (int i=0; i<N; ++i){
	Z[i] = stream->rtnorm_combo(Xbeta[i], 1.0, gamma[Y[i]-1], gamma[Y[i]]);
      }
      
      
      // [beta|Z, gamma]
      const Matrix<double> XpZ = t(X) * Z;
      beta = NormNormregress_beta_draw(XpX, XpZ, b0, B0, 1.0, stream);      
      Xbeta = X * beta;
      
      
      // store values in matrices
      if (iter >= *burnin && ((iter % *thin)==0)){ 
	for (int j=0; j<k; ++j)
	  storemat(count, j) = beta[j];
	for (int j=0; j<(ncat+1); ++j)
	  storemat(count, j+k) = gamma[j];	    
	++count;
      }
      
      // print output to stdout
      if(*verbose == 1 && iter % 500 == 0){
	Rprintf("\n\nMCMCoprobit iteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("beta = \n");
	for (int j=0; j<k; ++j)
	  Rprintf("%10.5f\n", beta[j]);
	Rprintf("Metropolis acceptance rate for gamma = %3.5f\n\n", 
		static_cast<double>(accepts) / 
		static_cast<double>(iter+1));		
      }
     
      
      void R_CheckUserInterrupt(void); // allow user interrupts     
    }

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
