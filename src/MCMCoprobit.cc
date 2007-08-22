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
// updated to Scythe 1.0.X 7/10/2007 ADM 
//
#ifndef MCMCOPROBIT_CC
#define MCMCOPROBIT_CC

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

/* MCMCoprobit implementation.  Takes Matrix<> reference which it
 * fills with the posterior.
 */
template <typename RNGTYPE>
void MCMCoprobit_impl (rng<RNGTYPE>& stream, const int * Y,
    const Matrix<>& X, Matrix<>& beta, Matrix<>& gamma,
    const Matrix<>& b0, const Matrix<>& B0,
    const unsigned int burnin, const unsigned int mcmc, const unsigned int thin, 
    const unsigned int verbose, const double tune, Matrix<>& result) {

    // define constants and from cross-product matrices
    const unsigned int tot_iter = burnin + mcmc;  // total number of mcmc iterations
    const unsigned int nstore = mcmc / thin;      // number of draws to store
    const unsigned int k = X.cols();
    const unsigned int N = X.rows();
    const int ncat = gamma.rows() - 1;
    const Matrix<> XpX = crossprod(X);

    // storage matrix or matrices
    Matrix<> storemat(nstore, k+ncat+1);
      
    // initialize Z
    Matrix<> gamma_p = gamma;
    Matrix<> Z(N,1);
    Matrix<> Xbeta = X * beta;



    // Gibbs loop
    int count = 0;
    int accepts = 0;
    for (unsigned int iter = 0; iter < tot_iter; ++iter) {

      // [gamma | Z, beta]
      for (int i=2; i<(ncat); ++i){
	if (i==(ncat-1)){
	  gamma_p[i] = stream.rtbnorm_combo(gamma[i], ::pow(tune, 2.0), 
					    gamma_p[i-1]);
	}
	else {
	  gamma_p[i] = stream.rtnorm_combo(gamma[i], ::pow(tune, 2.0), 
					   gamma_p[i-1], 
					   gamma[i+1]);
	}
      }
      double loglikerat = 0.0;
      double loggendenrat = 0.0;
	
      // loop over observations and construct the acceptance ratio
      for (unsigned int i=0; i<N; ++i){
	if (Y[i] == ncat){
	  loglikerat = loglikerat 
	    + log(1.0  - 
		  pnorm(gamma_p[Y[i]-1] - Xbeta[i], 0.0, 1.0) ) 
	    - log(1.0 - 
		  pnorm(gamma[Y[i]-1] - Xbeta[i], 0.0, 1.0) );
	}
	else if (Y[i] == 1){
	  loglikerat = loglikerat 
	    + log(pnorm(gamma_p[Y[i]] - Xbeta[i], 0.0, 1.0)  ) 
	    - log(pnorm(gamma[Y[i]] - Xbeta[i], 0.0, 1.0) );
	}
	else{
	  loglikerat = loglikerat 
	    + log(pnorm(gamma_p[Y[i]] - Xbeta[i], 0.0, 1.0) - 
		  pnorm(gamma_p[Y[i]-1] - Xbeta[i], 0.0, 1.0) ) 
	    - log(pnorm(gamma[Y[i]] - Xbeta[i], 0.0, 1.0) - 
		  pnorm(gamma[Y[i]-1] - Xbeta[i], 0.0, 1.0) );
	}
      }
      //		  for (int j=2; j<ncat; ++j){	   
      //     		loggendenrat = loggendenrat
      //		  + log(pnorm(gamma[j+1], gamma[j], tune) - 
      //			pnorm(gamma_p[j-1], gamma[j], tune) )  
      //		  - log(pnorm(gamma_p[j+1], gamma_p[j], tune) - 
      //			pnorm(gamma[j-1], gamma_p[j], tune) );	  
      //		  }
      double logacceptrat = loglikerat + loggendenrat;
      if (stream.runif() <= exp(logacceptrat)){
	gamma = gamma_p;
	++accepts;
      }
          
	
      // [Z| gamma, beta, y] 
      //Matrix<double> Z_mean = X * beta;
      for (unsigned int i=0; i<N; ++i){
	Z[i] = stream.rtnorm_combo(Xbeta[i], 1.0, gamma[Y[i]-1], gamma[Y[i]]);
      }
      
      // [beta|Z, gamma]
      const Matrix<> XpZ = t(X) * Z;
      beta = NormNormregress_beta_draw(XpX, XpZ, b0, B0, 1.0, stream);      
      Xbeta = X * beta;
		  
		  
      // store values in matrices
      if (iter >= burnin && ((iter % thin)==0)){ 
	for (unsigned int j=0; j<k; ++j)
	  storemat(count, j) = beta[j];
	for (int j=0; j<(ncat+1); ++j)
	  storemat(count, j+k) = gamma[j];	    
	++count;
      }


      
      
      // print output to stdout
      if(verbose > 0 && iter % verbose == 0){
	Rprintf("\n\nMCMCoprobit iteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("beta = \n");
	for (unsigned int j=0; j<k; ++j)
	  Rprintf("%10.5f\n", beta[j]);
	Rprintf("Metropolis acceptance rate for gamma = %3.5f\n\n", 
		static_cast<double>(accepts) / 
		static_cast<double>(iter+1));		
      }
     
      R_CheckUserInterrupt(); // allow user interrupts           
    }
    result = storemat;
    

    Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    Rprintf("The Metropolis acceptance rate for beta was %3.5f", 
	    static_cast<double>(accepts) / static_cast<double>(tot_iter));
    Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

}


extern "C"{
  
  void MCMCoprobit(double *sampledata, const int *samplerow, 
		   const int *samplecol, const int *Y, 
		   const double *Xdata, 
		   const int *Xrow, const int *Xcol, const int *burnin, 
		   const int *mcmc, const int *thin, const double* tune,
		   const int *uselecuyer, const int *seedarray, 
		   const int *lecuyerstream, const int *verbose, 
		   const double *betadata, const int *betarow, 
		   const int *betacol, const double* gammadata, 
		   const int* gammarow, const int* gammacol,
		   const double *b0data, const int *b0row, const int *b0col, 
		   const double *B0data, const int *B0row, const int *B0col) {  
    
    // pull together Matrix objects 
    const Matrix <> X(*Xrow, *Xcol, Xdata);
    Matrix <> beta(*betarow, *betacol, 
				    betadata);
    Matrix <> gamma(*gammarow, *gammacol, 
				     gammadata);
    const Matrix <> b0(*b0row, *b0col, b0data);
    const Matrix <> B0(*B0row, *B0col, B0data);

    Matrix<> storagematrix;
    MCMCPACK_PASSRNG2MODEL(MCMCoprobit_impl, Y, X, beta, gamma, b0, B0, 
                            *burnin, *mcmc, *thin, *verbose, *tune, 
                            storagematrix);
    
    const unsigned int size = *samplerow * *samplecol;
    for (unsigned int i = 0; i < size; ++i)
       sampledata[i] = storagematrix(i);
    }
}

#endif
