// MCMCordfactanal.cc is C++ code to estimate an ordinal data 
// factor analysis model
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
// revised version of older MCMCordfactanal 
// 7/16/2004 KQ
// updated to new version of Scythe ADM 7/24/2004
// fixed a bug pointed out by Alexander Raach 1/16/2005 KQ
//

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


extern "C"{

// function called by R to fit model
void
ordfactanalpost (double* sampledata, const int* samplerow, 
		 const int* samplecol,
		 const int* Xdata, const int* Xrow, const int* Xcol,
		 const int* burnin, const int* mcmc,  const int* thin,
		 const double* tune, const int *lecuyer, const int *seedarray,
		 const int *lecuyerstream, const int* verbose, 
		 const double* Lamstartdata, const int* Lamstartrow, 
		 const int* Lamstartcol, 
		 const double* gamdata, const int* gamrow, const int* gamcol,
		 const int* ncatdata, const int* ncatrow, const int* ncatcol,
		 const double* Lameqdata, const int* Lameqrow, 
		 const int* Lameqcol,
		 const double* Lamineqdata, const int* Lamineqrow, 
		 const int* Lamineqcol,
		 const double* Lampmeandata, const int* Lampmeanrow, 
		 const int* Lampmeancol,
		 const double* Lampprecdata, const int* Lampprecrow,
		 const int* Lamppreccol, const int* storelambda,
		 const int* storescores,
		 int* acceptsdata, const int* acceptsrow, 
		 const int* acceptscol, const int* outswitch
		 ) {

  // put together matrices
  const Matrix<int> X = r2scythe(*Xrow, *Xcol, Xdata);
  Matrix<double> Lambda = r2scythe(*Lamstartrow, *Lamstartcol, Lamstartdata);
  Matrix<double> gamma = r2scythe(*gamrow, *gamcol, gamdata);
  const Matrix<double> ncateg = r2scythe(*ncatrow, *ncatcol, ncatdata);
  const Matrix<double> Lambda_eq = r2scythe(*Lameqrow, *Lameqcol, Lameqdata);
  const Matrix<double> Lambda_ineq = r2scythe(*Lamineqrow, *Lamineqcol, 
					      Lamineqdata);
  const Matrix<double> Lambda_prior_mean = r2scythe(*Lampmeanrow, 
						    *Lampmeancol, 
						    Lampmeandata);
  const Matrix<double> Lambda_prior_prec = r2scythe(*Lampprecrow, 
						    *Lamppreccol, 
						    Lampprecdata);  
  Matrix<int> accepts = r2scythe(*acceptsrow, *acceptscol, acceptsdata);

  
  // initialize rng stream
  rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
  
  // constants 
  const int K = X.cols();  // number of manifest variables
  const int N = X.rows();  // number of observations
  const int D = Lambda.cols();  // number of factors (including constant)
  const int tot_iter = *burnin + *mcmc;  
  const int nsamp = *mcmc / *thin;
  const Matrix<double> I = eye<double>(D-1);
  const Matrix<double> Lambda_free_indic = Matrix<double>(K, D);
  for (int i=0; i<(K*D); ++i){
    if (Lambda_eq[i] == -999) Lambda_free_indic[i] = 1.0;
  }
  const Matrix<double> Psi = eye<double>(K);
  const Matrix<double> Psi_inv = eye<double>(K);
  

  // starting values for phi, Xstar, and gamma_p
  Matrix<double> phi = Matrix<double>(N, D-1);
  //Matrix<double> phi = stream->rnorm(N, D-1);
  phi = cbind(ones<double>(N,1), phi);
  Matrix<double> Xstar = Matrix<double>(N, K);

  // storage matrices (row major order)
  Matrix<double> Lambda_store;
  if (storelambda[0]==1){
    Lambda_store = Matrix<double>(nsamp, K*D);
  }
  Matrix<double> gamma_store = Matrix<double>(nsamp, *gamrow * *gamcol);
  Matrix<double> phi_store;
  if (*storescores==1){
    phi_store = Matrix<double>(nsamp, N*D);
  }
 
  ///////////////////
  // Gibbs Sampler //
  ///////////////////
  int count = 0;  
  for (int iter=0; iter < tot_iter; ++iter){

    // sample Xstar
    for (int i=0; i<N; ++i){
      Matrix<double> X_mean = Lambda * t(phi(i,_));
      for (int j=0; j<K; ++j){
	if (X(i,j) == -999){ // if missing
	  Xstar(i,j) = stream->rnorm(X_mean[j], 1.0);
	}
	else { // if not missing
	  Xstar(i,j) = stream->rtnorm_combo(X_mean[j], 1.0, 
			      gamma(X(i,j)-1, j), gamma(X(i,j), j));
	}
      }
    }


    // sample phi
    Matrix<double> Lambda_const = Lambda(_,0);
    Matrix<double> Lambda_rest = Lambda(0, 1, K-1, D-1);
    Matrix<double> phi_post_var = invpd(I + crossprod(Lambda_rest) );
    Matrix<double> phi_post_C = cholesky(phi_post_var);
    for (int i=0; i<N; ++i){
      Matrix<double> phi_post_mean = phi_post_var * 
	(t(Lambda_rest)  * (t(Xstar(i,_))-Lambda_const));
      Matrix<double> phi_samp = gaxpy(phi_post_C, stream->rnorm(D-1, 1), 
      			      phi_post_mean);
      for (int j=0; j<(D-1); ++j)
	phi(i,j+1) = phi_samp[j];
    }
        
    // sample Lambda
    NormNormfactanal_Lambda_draw(Lambda, Lambda_free_indic,
				 Lambda_prior_mean, 
				 Lambda_prior_prec,
				 phi, Xstar, Psi_inv, Lambda_ineq, 
				 D, K, stream);

    // sample gamma
    for (int j=0; j<K; ++j){ // do the sampling for each manifest var
      Matrix<double> gamma_p = gamma(_,j);
      Matrix<double> X_mean = phi * t(Lambda(j,_));
      for (int i=2; i<(ncateg[j]); ++i){
	if (i==(ncateg[j]-1)){
	  gamma_p[i] = stream->rtbnorm_combo(gamma(i,j), ::pow(tune[j], 2.0), 
				     gamma_p[i-1]);
	}
	else {
	  gamma_p[i] = stream->rtnorm_combo(gamma(i,j), ::pow(tune[j], 2.0), 
				    gamma_p[i-1], 
				    gamma(i+1, j));
	}
      }
      double loglikerat = 0.0;
      double loggendenrat = 0.0;
      
      
      // loop over observations and construct the acceptance ratio
      for (int i=0; i<N; ++i){
	if (X(i,j) != -999){
	  if (X(i,j) == ncateg[j]){
	    loglikerat = loglikerat 
	      + log(1.0  - 
		    pnorm(gamma_p[X(i,j)-1] - X_mean[i]) ) 
	      - log(1.0 - 
		    pnorm(gamma(X(i,j)-1,j) - X_mean[i]) );
	  }
	  else if (X(i,j) == 1){
	    loglikerat = loglikerat 
	      + log(pnorm(gamma_p[X(i,j)] - X_mean[i])  ) 
	      - log(pnorm(gamma(X(i,j), j) - X_mean[i]) );
	  }
	  else{
	    loglikerat = loglikerat 
	      + log(pnorm(gamma_p[X(i,j)] - X_mean[i]) - 
		    pnorm(gamma_p[X(i,j)-1] - X_mean[i]) ) 
	      - log(pnorm(gamma(X(i,j), j) - X_mean[i]) - 
		    pnorm(gamma(X(i,j)-1, j) - X_mean[i]) );
	  }
	}
      }
      for (int k=2; k<ncateg[j]; ++k){
	loggendenrat = loggendenrat 
	    + log(pnorm(gamma(k+1,j), gamma(k,j), tune[j]) - 
		  pnorm(gamma_p[k-1], gamma(k,j), tune[j]) )  
	    - log(pnorm(gamma_p[k+1], gamma_p[k], tune[j]) - 
		  pnorm(gamma(k-1,j), gamma_p[k], tune[j]) );
      }
      double logacceptrat = loglikerat + loggendenrat;
      if (stream->runif() <= exp(logacceptrat)){
	for (int i=0; i<gamrow[0]; ++i){
	  if (gamma(i,j) == 300) break;
	  gamma(i,j) = gamma_p[i];
	}
	++accepts[j];
      }
    }
  
    
    // print results to screen
    if (verbose[0] > 0 && iter % verbose[0] == 0 && *outswitch == 1){
      Rprintf("\n\nMCMCordfactanal iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("Lambda = \n");
      for (int i=0; i<K; ++i){
	for (int j=0; j<D; ++j){
	  Rprintf("%10.5f", Lambda(i,j));
	}
	Rprintf("\n");
      }
      Rprintf("\nMetropolis-Hastings acceptance rates = \n");
	for (int j = 0; j<K; ++j){
	  Rprintf("%6.2f", 
		  static_cast<double>(accepts[j]) / 
		  static_cast<double>((iter+1)));
	} 
    }
    if (verbose[0] > 0 && iter % verbose[0] == 0 && *outswitch == 2){
      Rprintf("\n\nMCMCirtKd iteration %i of %i \n", (iter+1), tot_iter);
    }
    
    
    // store results
    if ((iter >= burnin[0]) && ((iter % thin[0]==0))) {      
      // store Lambda
      if (storelambda[0]==1){
	if (*outswitch==2){
	  for(int l=0; l<K; ++l){
	    Lambda(l,0) = Lambda(l,0) * -1.0;
	  }	    
	}
	Matrix<double> Lambda_store_vec = reshape(Lambda,1,K*D);
	for (int l=0; l<K*D; ++l)
	  Lambda_store(count, l) = Lambda_store_vec[l];
      }
      
      // store gamma
      Matrix<double> gamma_store_vec = reshape(gamma, 1, gamrow[0]*gamcol[0]);
      for (int l=0; l<gamrow[0]*gamcol[0]; ++l)
	gamma_store(count, l) = gamma_store_vec[l];

      // store phi
      if (storescores[0]==1){
	Matrix<double> phi_store_vec = reshape(phi, 1, N*D);
	for (int l=0; l<N*D; ++l)
	  phi_store(count, l) = phi_store_vec[l];
      }
      count++;
    }

    // allow user interrupts
    R_CheckUserInterrupt();    

  } // end MCMC loop
  
     delete stream; // clean up random number stream  
  
  // return output
  Matrix<double> output;
  if (*storelambda == 1){
    output = cbind(Lambda_store, gamma_store);
  }
  else {
    output = gamma_store;
  }
  if(*storescores == 1) {
    output = cbind(output, phi_store);
  }
  int size = *samplerow * *samplecol;
  for (int i=0; i<size; ++i)
    sampledata[i] = output[i];

  for (int j=0; j<K; ++j)
    acceptsdata[j] = accepts[j];
  
}

}


