// MCMCfactanal.cc is C++ code to estimate a factor analysis model
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
// revised version of older MCMCfactanal 5/11/2004 KQ
// updated to new verion of scythe 7/25/2004 ADM


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

// ADD USER-DEFINED FUNCTIONS HERE

extern "C" {
  
  // BRIEF FUNCTION DESCRIPTION
  void MCMCfactanal(double *sampledata, const int *samplerow, 
		    const int *samplecol, const double *Xdata, 
		    const int *Xrow, const int *Xcol, const int *burnin, 
		    const int *mcmc, const int *thin, const int *lecuyer,
		    const int *seedarray, const int *lecuyerstream, 
		    const int *verbose, const double *Lambdadata, 
		    const int *Lambdarow, const int *Lambdacol, 
		    const double *Psidata, const int *Psirow, 
		    const int *Psicol, const double *Lameqdata, 
		    const int *Lameqrow, const int *Lameqcol, 
		    const double *Lamineqdata, const int *Lamineqrow, 
		    const int *Lamineqcol, const double *Lampmeandata, 
		    const int *Lampmeanrow, const int *Lampmeancol, 
		    const double *Lampprecdata, const int *Lampprecrow, 
		    const int *Lamppreccol, const double *a0data, 
		    const int *a0row, const int *a0col, 
		    const double *b0data, const int *b0row, 
		    const int *b0col, const int *storescores) {
    
    // pull together Matrix objects
    const Matrix <double> X = r2scythe(*Xrow, *Xcol, Xdata);
    Matrix <double> Lambda = r2scythe(*Lambdarow, *Lambdacol, Lambdadata);
    Matrix <double> Psi = r2scythe(*Psirow, *Psicol, Psidata);
    Matrix <double> Psi_inv = invpd(Psi);
    const Matrix <double> Lambda_eq = r2scythe(*Lameqrow, *Lameqcol, 
					       Lameqdata);
    const Matrix <double> Lambda_ineq = r2scythe(*Lamineqrow, *Lamineqcol, 
						 Lamineqdata);
    const Matrix <double> Lambda_prior_mean = r2scythe(*Lampmeanrow, 
						       *Lampmeancol, 
						       Lampmeandata);
    const Matrix <double> Lambda_prior_prec = r2scythe(*Lampprecrow, 
						       *Lamppreccol, 
						       Lampprecdata);
    const Matrix <double> a0 = r2scythe(*a0row, *a0col, a0data);
    const Matrix <double> b0 = r2scythe(*b0row, *b0col, b0data);

    // initialize rng stream
    rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);    
    
    // constants 
    const int K = X.cols();  // number of manifest variables
    const int N = X.rows();  // number of observations
    const int D = Lambda.cols();  // number of factors
    const int tot_iter = *burnin + *mcmc;  
    const int nsamp = *mcmc / *thin;
    const Matrix<double> I = eye<double>(D);
    const Matrix<double> Lambda_free_indic = Matrix<double>(K, D);
    for (int i=0; i<(K*D); ++i){
      if (Lambda_eq[i] == -999) Lambda_free_indic[i] = 1.0;
    }
    
    // starting value for phi
    Matrix<double> phi = Matrix<double>(N,D);
    
    // storage matrices (row major order)
    Matrix<double> Lambda_store = Matrix<double>(nsamp, K*D);
    Matrix<double> Psi_store = Matrix<double>(nsamp, K);
    Matrix<double> phi_store;
    if (*storescores==1){
      phi_store = Matrix<double>(nsamp, N*D);
    }
    

    int count    = 0; 
    // sampling begins here  
    for (int iter=0; iter < tot_iter; ++iter){

      // sample phi
      NormNormfactanal_phi_draw(phi, I, Lambda, Psi_inv, X, N, D, stream);

      // sample Lambda
      NormNormfactanal_Lambda_draw(Lambda, Lambda_free_indic,
					     Lambda_prior_mean, 
					     Lambda_prior_prec,
					     phi, X, Psi_inv, Lambda_ineq, 
					     D, K, stream);
      // sample Psi
      NormIGfactanal_Psi_draw(Psi, X, phi, Lambda, a0, b0, K, N, stream);
      for (int i=0; i<K; ++i)
	Psi_inv(i,i) = 1.0 / Psi(i,i);

      
      // print results to screen
      if (iter % verbose[0] == 0 && verbose[0] > 0){
	Rprintf("\n\nMCMCfactanal iteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("Lambda = \n");
	for (int i=0; i<K; ++i){
	  for (int j=0; j<D; ++j){
	    Rprintf("%10.5f", Lambda(i,j));
	  }
	  Rprintf("\n");
	}
	Rprintf("diag(Psi) = \n");
	for (int i=0; i<K; ++i){
	  Rprintf("%10.5f", Psi(i,i));
	}
	Rprintf("\n");
      }
      
      // store results
      if ((iter % thin[0])==0 && iter >= burnin[0]  ) {	
	// store Lambda
	Matrix<double> Lambda_store_vec = reshape(Lambda,1,K*D);
	for (int l=0; l<K*D; ++l)
	  Lambda_store(count, l) = Lambda_store_vec[l];
	// store Psi
	for (int i=0; i<K; ++i)
	  Psi_store(count, i) = Psi(i,i);
	// stop phi
	if (*storescores==1){
	  Matrix<double> phi_store_vec = reshape(phi, 1, N*D);
	  for (int l=0; l<N*D; ++l)
	    phi_store(count, l) = phi_store_vec[l];
	}
	count++;
      }

      // allow user interrupts
      R_CheckUserInterrupt();
            
    } // end Gibbs loop
    
      delete stream; // clean up random number stream   
    
    // return output
    
    Matrix<double> output = cbind(Lambda_store, Psi_store);
    if(*storescores == 1) {
      output = cbind(output, phi_store);
    }
    
    const int size = *samplerow * *samplecol;
    for (int i=0; i<size; ++i)
      sampledata[i] = output[i];
    
  }
  
}







