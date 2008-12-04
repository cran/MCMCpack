//////////////////////////////////////////////////////////////////////////
// MCMCmixfactanal.cc is C++ code to estimate a mixed data 
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
// revised version of older MCMCordfactanal 
// 7/20/2004 KQ
// updated to new version of Scythe 7/25/2004
// fixed a bug pointed out by Alexander Raach 1/16/2005 KQ
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////



#ifndef MCMCMIXFACTANAL_CC
#define MCMCMIXFACTANAL_CC

#include <iostream>

#include "matrix.h"
#include "algorithm.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

typedef Matrix<double,Row,View> rmview;

using namespace std;
using namespace scythe;


template <typename RNGTYPE>
void MCMCmixfactanal_impl(rng<RNGTYPE>& stream,
			  const Matrix<int>& X, 
			  const Matrix<>& Xstar,
			  Matrix<>& Psi,
			  Matrix<>& Psi_inv,
			  const Matrix<>& a0,
			  const Matrix<>& b0,
			  Matrix<>& Lambda,
			  Matrix<>& gamma, const Matrix<>& ncateg,
			  const Matrix<>& Lambda_eq,
			  const Matrix<>& Lambda_ineq,
			  const Matrix<>& Lambda_prior_mean,
			  const Matrix<>& Lambda_prior_prec,
                          const double* tune,
			  bool storelambda, bool storescores,
			  unsigned int burnin,
			  unsigned int mcmc, unsigned int thin,
			  unsigned int verbose, Matrix<int>& accepts,
			  Matrix<>& output
			  ){

  
  
    // constants 
    const unsigned int K = X.cols();  // number of manifest variables
    int n_ord_ge3 = 0; // number of ordinal varibles with >= 3 categories
    for (unsigned int i=0; i<K; ++i)
      if (ncateg[i] >= 3) ++n_ord_ge3;
    const unsigned int N = X.rows();  // number of observations
    const unsigned int D = Lambda.cols(); // number of factors 
                                          // (including constant)
    const unsigned int tot_iter = burnin + mcmc;  
    const unsigned int nsamp = mcmc / thin;
    const Matrix<> I = eye<double>(D-1);
    const Matrix<> Lambda_free_indic(K, D);
    for (unsigned int i=0; i<(K*D); ++i){
      if (Lambda_eq(i) == -999) Lambda_free_indic(i) = 1.0;
    }

    // starting values for phi  and gamma_p
    Matrix<> phi(N,D-1);
    phi = cbind(ones<double>(N,1), phi);
  
    // storage matrices (row major order)
    Matrix<> Lambda_store;
    if (storelambda==1){
      Lambda_store = Matrix<>(nsamp,K*D);
    }
    Matrix<> gamma_store(nsamp,  gamma.size());
    Matrix<> phi_store;
    if (storescores==1){
      phi_store = Matrix<>(nsamp, N*D);
    }
    Matrix<> Psi_store(nsamp, K);
    

    ///////////////////
    // Gibbs Sampler //
    ///////////////////
    int count = 0;  
    for (unsigned int iter=0; iter < tot_iter; ++iter){

      // sample Xstar
      for (unsigned int i=0; i<N; ++i){
	const Matrix<> X_mean = Lambda * t(phi(i,_));
	for (unsigned int j=0; j<K; ++j){
	  if (ncateg[j] >= 2){ // ordinal data
	    if (X(i,j) == -999){ // if missing
	      Xstar(i,j) = stream.rnorm(X_mean[j], 1.0);
	    }
	    else { // if not missing
	      Xstar(i,j) = stream.rtnorm_combo(X_mean[j], 1.0, 
					       gamma(X(i,j)-1, j), 
					       gamma(X(i,j), j));
	    }
	  }
	  else { // continuous data
	    if (X(i,j) == -999){ // if missing
	      Xstar(i,j) = stream.rnorm(X_mean[j], std::sqrt(Psi(j,j)));
	    }
	  }
	}
      }

      // sample phi
      const Matrix<> Lambda_const = Lambda(_,0);
      const Matrix<> Lambda_rest = Lambda(0, 1, K-1, D-1);
      // if Psi_inv is *not* diagonal then use:
      //Matrix<double> phi_post_var = invpd(I + t(Lambda_rest) * 
      //					Psi_inv * Lambda_rest);
      // instead of the following 2 lines:
      const Matrix<> AAA = scythe::sqrt(Psi_inv) * Lambda_rest;
      const Matrix<> phi_post_var = invpd(I + crossprod(AAA)); 
      // /////////////////////////////////////////////////////
      const Matrix<> phi_post_C = cholesky(phi_post_var);
      for (unsigned int i=0; i<N; ++i){
	const Matrix<> phi_post_mean = phi_post_var * 
	  (t(Lambda_rest)  * Psi_inv * (t(Xstar(i,_))-Lambda_const));
	const Matrix<> phi_samp = gaxpy(phi_post_C, 
					stream.rnorm(D-1, 1, 0.0, 1.0), 
					phi_post_mean);
	for (unsigned int j=0; j<(D-1); ++j)
	  phi(i,j+1) = phi_samp[j];
      }
    
      // sample Lambda
      NormNormfactanal_Lambda_draw(Lambda, Lambda_free_indic,
				   Lambda_prior_mean, 
				   Lambda_prior_prec,
				   phi, Xstar, Psi_inv, Lambda_ineq, 
				   D, K, stream);

      // sample Psi (assumes diagonal Psi)
      for (unsigned int i=0; i<K; ++i){
	if (ncateg[i] < 2){ // continuous data
	  const Matrix<> epsilon = gaxpy(phi, -1*(t(Lambda(i,_))), 
					 Xstar(_,i));      
	  const Matrix<>  SSE = crossprod(epsilon);
	  const double new_a0 = (a0[i] + N)*0.5;
	  const double new_b0 = (b0[i] + SSE[0])*0.5;
	  Psi(i,i) = stream.rigamma(new_a0, new_b0);
	  Psi_inv(i,i) = 1.0 / Psi(i,i); 
	}
      }

      // sample gamma
      for (unsigned int j=0; j<K; ++j){ // do the sampling for each categ. var
	Matrix<> gamma_p = gamma(_,j);
	if (ncateg[j] <= 2){ 
	  ++accepts[j]; 
	}
	if (ncateg[j] > 2){
	  const Matrix<> X_mean = phi * t(Lambda(j,_));
	  for (int i=2; i<(ncateg[j]); ++i){
	    if (i==(ncateg[j]-1)){
	      gamma_p(i) = stream.rtbnorm_combo(gamma(i,j), 
						::pow(tune[j], 2.0), 
						gamma_p[i-1]);
	    }
	    else {
	      gamma_p(i) = stream.rtnorm_combo(gamma(i,j), 
					       ::pow(tune[j], 2.0), 
					       gamma_p[i-1], 
					       gamma(i+1, j));
	    }
	  }
	  double loglikerat = 0.0;
	  double loggendenrat = 0.0;
		
	  // loop over observations and construct the acceptance ratio
	  for (unsigned int i=0; i<N; ++i){
	    if (X(i,j) != -999){
	      if (X(i,j) == ncateg[j]){
		loglikerat = loglikerat 
		  + log(1.0  - 
			pnorm(gamma_p[X(i,j)-1] - X_mean[i], 0.0, 1.0) ) 
		  - log(1.0 - 
			pnorm(gamma(X(i,j)-1,j) - X_mean[i], 0.0, 1.0) );
	      }
	      else if (X(i,j) == 1){
		loglikerat = loglikerat 
		  + log(pnorm(gamma_p[X(i,j)] - X_mean[i], 0.0, 1.0)  ) 
		  - log(pnorm(gamma(X(i,j), j) - X_mean[i], 0.0, 1.0) );
	      }
	      else{
		loglikerat = loglikerat 
		  + log(pnorm(gamma_p[X(i,j)] - X_mean[i], 0.0, 1.0) - 
			pnorm(gamma_p[X(i,j)-1] - X_mean[i], 0.0, 1.0) ) 
		  - log(pnorm(gamma(X(i,j), j) - X_mean[i], 0.0, 1.0) - 
			pnorm(gamma(X(i,j)-1, j) - X_mean[i], 0.0, 1.0) );
	      }
	    }
	  }
	  for (unsigned int k=2; k<ncateg[j]; ++k){
	    loggendenrat = loggendenrat 
	      + log(pnorm(gamma(k+1,j), gamma(k,j), tune[j]) - 
		    pnorm(gamma_p[k-1], gamma(k,j), tune[j]) )  
	      - log(pnorm(gamma_p[k+1], gamma_p[k], tune[j]) - 
		    pnorm(gamma(k-1,j), gamma_p[k], tune[j]) );
	  }
	  double logacceptrat = loglikerat + loggendenrat;
	  
	  if (stream() <= exp(logacceptrat)){
	    for (unsigned int i=0; i<gamma.rows(); ++i){
	      if (gamma(i,j) == 300) break;
	      gamma(i,j) = gamma_p[i];
	    }
	    ++accepts[j];
	  }
	}
      }

      // print results to screen
      if (verbose > 0 && iter % verbose == 0){
	Rprintf("\n\nMCMCmixfactanal iteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("Lambda = \n");
	for (unsigned int i=0; i<K; ++i){
	  for (unsigned int j=0; j<D; ++j){
	    Rprintf("%10.5f", Lambda(i,j));
	  }
	  Rprintf("\n");
	}
	Rprintf("diag(Psi) = \n");
	for (unsigned int i=0; i<K; ++i){
	  Rprintf("%10.5f", Psi(i,i));
	}
	Rprintf("\n");
	Rprintf("\nMetropolis-Hastings acceptance rates = \n");
	for (unsigned int j = 0; j<K; ++j){
	  Rprintf("%6.2f", 
		  static_cast<double>(accepts[j]) / 
		  static_cast<double>((iter+1)));
	} 
      }
    

      // store results
      if ((iter >= burnin) && ((iter % thin==0))) {
      
	// store Lambda
	if (storelambda==1){
	  //Matrix<double> Lambda_store_vec = reshape(Lambda,1,K*D);
	  //for (int l=0; l<K*D; ++l)
	  //  Lambda_store(count, l) = Lambda_store_vec[l];
	  rmview(Lambda_store(count, _)) = Lambda;
	}
      
	// store gamma
	//Matrix<double> gamma_store_vec = reshape(gamma, 1, *gamrow* *gamcol);
	//for (int l=0; l<*gamrow* *gamcol; ++l)
	//  gamma_store(count, l) = gamma_store_vec[l];
	rmview(gamma_store(count, _)) = gamma;

	// store Psi
	for (unsigned int l=0; l<K; ++l)
	  Psi_store(count, l) = Psi(l,l);
      
	// store phi
	if (storescores==1){
	  //Matrix<double> phi_store_vec = reshape(phi, 1, N*D);
	  //for (int l=0; l<N*D; ++l)
	  //  phi_store(count, l) = phi_store_vec[l];
	  rmview(phi_store(count, _)) = phi;
	}
	count++;
      }
    
      // allow user interrupts
      R_CheckUserInterrupt();    
    
    } // end Gibbs loop
  

    //    Matrix<double> output;
    if (storelambda==1){
      output = cbind(Lambda_store, gamma_store);
    }
    else {
      output = gamma_store;
    }
    if(storescores == 1) {
      output = cbind(output, phi_store);
    }
    output = cbind(output, Psi_store);

}




extern "C"{

  // function called by R to fit model
  void
  mixfactanalpost (double* sampledata, const int* samplerow, 
		   const int* samplecol,
		   const double* Xdata, const int* Xrow, const int* Xcol,
		   const int* burnin, const int* mcmc,  const int* thin,
		   const double* tune, const int *uselecuyer, 
		   const int *seedarray,
		   const int *lecuyerstream, const int* verbose, 
		   const double* Lamstartdata, const int* Lamstartrow, 
		   const int* Lamstartcol, 
		   const double* gamdata, const int* gamrow, const int* gamcol,
		   const double* Psistartdata, 
		   const int* Psistartrow, const int* Psistartcol,
		   const int* ncatdata, const int* ncatrow, const int* ncatcol,
		   const double* Lameqdata, const int* Lameqrow, 
		   const int* Lameqcol,
		   const double* Lamineqdata, const int* Lamineqrow, 
		   const int* Lamineqcol,
		   const double* Lampmeandata, const int* Lampmeanrow, 
		   const int* Lampmeancol,
		   const double* Lampprecdata, const int* Lampprecrow,
		   const int* Lamppreccol, 
		   const double* a0data, const int* a0row, const int* a0col,
		   const double* b0data, const int* b0row, const int* b0col,
		   const int* storelambda,
		   const int* storescores,
		   int* acceptsdata, const int* acceptsrow, 
		   const int* acceptscol
		   ) {
    
    // put together matrices
    const Matrix<> Xstar(*Xrow, *Xcol, Xdata);
    const Matrix<int> X = Matrix<int>(*Xrow, *Xcol);
    for (int i=0; i<(*Xrow * *Xcol); ++i)
      X[i] = static_cast<int>(Xstar[i]);
  
    Matrix<> Lambda(*Lamstartrow, *Lamstartcol, Lamstartdata);
    Matrix<> gamma(*gamrow, *gamcol, gamdata);
    Matrix<> Psi(*Psistartrow, *Psistartcol, Psistartdata);
    Matrix<> Psi_inv = invpd(Psi);
    const Matrix<int> ncateg(*ncatrow, *ncatcol, ncatdata);
    const Matrix<> Lambda_eq(*Lameqrow, *Lameqcol, Lameqdata);
    const Matrix<> Lambda_ineq(*Lamineqrow, *Lamineqcol, Lamineqdata);
    const Matrix<> Lambda_prior_mean(*Lampmeanrow, *Lampmeancol, Lampmeandata);
    const Matrix<> Lambda_prior_prec(*Lampprecrow, *Lamppreccol, Lampprecdata);
    const Matrix <> a0(*a0row, *a0col, a0data);
    const Matrix <> b0(*b0row, *b0col, b0data);
    Matrix<int> accepts(*acceptsrow, *acceptscol, acceptsdata);
    
  
    // return output
    Matrix<double> output;
    MCMCPACK_PASSRNG2MODEL(MCMCmixfactanal_impl, X, 
			   Xstar, Psi, Psi_inv, a0, b0, 
			   Lambda, gamma,
			   ncateg, Lambda_eq, Lambda_ineq, Lambda_prior_mean,
			   Lambda_prior_prec, tune, *storelambda, 
			   *storescores,
			   *burnin, *mcmc, *thin, *verbose, accepts, output);



    const unsigned int size = (unsigned int) *samplerow * *samplecol;
    for (unsigned int i=0; i<size; ++i)
      sampledata[i] = output[i];
    
    for (unsigned int j=0; j<X.cols(); ++j)
      acceptsdata[j] = accepts[j];
    
  
  }
  
}


#endif

