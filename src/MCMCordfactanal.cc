//////////////////////////////////////////////////////////////////////////
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
// revised version of older MCMCordfactanal 
// 7/16/2004 KQ
// updated to new version of Scythe ADM 7/24/2004
// fixed a bug pointed out by Alexander Raach 1/16/2005 KQ
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef MCMCORDFACTANAL_CC
#define MCMCORDFACTANAL_CC

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
void MCMCordfactanal_impl(rng<RNGTYPE>& stream,
                          const Matrix<int>& X, Matrix<>& Lambda,
			  Matrix<>& gamma, const Matrix<>& ncateg,
			  const Matrix<>& Lambda_eq,
			  const Matrix<>& Lambda_ineq,
			  const Matrix<>& Lambda_prior_mean,
			  const Matrix<>& Lambda_prior_prec,
                          const double* tune,
			  bool storelambda, bool storescores,
			  int outswitch, unsigned int burnin,
			  unsigned int mcmc, unsigned int thin,
			  unsigned int verbose, Matrix<int>& accepts,
			  Matrix<>& output)
{
  // constants 
  const unsigned int K = X.cols();  // number of manifest variables
  const unsigned int N = X.rows();  // number of observations
  const unsigned int D = Lambda.cols();  // # of factors (incl constant)
  const unsigned int tot_iter = burnin + mcmc;  
  const unsigned int nsamp = mcmc / thin;
  const Matrix<> I = eye<double>(D-1);
  const Matrix<bool> Lambda_free_indic(K, D);
  for (unsigned int i = 0; i < (K * D); ++i) 
    if (Lambda_eq(i) == -999) 
      Lambda_free_indic(i) = true;

  const Matrix<> Psi = eye<double>(K);
  const Matrix<> Psi_inv = eye<double>(K);

  //Rprintf("Switches are %i %i %i\n", storelambda, storescores, outswitch);

  // starting values for phi, Xstar, and gamma_p
  Matrix<> phi(N, D-1);
  //Matrix<double> phi = stream->rnorm(N, D-1);
  phi = cbind(ones<double>(N,1), phi);
  Matrix<> Xstar(N, K);

  // storage matrices (row major order)
  Matrix<> Lambda_store;
  if (storelambda){
    Lambda_store = Matrix<double>(nsamp, K*D);
  }
  Matrix<> gamma_store(nsamp, gamma.size());
  Matrix<> phi_store;
  if (storescores){
    phi_store = Matrix<>(nsamp, N*D);
  }
 
  ///////////////////
  // Gibbs Sampler //
  ///////////////////
  int count = 0;  
  for (unsigned int iter = 0; iter < tot_iter; ++iter) {

    // sample Xstar
    for (unsigned int i = 0; i < N; ++i) {
      Matrix<> X_mean = Lambda * t(phi(i,_));
      for (unsigned int j = 0; j < K; ++j) {
	if (X(i,j) == -999) { // if missing
	  Xstar(i,j) = stream.rnorm(X_mean[j], 1.0);
	} else { // if not missing
	  Xstar(i,j) = stream.rtnorm_combo(X_mean[j], 1.0, 
					   gamma(X(i,j)-1, j), gamma(X(i,j), j));
	}
      }
    }

    // sample phi
    Matrix<> Lambda_const = Lambda(_,0);
    Matrix<> Lambda_rest = Lambda(0, 1, K-1, D-1);
    Matrix<> phi_post_var = invpd(I + crossprod(Lambda_rest) );
    Matrix<> phi_post_C = cholesky(phi_post_var);
    for (unsigned int i = 0; i < N; ++i) {
      Matrix<> phi_post_mean = phi_post_var * (t(Lambda_rest)  
					       * (t(Xstar(i,_))-Lambda_const));
      Matrix<> phi_samp = gaxpy(phi_post_C, stream.rnorm(D-1, 1, 0, 1), 
				phi_post_mean);
      for (unsigned int j = 0; j < (D-1); ++j)
	phi(i,j+1) = phi_samp(j);
    }
				
    // sample Lambda
    NormNormfactanal_Lambda_draw(Lambda, Lambda_free_indic, 
				 Lambda_prior_mean, Lambda_prior_prec,
				 phi, Xstar, Psi_inv, Lambda_ineq, D, K,
				 stream);

    // sample gamma
    for (unsigned int j = 0; j < K; ++j) { 
      // do the sampling for each manifest var
      Matrix<> gamma_p = gamma(_,j);
      Matrix<> X_mean = phi * t(Lambda(j,_));
      for (unsigned int i = 2; i < (ncateg(j)); ++i) {
	if (i == (ncateg(j)-1)) {
	  gamma_p(i) = stream.rtbnorm_combo(gamma(i,j), 
					    std::pow(tune[j], 2.0), gamma_p[i-1]);
	} else {
	  gamma_p[i] = stream.rtnorm_combo(gamma(i,j), 
					   std::pow(tune[j], 2.0), gamma_p[i-1], gamma(i+1, j));
	}
      }
      double loglikerat = 0.0;
      double loggendenrat = 0.0;
			
			
      // loop over observations and construct the acceptance ratio
      for (unsigned int i = 0; i < N; ++i) {
	if (X(i,j) != -999) {
	  if (X(i,j) == ncateg(j)) {
	    loglikerat = loglikerat + 
	      log(1.0  - pnorm(gamma_p[X(i,j)-1] - X_mean[i], 0, 1) ) 
	      - log(1.0 - pnorm(gamma(X(i,j)-1,j) - X_mean[i], 0, 1) );
	  } else if (X(i,j) == 1) { 
	    loglikerat = loglikerat + 
	      log(pnorm(gamma_p[X(i,j)] - X_mean[i], 0, 1)  ) 
	      - log(pnorm(gamma(X(i,j), j) - X_mean[i], 0, 1) );
	  } else { 
	    loglikerat = loglikerat + 
	      log(pnorm(gamma_p[X(i,j)] - X_mean[i], 0, 1) 
		  - pnorm(gamma_p[X(i,j)-1] - X_mean[i], 0, 1) ) 
	      - log(pnorm(gamma(X(i,j), j) - X_mean[i], 0, 1) - 
		    pnorm(gamma(X(i,j)-1, j) - X_mean[i], 0, 1) ); 
	  }
	}
      }

      for (unsigned int k = 2; k < ncateg(j); ++k) {
	loggendenrat = loggendenrat 
	  + log(pnorm(gamma(k+1,j), gamma(k,j), tune[j]) - 
		pnorm(gamma_p[k-1], gamma(k,j), tune[j]) )  - 
	  log(pnorm(gamma_p[k+1], gamma_p[k], tune[j]) -
	      pnorm(gamma(k-1,j), gamma_p[k], tune[j]) );
      }
      double logacceptrat = loglikerat + loggendenrat;
      if (stream() <= exp(logacceptrat)) { 
	for (unsigned int i = 0; i < gamma.rows(); ++i) {
	  if (gamma(i,j) == 300) break;
	  gamma(i,j) = gamma_p[i];
	}
	++accepts(j);
      }
    }
	
		
    // print results to screen
    if (verbose > 0 && iter % verbose == 0 && outswitch == 1) {
      Rprintf("\n\nMCMCordfactanal iteration %i of %i \n", (iter+1),
	      tot_iter);
      Rprintf("Lambda = \n");
      for (unsigned int i = 0; i < K; ++i) {
	for (unsigned int j = 0; j < D; ++j) {
	  Rprintf("%10.5f", Lambda(i,j));
	}
	Rprintf("\n");
      }
      Rprintf("\nMetropolis-Hastings acceptance rates = \n");
      for (unsigned int j = 0; j < K; ++j) { 
	Rprintf("%6.2f", static_cast<double>(accepts[j]) / 
		static_cast<double>((iter+1))); 
      }
    }
    if (verbose > 0 && iter % verbose == 0 && outswitch == 2) {
      Rprintf("\n\nMCMCirtKd iteration %i of %i \n", (iter+1),
	      tot_iter);
    }
		
    // store results
    if ((iter >= burnin) && ((iter % thin==0))) {      
      // store Lambda
      if (storelambda) { 
        rmview(Lambda_store(count, _)) = Lambda;
      }
			
      // store gamma
      //Matrix<> gamma_store_vec = reshape(gamma, 1, gamma.size());
      //for (unsigned int l = 0; l < gamma.size(); ++l) 
      //	gamma_store(count, l) = gamma_store_vec(l);
      rmview(gamma_store(count, _)) = gamma;

      // store phi
      if (storescores) {
	//Matrix<> phi_store_vec = reshape(phi, 1, N*D);
	//for (unsigned int l = 0; l < N * D; ++l)
	//	phi_store(count, l) = phi_store_vec(l);
        rmview(phi_store(count, _)) = phi;
      }
      count++;
    }

    // allow user interrupts
    R_CheckUserInterrupt();    

  } // end MCMC loop

  if (storelambda) {
    output = cbind(Lambda_store, gamma_store);
  } else {
    output = gamma_store;
  }
  if(storescores) {
    output = cbind(output, phi_store);
  }
}

extern "C"{

  // function called by R to fit model
  void
  ordfactanalpost (double* sampledata, const int* samplerow, 
		   const int* samplecol,
		   const int* Xdata, const int* Xrow, const int* Xcol,
		   const int* burnin, const int* mcmc,  const int* thin,
		   const double* tune, const int *uselecuyer, 
		   const int *seedarray,
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
		   const int* acceptscol, const int* outswitch) 
  {

    // put together matrices
    const Matrix<int> X(*Xrow, *Xcol, Xdata);
    Matrix<> Lambda(*Lamstartrow, *Lamstartcol, Lamstartdata);
    Matrix<> gamma(*gamrow, *gamcol, gamdata);
    const Matrix<> ncateg(*ncatrow, *ncatcol, ncatdata);
    const Matrix<> Lambda_eq(*Lameqrow, *Lameqcol, Lameqdata);
    const Matrix<> Lambda_ineq(*Lamineqrow, *Lamineqcol, Lamineqdata);
    const Matrix<> Lambda_prior_mean(*Lampmeanrow, *Lampmeancol, 
				     Lampmeandata);
    const Matrix<> Lambda_prior_prec(*Lampprecrow, *Lamppreccol,
				     Lampprecdata);  
    Matrix<int> accepts(*acceptsrow, *acceptscol, acceptsdata);

			
			
    // return output
    Matrix<double> output;
    MCMCPACK_PASSRNG2MODEL(MCMCordfactanal_impl, X, Lambda, gamma,
			   ncateg, Lambda_eq, Lambda_ineq, Lambda_prior_mean,
			   Lambda_prior_prec, tune, *storelambda, 
			   *storescores, *outswitch,
			   *burnin, *mcmc, *thin, *verbose, accepts, output);


    for (unsigned int i = 0; i < output.size(); ++i)
      sampledata[i] = output(i);

    for (unsigned int j = 0; j < X.cols(); ++j)
      acceptsdata[j] = accepts(j);
  }
}

#endif
