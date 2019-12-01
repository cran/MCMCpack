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


#ifndef MCMCIRTKDHET_CC
#define MCMCIRTKDHET_CC

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
void MCMCirtKdHet_impl(rng<RNGTYPE>& stream,
                const Matrix<int>& X, Matrix<>& Lambda,
			  const Matrix<>& Lambda_eq,
			  const Matrix<>& Lambda_ineq,
			  const Matrix<>& Lambda_prior_mean,
			  const Matrix<>& Lambda_prior_prec,
			  const double sigmapriorc,
			  const double sigmapriord,
                bool storelambda, bool storescores,
			  bool storesigma, unsigned int burnin,
			  unsigned int mcmc, unsigned int thin,
			  unsigned int verbose, 
			  Matrix<>& output) {
	
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
  const double c0 = sigmapriorc;
  const double d0 = sigmapriord;

  //Rprintf("Switches are %i %i %i\n", storelambda, storescores, storesigma);

  // starting values for phi, Xstar, and sigma
  Matrix<> phi(N, D-1);
  phi = cbind(ones<double>(N,1), phi);
  Matrix<> sigma2 = ones<double>(N,1); 
  Matrix<> sigma = ones<double>(N,1);
  Matrix<> Xstar(N, K);
  double sigma_norm = 1;

  // storage matrices (row major order)
  Matrix<> Lambda_store;
  if (storelambda){
    Lambda_store = Matrix<>(nsamp, K*D);
  }
  Matrix<> phi_store;
  if (storescores){
    phi_store = Matrix<>(nsamp, N*D);
  }
  Matrix<> sigma_store;
  if (storesigma){
    sigma_store = Matrix<>(nsamp, N*1);
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
	   if (X(i,j) == -999) // if missing
		 Xstar(i,j) = stream.rnorm(X_mean[j], sigma(i,0)); 
	   if (X(i,j) == 0) // if not missing
	     Xstar(i,j) = stream.rtnorm_combo(X_mean[j], sigma2(i,0), -300, 0);
	   if (X(i,j) == 1)  // if not missing
	     Xstar(i,j) = stream.rtnorm_combo(X_mean[j], sigma2(i,0), 0, 300);
      }
    }

    // sample phi
    Matrix<> Lambda_const = Lambda(_,0);
    Matrix<> Lambda_rest = Lambda(0, 1, K-1, D-1);

    for (unsigned int i=0; i<N; ++i){
      Matrix<> phi_post_var = invpd(I + crossprod(Lambda_rest) / sigma2(i,0)  );   
      Matrix<> phi_post_C = cholesky(phi_post_var);
      Matrix<> phi_post_mean = phi_post_var * 
	(t(Lambda_rest)  * (t(Xstar(i,_))-Lambda_const)) / sigma2(i,0);  
      Matrix<> phi_samp = gaxpy(phi_post_C,  stream.rnorm(D-1, 1, 0, 1),  
      			      phi_post_mean);
      for (unsigned int j=0; j<(D-1); ++j)
	phi(i,j+1) = phi_samp[j];
    }

				
    // sample Lambda
    
 for (unsigned int i=0; i<K; ++i){
      const Matrix<bool> free_indic = t(Lambda_free_indic(i,_));
      const Matrix<bool> not_free_indic = 1 - free_indic;
      if (sumc(free_indic)[0] > 0 && 
	  sumc(not_free_indic)[0] > 0){ // both constrnd & unconstrnd
	const Matrix<> phifree_i =  t(selif(t(phi), free_indic));
	const Matrix<> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
					   free_indic); // prior mean
	const Matrix<> hold = selif(t(Lambda_prior_prec(i,_)), 
				    free_indic);
	Matrix<> sig2lamfree_inv_i = 
	  eye<double>(hold.rows());   // prior prec
	for (unsigned int j=0; j<(hold.rows()); ++j)
	  sig2lamfree_inv_i(j,j) = hold[j];
	const Matrix<> Lambdacon_i = 
	  selif(t(Lambda(i,_)), not_free_indic);
	const Matrix<> phicon_i  = t(selif(t(phi), not_free_indic));
	const Matrix<> X_i = Xstar(_,i);
	const Matrix<> newX_i = gaxpy((-1.0*phicon_i), Lambdacon_i, X_i);

	for (unsigned int R=0; R<N; ++R){  // multiply through by inverse sigma
		for (unsigned int Q=0; Q<D; ++Q){
		phifree_i(R,Q) = phifree_i(R,Q) / sigma(R,0);
		}
		newX_i(R,0) = newX_i(R,0) /sigma(R,0);
	}										
	
	const Matrix<> Lam_post_var = invpd(sig2lamfree_inv_i + 
						  Psi_inv(i,i) * crossprod(phifree_i));  
	const Matrix<> Lam_post_C = cholesky(Lam_post_var);
	const Matrix<> Lam_post_mean = Lam_post_var * 
	  (sig2lamfree_inv_i * mulamfree_i + Psi_inv(i,i) * 
	   t(phifree_i) * newX_i);       
	
	Matrix<double> Lambdafree_i = 
	  gaxpy(Lam_post_C, stream.rnorm(hold.rows(), 1, 0, 1), Lam_post_mean);
	
	// check to see if inequality constraints hold
	const Matrix<> Lambda_ineq_vec = Lambda_ineq(i,_);
	double ineq_holds = 0;
	int Lam_count = 0;
	for (unsigned int j=0; j<D; ++j){
	  if (free_indic[j]==1)
	    ineq_holds = std::min(ineq_holds, 
				  Lambda_ineq_vec[j] * 
				  Lambdafree_i[Lam_count]);
	  ++Lam_count;
	}
	while (ineq_holds < 0){
	  Lambdafree_i = 
	    gaxpy(Lam_post_C, stream.rnorm(hold.rows(), 1, 0, 1), Lam_post_mean);
	  Lam_count = 0;
	  double test = 0;
	  for (unsigned int j=0; j<D; ++j){
	    if (free_indic[j]==1){
	      Matrix<> prodcheck =
		Lambda_ineq_vec[j]*Lambdafree_i[Lam_count];
	      test = std::min(test, prodcheck[0]); 	      
	      ++Lam_count;
	    }
	  }
	  ineq_holds = test;
	}
	
	// put draw into Lambda 
	Lam_count = 0;
	for (unsigned int j=0; j<D; ++j){
	  if (free_indic[j] == 1){
	    Lambda(i,j) = Lambdafree_i[Lam_count];
	    ++Lam_count;
	  }
	}
      }
      else  if (sumc(free_indic)[0] > 0){ // just unconstrained
	const Matrix<> phifree_i =  t(selif(t(phi), free_indic));
	const Matrix<> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
						 free_indic); // prior mean
	const Matrix<> hold = selif(t(Lambda_prior_prec(i,_)), 
					  free_indic);
	Matrix<> sig2lamfree_inv_i = 
	  eye<double>(hold.rows());  // prior prec
	for (unsigned int j=0; j<hold.rows(); ++j)
	  sig2lamfree_inv_i(j,j) = hold[j];

	Matrix<> sX_i = ones<double>(N,1);
	for (unsigned int R=0; R<N; ++R){  // multiply through by inverse sigma
		for (unsigned int Q=0; Q<D; ++Q){
			phifree_i(R,Q) = phifree_i(R,Q) / sigma(R,0);
		}
		sX_i(R,0) = Xstar(R,i) / sigma(R,0);
	}			
	
	const Matrix<> Lam_post_var = invpd(sig2lamfree_inv_i + 
						  Psi_inv(i,i) * crossprod(phifree_i));    // 
						                                            
	const Matrix<> Lam_post_C = cholesky(Lam_post_var);
	const Matrix<> Lam_post_mean = Lam_post_var * 
	  (sig2lamfree_inv_i * mulamfree_i + Psi_inv(i,i) * 
	   t(phifree_i) * sX_i);          // 
	Matrix<> Lambdafree_i = 
	  gaxpy(Lam_post_C,  stream.rnorm(hold.rows(), 1, 0, 1), Lam_post_mean);
	
	// check to see if inequality constraints hold
	Matrix<> Lambda_ineq_vec = Lambda_ineq(i,_);
	double ineq_holds = 0;
	for (unsigned int j=0; j<D; ++j){
	  ineq_holds = 
	    std::min(ineq_holds, Lambda_ineq_vec[j]*Lambdafree_i[j]);
	}
	while (ineq_holds < 0){
	  Lambdafree_i = 
	    gaxpy(Lam_post_C, stream.rnorm(hold.rows(), 1, 0, 1), Lam_post_mean);
	  double test = 0;
	  for (unsigned int j=0; j<D; ++j){
	    //if (free_indic[j]==1)
	    double prodcheck = Lambda_ineq_vec[j]*Lambdafree_i[j];
	    test = std::min(test, prodcheck); 
	  }
	  ineq_holds = test;
	}
	
	// put draw into Lambda
	for (unsigned int j=0; j<D; ++j){
	  Lambda(i,j) = Lambdafree_i[j];
	}
      }	
    }      

			 
	// sample sigma 
	for (unsigned int i=0; i<N; ++i){
		Matrix<> phi_i = phi(i,_);
		Matrix<> Xstar_i = Xstar(i,_);
		const Matrix<> e = gaxpy(phi_i, (-1*t(Lambda)), Xstar_i);
		const Matrix<> SSE = crossprod (t(e)); 
		const double c_post = (c0 + K) * 0.5;
		const double d_post = (d0 + SSE[0]) * 0.5;
		sigma2(i,0) = stream.rigamma(c_post, d_post);
	}	
	
	Matrix<> sigma = sqrt(sigma2);  // store as sigma
	sigma_norm = N / sum(pow(sigma,-1));  // renormalize inverse sigma to sum to N
	for (unsigned int i=0; i<N; ++i){
		sigma(i,0) = sigma(i,0) / sigma_norm;
		sigma2(i,0) = pow(sigma(i,0),2);
	}		

		
    // print results to screen
    if (verbose > 0 && iter % verbose == 0){
	  
	  //Rprintf("xstar 22, xstar 23, xstar 24, phi, sigma = \n");
	  //for (int i=0; i<N; ++i){
	  //Rprintf("%10.3f",Xstar(i,21));
	  //Rprintf("%10.3f",Xstar(i,22));
	  //Rprintf("%10.3f",Xstar(i,23));
	  //Rprintf("%10.3f",phi(i,1));
	  //Rprintf("%10.3f",sigma(i,0));
	  //Rprintf("\n");
	  //}
	  
      //Rprintf("Lambda = \n");
      //for (unsigned int i=0; i<K; ++i){
	//for (unsigned int j=0; j<D; ++j){
	  //Rprintf("%10.5f", Lambda(i,j));
	//}
	//Rprintf("\n");
      //}
     	
	//Rprintf("\nsigma renormalization factor");
	//Rprintf("%10.3f \n",sigma_norm);
	
	Rprintf("\n\nMCMCirtKdHet iteration %i of %i \n", (iter+1), tot_iter);
	
    }
		
   // store results
    if ((iter >= burnin) && ((iter % thin==0))) {    
  
     // store Lambda
      if (storelambda==1)
	    rmview(Lambda_store(count, _)) = Lambda;
	    
     // store phi
      if (storescores==1)
	    rmview(phi_store(count, _)) = phi;
	
     // store sigma
      if (storesigma==1)
	    rmview(sigma_store(count, _)) = sigma;
      
      count++;

    }

    // allow user interrupts
    R_CheckUserInterrupt();    

  } // end MCMC loop

  if (storelambda == 1){
     output = Lambda_store;
     if (storescores == 1) 
       output = cbind(output, phi_store);  
     if (storesigma == 1) 
       output = cbind(output, sigma_store);     
  } 
  if (storelambda == 0) {
     if (storescores == 1){
       output = phi_store;
     if (storesigma == 1)
       output = cbind(output, sigma_store);
     } 
     if (storescores == 0) {
        output = sigma_store;
     }
  }

}


extern "C"{

  // function called by R to fit model
  void irtKdHetpost(double *samdata, const int *samrow, const int *samcol, const int *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *Lamstartdata, const int *Lamstartrow, const int *Lamstartcol, const double *Lameqdata, const int *Lameqrow, const int *Lameqcol, const double *Lamineqdata, const int *Lamineqrow, const int *Lamineqcol, const double *Lampmeandata, const int *Lampmeanrow, const int *Lampmeancol, const double *Lampprecdata, const int *Lampprecrow, const int *Lamppreccol, const int *storelambda, const int *storescores, const int *storesigma, const double *sigmapriorc, const double *sigmapriord) {

    // put together matrices
    const Matrix<int> X(*Xrow, *Xcol, Xdata);
    Matrix<> Lambda(*Lamstartrow, *Lamstartcol, Lamstartdata);
    const Matrix<> Lambda_eq(*Lameqrow, *Lameqcol, Lameqdata);
    const Matrix<> Lambda_ineq(*Lamineqrow, *Lamineqcol, Lamineqdata);
    const Matrix<> Lambda_prior_mean(*Lampmeanrow, *Lampmeancol, 
				     Lampmeandata);
    const Matrix<> Lambda_prior_prec(*Lampprecrow, *Lamppreccol,
				     Lampprecdata);  
    
			
			
    // return output
    Matrix<double> output;
    MCMCPACK_PASSRNG2MODEL(MCMCirtKdHet_impl, X, Lambda, Lambda_eq, Lambda_ineq, Lambda_prior_mean,
			   Lambda_prior_prec, *sigmapriorc, *sigmapriord, *storelambda, 
			   *storescores, *storesigma,
			   *burnin, *mcmc, *thin, *verbose, output);


    for (unsigned int i = 0; i < output.size(); ++i)
      samdata[i] = output(i);

  }
  
}

#endif
