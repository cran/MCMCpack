//////////////////////////////////////////////////////////////////////////
// MCMCpaircompare.cc is C++ code to estimate a pairwise comparison model. 
//
// KQ 3/18/2015
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////

#ifndef MCMCPAIRCOMPARE_CC
#define MCMCPAIRCOMPARE_CC

#include<vector>

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

/* MCMCpaircompare implementation. */
template <typename RNGTYPE>
void MCMCpaircompare_impl (rng<RNGTYPE>& stream,
			   const Matrix<unsigned int>& MD, 
			   Matrix<>& theta,
			   Matrix<>& alpha,
			   const Matrix<>& theta_eq, 
			   const Matrix<>& theta_ineq,
			   const double a0,
			   const double A0, 
			   const int alphafixed,
			   const unsigned int burnin,
			   const unsigned int mcmc, 
			   const unsigned int thin,
			   const unsigned int verbose,
			   const bool storealpha,
			   const bool storetheta, 
			   double* sampledata,
			   const unsigned int samplesize){


  // constants
  const unsigned int N = MD.rows();
  const unsigned int J = theta.rows();
  const unsigned int I = alpha.rows();
  const unsigned int tot_iter = burnin + mcmc;  
  const unsigned int nsamp = mcmc / thin;

  
  // starting values for Ystar
  Matrix<> Ystar(N,1);
  

  // pre-compute what can be pre-computed
  const double A0a0 = A0 * a0;
  Matrix<unsigned int> theta_n(J,1,true,0); // vector of 0s. Item j's total instances of being compared.
  Matrix<unsigned int> alpha_n(I,1,true,0); // vector of 0s. Judge i's total instances of making comparisons.
  for (unsigned int i = 0; i < N; ++i){
    alpha_n(MD(i,0)) += 1;
    theta_n(MD(i,1)) += 1;
    theta_n(MD(i,2)) += 1;
  }

  vector< vector< double* > > theta_Ystar_ptr;
  vector< vector< double* > > theta_alpha_ptr;
  vector< vector< double* > > theta_theta_ptr;
  vector< vector< double > > theta_sign;

  vector< vector< double* > > alpha_Ystar_ptr;
  vector< vector< double* > > alpha_theta1_ptr;
  vector< vector< double* > > alpha_theta2_ptr;

  theta_Ystar_ptr.reserve(J);
  theta_alpha_ptr.reserve(J);
  theta_theta_ptr.reserve(J);
  theta_sign.reserve(J);

  alpha_Ystar_ptr.reserve(I);
  alpha_theta1_ptr.reserve(I);
  alpha_theta2_ptr.reserve(I);

  for (unsigned int j = 0; j < J; ++j){
    vector< double* > Ystar_j_ptr;
    vector< double* > alpha_j_ptr;
    vector< double* > theta_j_ptr;
    vector< double > sign_j;
	 
    Ystar_j_ptr.reserve(theta_n(j));
    alpha_j_ptr.reserve(theta_n(j));
    theta_j_ptr.reserve(theta_n(j));
    sign_j.reserve(theta_n(j));
	/*
	Allocate proper length for the empty vectors (Ystar_j_ptr, alpha_j_ptr, theta_j_ptr, sign_j)
	Then we save these pre-allocated vectors inside the outer layer vectors.
    */
	theta_Ystar_ptr.push_back(Ystar_j_ptr);
    theta_alpha_ptr.push_back(alpha_j_ptr);
    theta_theta_ptr.push_back(theta_j_ptr);
    theta_sign.push_back(sign_j);
  }

  for (unsigned int i = 0; i < I; ++i){
    vector< double* > Ystar_i_ptr;
    vector< double* > theta1_i_ptr;
    vector< double* > theta2_i_ptr;

    Ystar_i_ptr.reserve(alpha_n(i));
    theta1_i_ptr.reserve(alpha_n(i));
    theta2_i_ptr.reserve(alpha_n(i));
	/*
	Allocate proper length for the empty vectors (Ystar_j_ptr, theta1_i_ptr, theta2_i_ptr)
	Then we save these pre-allocated vectors inside the outer layer vectors.
	*/
    alpha_Ystar_ptr.push_back(Ystar_i_ptr);
    alpha_theta1_ptr.push_back(theta1_i_ptr);
    alpha_theta2_ptr.push_back(theta2_i_ptr);
  }

  
  for (unsigned int i = 0; i < N; ++i){
    unsigned int resp = MD(i,0); 
    unsigned int c1 = MD(i,1);
    unsigned int c2 = MD(i,2);
	/*
	Assign the starting values of Ystar, alpha, theta and sign to the above vectors.
	This is convenient for accessving wanted objects in later updation.
	*/
    theta_Ystar_ptr[c1].push_back(&Ystar(i,0));
    theta_alpha_ptr[c1].push_back(&alpha(MD(i,0)));
    theta_theta_ptr[c1].push_back(&theta(MD(i,2)));
    theta_sign[c1].push_back(1.0);
    
    theta_Ystar_ptr[c2].push_back(&Ystar(i,0));
    theta_alpha_ptr[c2].push_back(&alpha(MD(i,0)));
    theta_theta_ptr[c2].push_back(&theta(MD(i,1)));
    theta_sign[c2].push_back(-1.0);

    alpha_Ystar_ptr[resp].push_back(&Ystar(i,0));
    alpha_theta1_ptr[resp].push_back(&theta(MD(i,1)));
    alpha_theta2_ptr[resp].push_back(&theta(MD(i,2)));
  }


  
  
  
  // storage matrices (col major order)
  Matrix<> theta_store;
  Matrix<> alpha_store;
  if (storetheta)
    theta_store = Matrix<>(nsamp, J);
  
  if (storealpha)
    alpha_store = Matrix<>(nsamp, I);



  unsigned int count = 0;
  // MCMC sampling occurs in this for loop
  for (unsigned int iter = 0; iter < tot_iter; ++iter){

	// The following three functions are defined in MCMCfcds.h
    // sample Ystar
    paircompare_Ystar_update(Ystar, MD, theta, alpha, stream);

    // sample theta
    paircompare_theta_update(theta, Ystar, MD, alpha, theta_n, theta_eq,
			     theta_ineq,
			     theta_Ystar_ptr,
			     theta_alpha_ptr,
			     theta_theta_ptr,
			     theta_sign,
			     stream);

    if (alphafixed == 0){
    // sample alpha
    paircompare_alpha_update(alpha, Ystar, MD, theta, A0, A0a0,
			     alpha_n, 
			     alpha_Ystar_ptr,
			     alpha_theta1_ptr,
			     alpha_theta2_ptr,
			     stream);
    }
    
    // print results to screen
    if (verbose > 0 && iter % verbose == 0) {
      Rprintf("\n\nMCMCpaircompare iteration %i of %i \n", (iter+1), tot_iter);
      //Rprintf("theta = \n");
      //for (int j=0; j<J; ++j)
      //Rprintf("%10.5f\n", theta[j]);    
    }
    
    // store results
    if ((iter >= burnin) && ((iter % thin == 0))) {
      
      // store theta
      if (storetheta)
        theta_store(count, _) = theta;

      // store alpha
      if (storealpha)
        alpha_store(count, _) = alpha;

      count++;	
    }
    
    R_CheckUserInterrupt(); // allow user interrupts
       
  } // end MCMC sampling loop
 


  
  // put output back into sampledata
  Matrix<> output;
  if(storetheta && ! storealpha) {
    output = theta_store;
  } else if (storealpha && ! storetheta){
    output = alpha_store;
  } else {
    output = cbind(theta_store, alpha_store);
  }

  for (unsigned int i = 0; i < samplesize; ++i)
    sampledata[i] = output[i];


} // end MCMCpaircompare implementation









extern "C" {

  void
  cMCMCpaircompare(double* sampledata, 
		  const int* samplerow, 
		  const int* samplecol,
		  const unsigned int* MDdata, 
		  const int* MDrow, 
		  const int* MDcol,
		  const int* alphafixed,
		  const int* burnin, 
		  const int* mcmc,  
		  const int* thin,
		  const int *uselecuyer, 
		  const int *seedarray, 
		  const int *lecuyerstream, 
		  const int* verbose, 
		  const double* thetastartdata,
		  const int* thetastartrow, 
		  const int* thetastartcol, 
		  const double* astartdata, 
		  const int* astartrow, 
		  const int* astartcol,
		  const double* a0, 
		  const double* A0,	
		  const double* thetaeqdata, 
		  const int* thetaeqrow, 
		  const int* thetaeqcol,
		  const double* thetaineqdata, 
		  const int* thetaineqrow, 
		  const int* thetaineqcol, 
		  const int* storealpha, 
		  const int* storetheta){


    // put together matrices
    const Matrix<unsigned int> MD(*MDrow, *MDcol, MDdata);
    Matrix<> theta(*thetastartrow, *thetastartcol, thetastartdata);
    Matrix<> alpha(*astartrow, *astartcol, astartdata);
    const Matrix<> theta_eq(*thetaeqrow, *thetaeqcol, thetaeqdata);
    const Matrix<> theta_ineq(*thetaineqrow, *thetaineqcol, thetaineqdata);
    const int samplesize = (*samplerow) * (*samplecol);

    
    MCMCPACK_PASSRNG2MODEL(MCMCpaircompare_impl, MD, theta, alpha, 
			   theta_eq, theta_ineq, *a0, *A0,
			   *alphafixed,
			   *burnin, *mcmc, *thin, 
			   *verbose, *storealpha, *storetheta, 
			   sampledata, samplesize);
  }
}


#endif










