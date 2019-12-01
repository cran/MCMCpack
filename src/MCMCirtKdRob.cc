//////////////////////////////////////////////////////////////////////////
// MCMCirtKdRob.cc is C++ code to estimate a robust K-dimensional
// item response model 
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
// 2/18/2005 KQ
// 8/1/2007  ported to Scythe 1.0.2 KQ
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////



#ifndef MCMCIRTKDROB_CC
#define MCMCIRTKDROB_CC

#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"
#include "MCMCmnl.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>



typedef Matrix<double,Row,View> rmview;

using namespace std;
using namespace scythe;


/* Equal probability sampling; without-replacement case */
// pulled from R-2.0.1/src/main/random.c lines 352-364
// slightly modified by KQ 2/21/2005
// x: n array of original indices running from 0 to (n-1)
// y: k array of samled indices
// k: length of y (must be <= n)
// n: length of x (must be >= k)
template <typename RNGTYPE>
static void SampleNoReplace(const int k, int n, int *y, int *x, 
			    rng<RNGTYPE>& stream){
  
  for (int i = 0; i < n; i++)
    x[i] = i;
  for (int i = 0; i < k; i++) {
    int j = static_cast<int>(n * stream.runif());
    y[i] = x[j];
    x[j] = x[--n];
  }
}



// full conditional distribution for a single element of Lambda
// this single element is Lambda(rowindex, colindex)
static double Lambda_logfcd(const double& lam_ij,
			    const Matrix<int>& X,
			    const Matrix<>& Lambda,
			    const Matrix<>& theta, 
			    const double& delta0,
			    const double& delta1, 
			    const Matrix<>& Lambda_prior_mean, 
			    const Matrix<>& Lambda_prior_prec,
			    const Matrix<>& Lambda_ineq,
			    const Matrix<>& theta_ineq,
			    const double& k0,
			    const double& k1,
			    const double& c0,
			    const double& d0,
			    const double& c1,
			    const double& d1,
			    const int& rowindex,
			    const int& colindex){


  const int D = Lambda.cols();
  
  // check to see if inequality constraint is satisfied and 
  // evaluate prior
    
  double logprior = 0.0;  
  if (Lambda_ineq(rowindex,colindex) * lam_ij < 0){
    return log(0.0);
  } 
  if (Lambda_prior_prec(rowindex,colindex) != 0){
    logprior += lndnorm(lam_ij, 
			Lambda_prior_mean(rowindex,colindex), 
			sqrt(1.0 / Lambda_prior_prec(rowindex,colindex)));
  } 
  
  // prior is uniform on hypersphere with radius 10
  /*
    if (Lambda_ineq(rowindex,colindex) * lam_ij < 0){
    return log(0.0);
    } 
    double absqdist = 0.0;
    for (int i=0; i<D; ++i){
    if (i==colindex){
    absqdist += ::pow(lam_ij, 2);
    }
    else{
    absqdist += ::pow(Lambda(rowindex,i), 2);
    }
    }
    if (absqdist > 100.0){
    return log(0.0);
    }
    const double logprior = 0.0;
  */

  // likelihood
  double loglike = 0.0;
  const int N = X.rows();
  for (int i=0; i<N; ++i){
    if (X(i,rowindex) != -999){
      double eta = 0.0;
      for (int j=0; j<D; ++j){
	if (j==colindex){
	  eta += theta(i,j) * lam_ij;
	}
	else{
	  eta += theta(i,j) * Lambda(rowindex,j);
	}
      } 
      const double p = delta0 + (1 - delta0 - delta1) * 
	(1.0 / (1.0 + exp(-1*eta))); 
      loglike += X(i,rowindex) * log(p) + (1.0 - X(i,rowindex)) * log(1.0 - p);
    } 
  }
  
  return (loglike + logprior);
}



// full conditional for a single element of theta
// this single element is theta(rowindex, colindex)
static double theta_logfcd(const double& t_ij,
			   const Matrix<int>& X, 
			   const Matrix<>& Lambda,
			   const Matrix<>& theta, 
			   const double& delta0,
			   const double& delta1, 
			   const Matrix<>& Lambda_prior_mean, 
			   const Matrix<>& Lambda_prior_prec,
			   const Matrix<>& Lambda_ineq,
			   const Matrix<>& theta_ineq,
			   const double& k0,
			   const double& k1,
			   const double& c0,
			   const double& d0,
			   const double& c1,
			   const double& d1,
			   const int& rowindex,
			   const int& colindex){
 
  const int D = Lambda.cols();     
  // evaluate prior  
  if (theta_ineq(rowindex,colindex-1) * t_ij < 0){
    return log(0.0);
  } 
  const double logprior = lndnorm(t_ij, 0.0, 1.0); 
  
  // prior is uniform on unit circle
  /*
    if (theta_ineq(rowindex,colindex-1) * t_ij < 0){
    return log(0.0);
    } 
    double tsqdist = 0.0;
    for (int i=0; i<(D-1); ++i){
    if (i==(colindex-1)){
    tsqdist += ::pow(t_ij, 2);
    }
    else{
    tsqdist += ::pow(theta(rowindex,(i-1)), 2);
    }
    }
    if (tsqdist > 1.0){
    return log(0.0);
    }
    const double logprior = 1.0;
  */

  // likelihood
  double loglike = 0.0;
  const int K = X.cols();
  for (int i=0; i<K; ++i){
    if (X(rowindex,i) != -999){
      double eta = 0.0;
      for (int j=0; j<D; ++j){
	if (j==colindex){
	  eta += t_ij * Lambda(i,j);
	}
	else{
	  eta += theta(rowindex,j) * Lambda(i,j);
	}
      }       
      const double p = delta0 + (1 - delta0 - delta1) * 
	(1.0 / (1.0 + exp(-1*eta)));       
      loglike += X(rowindex,i) * log(p) + (1.0 - X(rowindex,i)) * log(1.0 - p);
    } 
  }
  
  return (loglike + logprior);
}



// full conditional for delta0
static double delta0_logfcd(const double& delta0,
			    const Matrix<int>& X, 
			    const Matrix<>& Lambda,
			    const Matrix<>& theta, 
			    const double& junk,
			    const double& delta1, 
			    const Matrix<>& Lambda_prior_mean, 
			    const Matrix<>& Lambda_prior_prec,
			    const Matrix<>& Lambda_ineq,
			    const Matrix<>& theta_ineq,
			    const double& k0,
			    const double& k1,
			    const double& c0,
			    const double& d0,
			    const double& c1,
			    const double& d1,
			    const int& rowindex,
			    const int& colindex){
  
  // evaluate prior
  if (delta0 >=k0 || delta0 <=0){
    return log(0.0);
  }
  const double logprior = lndbeta1(delta0 * (1.0/k0), c0, d0); 
  
  // likelihood
  double loglike = 0.0;
  const int D = Lambda.cols();
  const int N = X.rows();
  const int K = X.cols();
  for (int i=0; i<N; ++i){
    for (int k=0; k<K; ++k){
      if (X(i,k) != -999){
	double eta = 0.0;
	for (int j=0; j<D; ++j){
	  eta += theta(i,j) * Lambda(k,j);
	} 
	const double p = delta0 + (1 - delta0 - delta1) * 
	  (1.0 / (1.0 + exp(-1*eta))); 
	loglike += X(i,k) * log(p) + (1.0 - X(i,k)) * log(1.0 - p);
      }
    } 
  }
  
  return (loglike + logprior);
}



// full conditional for delta1
static double delta1_logfcd(const double& delta1, 
			    const Matrix<int>& X, 
			    const Matrix<>& Lambda,
			    const Matrix<>& theta, 
			    const double& delta0,
			    const double& junk, 
			    const Matrix<>& Lambda_prior_mean, 
			    const Matrix<>& Lambda_prior_prec,
			    const Matrix<>& Lambda_ineq,
			    const Matrix<>& theta_ineq,
			    const double& k0,
			    const double& k1,
			    const double& c0,
			    const double& d0,
			    const double& c1,
			    const double& d1,
			    const int& rowindex,
			    const int& colindex){
  
  // evaluate prior
  if (delta1 >=k1 || delta1 <=0){
    return log(0.0);
  }
  const double logprior = lndbeta1(delta1 * (1.0/k1), c1, d1); 
  
  // likelihood
  double loglike = 0.0;
  const int D = Lambda.cols();
  const int N = X.rows();
  const int K = X.cols();
  for (int i=0; i<N; ++i){
    for (int k=0; k<K; ++k){
      if (X(i,k) != -999){
	double eta = 0.0;
	for (int j=0; j<D; ++j){
	  eta += theta(i,j) * Lambda(k,j);
	} 
	const double p = delta0 + (1 - delta0 - delta1) * 
	  (1.0 / (1.0 + exp(-1*eta))); 
	loglike += X(i,k) * log(p) + (1.0 - X(i,k)) * log(1.0 - p);
      }
    } 
  }
  
  return (loglike + logprior);
}





// Radford Neal's (2000) doubling procedure coded for a logdensity
template<typename RNGTYPE>
static void doubling(double (*logfun)(const double&,
				      const Matrix<int>&,
				      const Matrix<>&,
				      const Matrix<>&,
				      const double&, 
				      const double&,
				      const Matrix<>&,
				      const Matrix<>&,
				      const Matrix<>&,
				      const Matrix<>&,
				      const double&, 
				      const double&,
				      const double&, 
				      const double&,
				      const double&, 
				      const double&,
				      const int&,
				      const int&), 
		     const Matrix<int>& X,
		     const Matrix<>& Lambda,
		     const Matrix<>& theta,
		     const double& delta0,
		     const double& delta1, 
		     const Matrix<>& Lambda_prior_mean, 
		     const Matrix<>& Lambda_prior_prec,
		     const Matrix<>& Lambda_ineq,
		     const Matrix<>& theta_ineq,
		     const double& k0,
		     const double& k1,
		     const double& c0,
		     const double& d0,
		     const double& c1, 
		     const double& d1,
		     const int& rowindex,
		     const int& colindex,
		     const double& z, 
		     const double& w, const int& p, 
		     rng<RNGTYPE>& stream, double& L, double& R, 
		     const int& param){
  
  const double U = stream.runif();
  double x0 = 0.0;
  if (param==0){ // Lambda
    x0 = Lambda(rowindex, colindex);
  }
  else if (param==1){ // theta
    x0 = theta(rowindex, colindex);
  }
  else if (param==2){ // delta0
    x0 = delta0;
  }
  else if (param==3){ // delta1
    x0 = delta1;
  }
  else {
    error("ERROR: param not in {0,1,2,3} in doubling().");
    //Rprintf("\nERROR: param not in {0,1,2,3} in doubling().\n");
    //exit(1);    
  }

  L = x0 - w*U;
  R = L + w;
  int K = p;
  while (K > 0 && 
	 (z < logfun(L, X, Lambda, theta, delta0, delta1, Lambda_prior_mean,
		     Lambda_prior_prec, Lambda_ineq, theta_ineq, 
		     k0, k1, c0, d0, c1, d1, rowindex, colindex) ||
	  z < logfun(R, X, Lambda, theta, delta0, delta1, Lambda_prior_mean,
		     Lambda_prior_prec, Lambda_ineq, theta_ineq, 
		     k0, k1, c0, d0, c1, d1, rowindex, colindex))){
    double V = stream.runif();
    if (V < 0.5){
      L = L - (R - L);
    }
    else {
      R = R + (R - L);
    }
    --K;
  }      
  
}


// Radford Neal's (2000) stepping out procedure coded for a logdensity
template<typename RNGTYPE>
static void StepOut(double (*logfun)(const double&,
				     const Matrix<int>&,
				     const Matrix<>&,
				     const Matrix<>&,
				     const double&, 
				     const double&,
				     const Matrix<>&,
				     const Matrix<>&,
				     const Matrix<>&,
				     const Matrix<>&,
				     const double&, 
				     const double&,
				     const double&, 
				     const double&,
				     const double&, 
				     const double&,
				     const int&,
				     const int&), 
		    const Matrix<int>& X,
		    const Matrix<>& Lambda,
		    const Matrix<>& theta,
		    const double& delta0,
		    const double& delta1, 
		    const Matrix<>& Lambda_prior_mean, 
		    const Matrix<>& Lambda_prior_prec,
		    const Matrix<>& Lambda_ineq,
		    const Matrix<>& theta_ineq,
		    const double& k0,
		    const double& k1,
		    const double& c0,
		    const double& d0,
		    const double& c1, 
		    const double& d1,
		    const int& rowindex,
		    const int& colindex,
		    const double& z, 
		    const double& w, const int& m, 
		    rng<RNGTYPE>& stream, double& L, double& R, 
		    const int& param){
  
  const double U = stream.runif();
  double x0 = 0.0;
  if (param==0){ // Lambda
    x0 = Lambda(rowindex, colindex);
  }
  else if (param==1){ // theta
    x0 = theta(rowindex, colindex);
  }
  else if (param==2){ // delta0
    x0 = delta0;
  }
  else if (param==3){ // delta1
    x0 = delta1;
  }
  else {
    error("ERROR: param not in {0,1,2,3} in StepOut().");
    //Rprintf("\nERROR: param not in {0,1,2,3} in StepOut().\n");
    //exit(1);    
  }
  

  L = x0 - w*U;
  R = L + w;
  const double V = stream.runif();
  int J = static_cast<int>(m*V);
  int K = (m-1) - J; 
  while (J > 0 && 
	 (z < logfun(L, X, Lambda, theta, delta0, delta1, Lambda_prior_mean,
		     Lambda_prior_prec, Lambda_ineq, theta_ineq, 
		     k0, k1, c0, d0, c1, d1, rowindex, colindex))){
    L = L - w;
    J = J - 1;	   
  }
  while (K > 0 && 
	 (z < logfun(R, X, Lambda, theta, delta0, delta1, Lambda_prior_mean,
		     Lambda_prior_prec, Lambda_ineq, theta_ineq, 
		     k0, k1, c0, d0, c1, d1, rowindex, colindex))){
    R = R + w;
    K = K - 1;	   
  }
  
}




// Radford Neal's (2000) Accept procedure coded for a logdensity
static bool Accept(double (*logfun)(const double&,
					  const Matrix<int>&,
					  const Matrix<>&,
					  const Matrix<>&,
					  const double&, 
					  const double&,
					  const Matrix<>&,
					  const Matrix<>&,
					  const Matrix<>&,
					  const Matrix<>&,
					  const double&, 
					  const double&,
					  const double&, 
					  const double&,
					  const double&, 
					  const double&,
					  const int&,
					  const int&), 
			 const Matrix<int>& X,
			 const Matrix<>& Lambda,
			 const Matrix<>& theta,
			 const double& delta0,
			 const double& delta1, 
			 const Matrix<>& Lambda_prior_mean, 
			 const Matrix<>& Lambda_prior_prec,
			 const Matrix<>& Lambda_ineq,
			 const Matrix<>& theta_ineq,
			 const double& k0,
			 const double& k1,
			 const double& c0,
			 const double& d0,
			 const double& c1, 
			 const double& d1, 
			 const int& rowindex,
			 const int& colindex,
			 const double& z,
			 const double& w, const double& x0, 
			 const double& x1, const double& L, const double& R){
  
  double Lhat = L;
  double Rhat = R;
  bool D = false;

  while ((Rhat - Lhat ) > 1.1 * w){
    double M = (Lhat + Rhat) / 2.0;
    if ( (x0 < M && x1 >= M) || (x0 >= M && x1 < M)){
      D = true;
    }
    if (x1 < M){
      Rhat = M;
    }
    else {
      Lhat = M;
    }
    
    if (D && z >= logfun(Lhat, X, Lambda, theta, delta0, delta1, 
			 Lambda_prior_mean, Lambda_prior_prec, 
			 Lambda_ineq, theta_ineq, k0, k1, c0, d0, 
			 c1, d1, rowindex, colindex) && 
	z >=  logfun(Rhat, X, Lambda, theta, delta0, delta1, Lambda_prior_mean,
		     Lambda_prior_prec, Lambda_ineq, theta_ineq, 
		     k0, k1, c0, d0, c1, d1, rowindex, colindex)){
      return(false);
    }    
  }
  return(true);
}


// Radford Neal's (2000) shrinkage procedure coded for a log density
// assumes the doubling procedure has been used to find L and R
template <typename RNGTYPE>
static double shrinkageDoubling(double (*logfun)(const double&,
						 const Matrix<int>&,
						 const Matrix<>&,
						 const Matrix<>&,
						 const double&, 
						 const double&,
						 const Matrix<>&,
						 const Matrix<>&,
						 const Matrix<>&,
						 const Matrix<>&,
						 const double&, 
						 const double&,
						 const double&, 
						 const double&,
						 const double&, 
						 const double&,
						 const int&,
						 const int&), 
				const Matrix<int>& X,
				const Matrix<>& Lambda,
				const Matrix<>& theta,
				const double& delta0,
				const double& delta1, 
				const Matrix<>& Lambda_prior_mean, 
				const Matrix<>& Lambda_prior_prec,
				const Matrix<>& Lambda_ineq,
				const Matrix<>& theta_ineq,
				const double& k0,
				const double& k1,
				const double& c0,
				const double& d0,
				const double& c1, 
				const double& d1,
				const int& rowindex,
				const int& colindex,
				const double& z, const double& w, 
				rng<RNGTYPE>& stream, const double& L, 
				const double& R,
				const int& param){
  
  double Lbar = L;
  double Rbar = R;
  double x0;
  if (param==0){ // Lambda
    x0 = Lambda(rowindex, colindex);
  }
  else if (param==1){ // theta
    x0 = theta(rowindex, colindex);
  }
  else if (param==2){ // delta0
    x0 = delta0;
  }
  else if (param==3){ // delta1
    x0 = delta1;
  }
  else {
    error("ERROR: param not in {0,1,2,3} in shrinkageDoubling().");
    //Rprintf("\nERROR: param not in {0,1,2,3} in shrinkageDoubling().\n");
    //exit(1);    
  }

  for (;;){
    const double U = stream.runif();
    const double x1 = Lbar + U*(Rbar - Lbar);
    if (z <= logfun(x1, X, Lambda, theta, delta0, delta1, Lambda_prior_mean,
		    Lambda_prior_prec, Lambda_ineq, theta_ineq, 
		    k0, k1, c0, d0, c1, d1, rowindex, colindex) &&
	Accept(logfun, X, Lambda, theta, delta0, delta1, Lambda_prior_mean, 
	       Lambda_prior_prec, Lambda_ineq, theta_ineq, 
	       k0, k1, c0, d0, c1, d1, rowindex, colindex, z, w, 
	       x0, x1, L, R)){
      return(x1);
    }
    if (x1 < x0){
      Lbar = x1;
    }
    else {
      Rbar = x1;
    }
  } // end infinite loop
}



// Radford Neal's (2000) shrinkage procedure coded for a log density
// assumes the stepping out procedure has been used to find L and R
template <typename RNGTYPE>
static double shrinkageStep(double (*logfun)(const double&,
					     const Matrix<int>&,
					     const Matrix<>&,
					     const Matrix<>&,
					     const double&, 
					     const double&,
					     const Matrix<>&,
					     const Matrix<>&,
					     const Matrix<>&,
					     const Matrix<>&,
					     const double&, 
					     const double&,
					     const double&, 
					     const double&,
					     const double&, 
					     const double&,
					     const int&,
					     const int&), 
			    const Matrix<int>& X,
			    const Matrix<>& Lambda,
			    const Matrix<>& theta,
			    const double& delta0,
			    const double& delta1, 
			    const Matrix<>& Lambda_prior_mean, 
			    const Matrix<>& Lambda_prior_prec,
			    const Matrix<>& Lambda_ineq,
			    const Matrix<>& theta_ineq,
			    const double& k0,
			    const double& k1,
			    const double& c0,
			    const double& d0,
			    const double& c1, 
			    const double& d1,
			    const int& rowindex,
			    const int& colindex,
			    const double& z, const double& w, 
			    rng<RNGTYPE>& stream, const double& L, 
			    const double& R,
			    const int& param){
  
  double Lbar = L;
  double Rbar = R;
  double x0;
  if (param==0){ // Lambda
    x0 = Lambda(rowindex, colindex);
  }
  else if (param==1){ // theta
    x0 = theta(rowindex, colindex);
  }
  else if (param==2){ // delta0
    x0 = delta0;
  }
  else if (param==3){ // delta1
    x0 = delta1;
  }
  else {
    error("ERROR: param not in {0,1,2,3} in shrinkageDoubling().");
    //Rprintf("\nERROR: param not in {0,1,2,3} in shrinkageDoubling().\n");
    //exit(1);    
  }

  for (;;){
    const double U = stream.runif();
    const double x1 = Lbar + U*(Rbar - Lbar);
    if (z <= logfun(x1, X, Lambda, theta, delta0, delta1, Lambda_prior_mean,
		    Lambda_prior_prec, Lambda_ineq, theta_ineq, 
		    k0, k1, c0, d0, c1, d1, rowindex, colindex) ){
      return(x1);
    }
    if (x1 < x0){
      Lbar = x1;
    }
    else {
      Rbar = x1;
    }
  } // end infinite loop
}













template <typename RNGTYPE>
void MCMCirtKdRob_impl(rng<RNGTYPE>& stream,
		       const Matrix<int>& X, 
		       Matrix<>& Lambda,
		       Matrix<>& theta,
		       const Matrix<>& Lambda_eq, const Matrix<>& Lambda_ineq,
		       const Matrix<>& theta_eq, const Matrix<>& theta_ineq,
		       const Matrix<>& Lambda_prior_mean,
		       const Matrix<>& Lambda_prior_prec,
		       const int* burnin, const int* mcmc,  const int* thin,
		       const int* verbose, 
		       const int* method_step,
		       const double* theta_w, const int* theta_p,
		       const double* lambda_w, const int* lambda_p,
		       const double* delta0_w, const int* delta0_p,
		       const double* delta1_w, const int* delta1_p, 
		       const double * delta0start, const double* delta1start,
		       const double* k0, const double* k1,
		       const double* c0, const double* c1,
		       const double* d0, const double* d1,
		       const int* storeitem,
		       const int* storeability,
		       double* sampledata, const int* samplerow, 
		       const int* samplecol
		       ){



    // constants 
    const int K = X.cols();  // number of items
    const int N = X.rows();  // number of subjects
    const int D = Lambda.cols();  // number of dimensions + 1
    const int tot_iter = *burnin + *mcmc;  
    const int nsamp = *mcmc / *thin;
    // const Matrix<double> Lambda_free_indic = Matrix<double>(K, D);
    //for (int i=0; i<(K*D); ++i){
    //  if (Lambda_eq[i] == -999) Lambda_free_indic[i] = 1.0;
    //}

    // starting values
    //  Matrix<double> theta = Matrix<double>(N, D-1);
    Matrix<> onesvec = ones<double>(N, 1);
    onesvec = onesvec * -1.0;
    theta = cbind(onesvec, theta);
    double delta0 = *delta0start;
    double delta1 = *delta1start;

    // index arrays (used for the random order of the sampling)
    // OLD
    //int K_array[K];
    //int N_array[N];
    //int D_array[D];
    //int Dm1_array[D-1];
    //int K_inds_perm[K];
    //int N_inds_perm[N];
    //int D_inds_perm[D];
    //int Dm1_inds_perm[D-1];
    // NEW
    int* K_array = new int[K];
    int* N_array = new int[N];
    int* D_array = new int[D];
    int* Dm1_array = new int[D-1];
    int* K_inds_perm = new int[K];
    int* N_inds_perm = new int[N];
    int* D_inds_perm = new int[D];
    int* Dm1_inds_perm = new int[D-1];


    // storage matrices (row major order)
    Matrix<> Lambda_store;
    if (storeitem[0]==1){
      Lambda_store = Matrix<double>(nsamp, K*D);
    }
    Matrix<> theta_store;
    if (*storeability==1){
      theta_store = Matrix<double>(nsamp, N*D);
    }
    Matrix<> delta0_store(nsamp, 1);
    Matrix<> delta1_store(nsamp, 1);

    ///////////////////
    // Slice Sampler //
    ///////////////////
    int count = 0;  
    for (int iter=0; iter < tot_iter; ++iter){

      double L, R, w, funval, z;
      int p;

      // sample theta
      int param = 1;    
      SampleNoReplace(N, N, N_inds_perm, N_array, stream);   
      for (int ii=0; ii<N; ++ii){
	int i = N_inds_perm[ii];
	SampleNoReplace(D-1, D-1, Dm1_inds_perm, Dm1_array, stream);    
	for (int jj=0; jj<(D-1); ++jj){
	  int j = Dm1_inds_perm[jj]+1;
	  if (theta_eq(i,j-1) == -999){
	    w = *theta_w;
	    p = *theta_p;
	    L = -1.0;
	    R = 1.0;
	    funval = theta_logfcd(theta(i,j), X, Lambda, 
				  theta, delta0,
				  delta1, Lambda_prior_mean, 
				  Lambda_prior_prec,
				  Lambda_ineq, theta_ineq, 
				  *k0, *k1, *c0, *d0,
				  *c1, *d1, i, j);
	    z = funval - stream.rexp(1.0);

	    if (*method_step == 1){
	      StepOut(&theta_logfcd, X, Lambda, theta, delta0, delta1, 
		      Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		      theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, i, j, z, 
		      w, p, stream, L, R, param);
	    
	      theta(i,j) = shrinkageStep(&theta_logfcd, X, Lambda, theta,
					 delta0, delta1, Lambda_prior_mean, 
					 Lambda_prior_prec, Lambda_ineq,
					 theta_ineq, *k0, *k1, 
					 *c0, *d0, *c1, 
					 *d1, i, j, z, w, stream, 
					 L, R, param);
	    }
	    else{
	      doubling(&theta_logfcd, X, Lambda, theta, delta0, delta1, 
		       Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		       theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, i, j, z, 
		       w, p, stream, L, R, param);
	      theta(i,j) = shrinkageDoubling(&theta_logfcd, X, Lambda, theta,
					     delta0, delta1, Lambda_prior_mean, 
					     Lambda_prior_prec, Lambda_ineq,
					     theta_ineq, *k0, *k1, 
					     *c0, *d0, *c1, 
					     *d1, i, j, z, w, stream, 
					     L, R, param);
	    }
	  }
	}
      }
    
    

      // sample Lambda
      param = 0;
      SampleNoReplace(K, K, K_inds_perm, K_array, stream);   
      for (int ii=0; ii<K; ++ii){
	int i = K_inds_perm[ii];
	SampleNoReplace(D, D, D_inds_perm, D_array, stream);   
	for (int jj=0; jj<D; ++jj){
	  int j = D_inds_perm[jj];
	  if (Lambda_eq(i,j) == -999){
	    w = *lambda_w;
	    p = *lambda_p;
	    L = -1.0;
	    R = 1.0;
	    funval = Lambda_logfcd(Lambda(i,j), X, Lambda, 
				   theta, delta0,
				   delta1, Lambda_prior_mean, 
				   Lambda_prior_prec,
				   Lambda_ineq, theta_ineq, *k0, *k1, *c0, *d0,
				   *c1, *d1, i, j);
	    z = funval - stream.rexp(1.0);
	    if (*method_step == 1){
	      StepOut(&Lambda_logfcd, X, Lambda, theta, delta0, delta1, 
		      Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		      theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, i, j, z, 
		      w, p, stream, L, R, param);
	      Lambda(i,j) = shrinkageStep(&Lambda_logfcd, X, Lambda, theta,
					  delta0, delta1, Lambda_prior_mean, 
					  Lambda_prior_prec, Lambda_ineq,
					  theta_ineq, *k0, *k1, *c0, *d0, 
					  *c1, *d1, i, j, z, w, stream, 
					  L, R, param);
	    }
	    else{
	      doubling(&Lambda_logfcd, X, Lambda, theta, delta0, delta1, 
		       Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		       theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, i, j, z, 
		       w, p, stream, L, R, param);
	      Lambda(i,j) = shrinkageDoubling(&Lambda_logfcd, X, Lambda, theta,
					      delta0, delta1, Lambda_prior_mean, 
					      Lambda_prior_prec, Lambda_ineq,
					      theta_ineq, *k0, *k1, *c0, *d0, 
					      *c1, *d1, i, j, z, w, stream, 
					      L, R, param);
	    }
	  }
	}
      }
    

    
      // sample delta0
      param = 2;
      w = *delta0_w;
      p = *delta0_p;
      L = -1.0;
      R =  1.0;    
      funval = delta0_logfcd(delta0, X, Lambda, 
			     theta, delta0,
			     delta1, Lambda_prior_mean, 
			     Lambda_prior_prec,
			     Lambda_ineq, theta_ineq, 
			     *k0, *k1, *c0, *d0,
			     *c1, *d1, 0, 0);
      z = funval - stream.rexp(1.0);
      if (*method_step == 1){
	StepOut(&delta0_logfcd, X, Lambda, theta, delta0, delta1, 
		Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, 0, 0, z, 
		w, p, stream, L, R, param);
	delta0 = shrinkageStep(&delta0_logfcd, X, Lambda, theta,
			       delta0, delta1, Lambda_prior_mean, 
			       Lambda_prior_prec, Lambda_ineq, theta_ineq,
			       *k0, *k1, *c0, *d0, *c1, *d1, 0, 0,
			       z, w, stream, L, R, param);    
      }
      else{
	doubling(&delta0_logfcd, X, Lambda, theta, delta0, delta1, 
		 Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		 theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, 0, 0, z, 
		 w, p, stream, L, R, param);
	delta0 = shrinkageDoubling(&delta0_logfcd, X, Lambda, theta,
				   delta0, delta1, Lambda_prior_mean, 
				   Lambda_prior_prec, Lambda_ineq, theta_ineq,
				   *k0, *k1, *c0, *d0, *c1, *d1, 0, 0,
				   z, w, stream, L, R, param);
      }



      // sample delta1
      param = 3;
      w = *delta1_w;
      p = *delta1_p;
      L = -1.0;
      R = 1.0; 
      funval = delta1_logfcd(delta1, X, Lambda, 
			     theta, delta0,
			     delta1, Lambda_prior_mean, 
			     Lambda_prior_prec,
			     Lambda_ineq, theta_ineq, *k0, *k1, *c0, *d0,
			     *c1, *d1, 0, 0);
      z = funval - stream.rexp(1.0);
      if (*method_step == 1){
	StepOut(&delta1_logfcd, X, Lambda, theta, delta0, delta1, 
		Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, 0, 0, z, 
		w, p, stream, L, R, param);
	delta1 = shrinkageStep(&delta1_logfcd, X, Lambda, theta,
			       delta0, delta1, Lambda_prior_mean, 
			       Lambda_prior_prec, Lambda_ineq, theta_ineq,
			       *k0, *k1, *c0, *d0, *c1, *d1, 0, 0,
			       z, w, stream, L, R, param);
      }
      else{
	doubling(&delta1_logfcd, X, Lambda, theta, delta0, delta1, 
		 Lambda_prior_mean, Lambda_prior_prec, Lambda_ineq,
		 theta_ineq, *k0, *k1, *c0, *d0, *c1, *d1, 0, 0, z, 
		 w, p, stream, L, R, param);
	delta1 = shrinkageDoubling(&delta1_logfcd, X, Lambda, theta,
				   delta0, delta1, Lambda_prior_mean, 
				   Lambda_prior_prec, Lambda_ineq, theta_ineq,
				   *k0, *k1, *c0, *d0, *c1, *d1, 0, 0,
				   z, w, stream, L, R, param);
      }
    
      // print results to screen
      if (verbose[0] > 0 && iter % verbose[0] == 0){
	Rprintf("\n\nMCMCirtKdRob iteration %i of %i \n", (iter+1), tot_iter);
      }
        
      // store results
      if ((iter >= burnin[0]) && ((iter % thin[0]==0))) {      
      
	// store Lambda
	if (storeitem[0]==1){
	  //Matrix<double> Lambda_store_vec = reshape(Lambda,1,K*D);
	  //for (int l=0; l<K*D; ++l)
	  //  Lambda_store(count, l) = Lambda_store_vec[l];
	  rmview(Lambda_store(count, _)) = Lambda;
	}
      
	// store theta
	if (storeability[0]==1){
	  //Matrix<double> theta_store_vec = reshape(theta, 1, N*D);
	  //for (int l=0; l<N*D; ++l)
	  //  theta_store(count, l) = theta_store_vec[l];
	  rmview(theta_store(count, _)) = theta;
	}
      
	// store delta0 and delta1
	delta0_store[count] = delta0;
	delta1_store[count] = delta1;

	count++;
      }
    
      // allow user interrupts
      R_CheckUserInterrupt();    
    } // end MCMC loop
  
 



    delete [] K_array;
    delete [] N_array;
    delete [] D_array;
    delete [] Dm1_array;
    delete [] K_inds_perm;
    delete [] N_inds_perm;
    delete [] D_inds_perm;
    delete [] Dm1_inds_perm;



  
    // return output
    Matrix<double> output = delta0_store;
    output = cbind(output, delta1_store);  
    if (*storeitem == 1){
      output = cbind(output, Lambda_store);
    }
    if(*storeability == 1) {
      output = cbind(output, theta_store);
    }

    int size = *samplerow * *samplecol;
    for (int i=0; i<size; ++i)
      sampledata[i] = output[i];
    

}







extern "C"{

  // function called by R to fit model
  void
  irtKdRobpost (double* sampledata, const int* samplerow, 
		const int* samplecol,
		const int* Xdata, const int* Xrow, const int* Xcol,
		const int* burnin, const int* mcmc,  const int* thin,
		const int *uselecuyer, const int *seedarray,
		const int *lecuyerstream, const int* verbose, 
		const int* method_step,
		const double* theta_w, const int* theta_p,
		const double* lambda_w, const int* lambda_p,
		const double* delta0_w, const int* delta0_p,
		const double* delta1_w, const int* delta1_p,		
		const double * delta0start, const double* delta1start,
		const double* Lamstartdata, const int* Lamstartrow, 
		const int* Lamstartcol, 
		const double* thetstartdata, const int* thetstartrow, 
		const int* thetstartcol, 
		const double* Lameqdata, const int* Lameqrow, 
		const int* Lameqcol,
		const double* Lamineqdata, const int* Lamineqrow, 
		const int* Lamineqcol,
		const double* theteqdata, const int* theteqrow, 
		const int* theteqcol,
		const double* thetineqdata, const int* thetineqrow, 
		const int* thetineqcol,
		const double* Lampmeandata, const int* Lampmeanrow, 
		const int* Lampmeancol,
		const double* Lampprecdata, const int* Lampprecrow,
		const int* Lamppreccol, 
		const double* k0, const double* k1,
		const double* c0, const double* c1,
		const double* d0, const double* d1,
		const int* storeitem,
		const int* storeability
		) {
    
    // put together matrices
    const Matrix<int> X(*Xrow, *Xcol, Xdata);
    Matrix<double> Lambda(*Lamstartrow, *Lamstartcol, Lamstartdata);
    Matrix<double> theta(*thetstartrow, *thetstartcol, thetstartdata);
    const Matrix<double> Lambda_eq(*Lameqrow, *Lameqcol, Lameqdata);
    const Matrix<double> Lambda_ineq(*Lamineqrow, *Lamineqcol, 
				     Lamineqdata);
    const Matrix<double> theta_eq(*theteqrow, *theteqcol, 
				  theteqdata);
    const Matrix<double> theta_ineq(*thetineqrow, *thetineqcol, 
				    thetineqdata);
    const Matrix<double> Lambda_prior_mean(*Lampmeanrow, 
					   *Lampmeancol, 
					   Lampmeandata);
    const Matrix<double> Lambda_prior_prec(*Lampprecrow, 
					   *Lamppreccol, 
					   Lampprecdata);  
    
  
   

    MCMCPACK_PASSRNG2MODEL(MCMCirtKdRob_impl, X, Lambda, theta, 
			   Lambda_eq, Lambda_ineq, theta_eq,
			   theta_ineq, Lambda_prior_mean,
			   Lambda_prior_prec, burnin, mcmc, thin,
			   verbose, 
			   method_step,
			   theta_w, theta_p,
			   lambda_w, lambda_p,
			   delta0_w, delta0_p,
			   delta1_w, delta1_p,
			   delta0start, delta1start,
			   k0,  k1,
			   c0,  c1,
			   d0,  d1,
			   storeitem,
			   storeability,
			   sampledata, samplerow, 
			   samplecol			  
			   );

  
  }

}





#endif
