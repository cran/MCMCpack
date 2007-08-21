// MCMCfcds.h is the header file for MCMCfcds.cc. It contains declarations 
// for a number of functions that produce draws from full conditional 
// distributions. 
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
// KQ 6/10/2004
// DBP 7/01/2007 [ported to scythe 1.0.x (partial)]


#ifndef MCMCFCDS_H
#define MCMCFCDS_H

#include "matrix.h"
#include "rng.h"
#include "stat.h"
#include "smath.h"
#include "ide.h"
#include "la.h"
#include "distributions.h"

#include <iostream>
using namespace std;

using namespace scythe;

// linear regression with Gaussian errors beta draw 
// (multivariate Normal prior)
// regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)
// XpX is X'X
// XpY is X'y
// b0 is the prior mean of beta
// B0 is the prior precision (the inverse variance) of beta
template <typename RNGTYPE>
Matrix<double> 
NormNormregress_beta_draw (const Matrix<>& XpX, const Matrix<>& XpY,
         const Matrix<>& b0, const Matrix<>& B0, double sigma2,
         rng<RNGTYPE>& stream)
{
  // this function gets the cross-product matrix X'X and the matrix X'Y
  // to minimize the amount of computation within the function
  const unsigned int k = XpX.cols ();
  const double sig2_inv = 1.0 / sigma2;
  const Matrix<> sig_beta = invpd (B0 + XpX * sig2_inv);
  const Matrix<> C = cholesky (sig_beta);
  const Matrix<> betahat = sig_beta * gaxpy(B0, b0, XpY*sig2_inv);

  return( gaxpy(C, stream.rnorm(k,1, 0, 1), betahat) );
}

// linear regression with Gaussian errors sigma2 draw 
// (inverse-Gamma  prior)
// regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)
// c0/2 is the prior shape parameter for sigma2
// d0/2 is the prior scale parameter for sigma2 
template <typename RNGTYPE>
double
NormIGregress_sigma2_draw (const Matrix <> &X, const Matrix <> &Y,
         const Matrix <> &beta, double c0, double d0,
         rng<RNGTYPE>& stream)
{
  
  const Matrix <> e = gaxpy(X, (-1*beta), Y);
  const Matrix <> SSE = crossprod (e); 
  const double c_post = (c0 + X.rows ()) * 0.5;
  const double d_post = (d0 + SSE[0]) * 0.5;

  return  stream.rigamma (c_post, d_post);   
}

// update latent data for standard item response models
// only works for 1 dimensional case
template <typename RNGTYPE>
void irt_Z_update1 (Matrix<>& Z, const Matrix<int>& X, 
        const Matrix<>& theta, const Matrix<>& eta, rng<RNGTYPE>& stream)
{
  // define constants
  const unsigned int J = theta.rows();
  const unsigned int K = eta.rows();
  
  // perform update from truncated Normal / standard Normals
  for (unsigned int i = 0; i < J; ++i) {
    for (unsigned int j = 0; j < K; ++j){
      const double Z_mean = -eta(j,0) + theta(i) * eta(j,1);
      if (X(i,j) == 1) {
        Z(i,j) = stream.rtbnorm_combo(Z_mean, 1.0, 0);
      } else if (X(i,j) == 0) {
        Z(i,j) = stream.rtanorm_combo(Z_mean, 1.0, 0);
      } else {
        Z(i,j) = stream.rnorm(Z_mean, 1.0);
      }
    }
  }
}

// update item (case, roll call)  parameters for item response model
// note: works only for one-dimensional case
template <typename RNGTYPE>
void irt_eta_update1 (Matrix<>& eta, const Matrix<>& Z,
     const Matrix<>& theta, const Matrix<>& AB0, const Matrix<>& AB0ab0, 
     rng<RNGTYPE>& stream)
{
  
  // define constants
  const unsigned int J = theta.rows();
  const unsigned int K = Z.cols();

  // perform update 
  //const Matrix<double> Ttheta_star = t(cbind(-1.0*ones<double>(J,1),theta)); // only needed for option 2
  const Matrix<> tpt(2,2);
  for (unsigned int i = 0; i < J; ++i) {
    const double theta_i = theta(i);
    tpt(0,1) -= theta_i;
    tpt(1,1) += std::pow(theta_i, 2.0);
  }
  tpt(1,0) = tpt(0,1);
  tpt(0,0) = J;
  const Matrix<> eta_post_var = invpd(tpt + AB0);
  const Matrix<> eta_post_C = cholesky(eta_post_var);
  
  for (unsigned int k = 0; k < K; ++k) {    
    const Matrix<> TZ(2, 1);
    for (unsigned int j = 0; j < J; ++j) {
      TZ[0] -= Z(j,k);
      TZ[1] += Z(j,k) * theta[j];
    }
    const Matrix<> eta_post_mean = eta_post_var * (TZ + AB0ab0);     
    const Matrix<> new_eta = gaxpy(eta_post_C, stream.rnorm(2, 1, 0, 1), 
                                   eta_post_mean);
     eta(k,0) = new_eta(0);
     eta(k,1) = new_eta(1);
  }
}

// update ability parameters (ideal points) for one dimensional 
// item response model 
// note: works only for one-dimensional case
template <typename RNGTYPE>
void irt_theta_update1 (Matrix<>& theta, const Matrix<>& Z, 
    const Matrix<>& eta, double t0, double T0, const Matrix<>& theta_eq,
    const Matrix<>& theta_ineq, rng<RNGTYPE>& stream)
{
  
  const unsigned int J = Z.rows();
  const unsigned int K = Z.cols();

  // perform update from multivariate Normal
  const double T0t0 = T0*t0;
  const Matrix<> alpha = eta(_, 0);
  const Matrix<> beta =  eta(_, 1);
  //const Matrix<double> tbeta = t(beta);  // only neede for option 2
  //const Matrix<double> talpha = t(alpha); // only needed for option 2

  // calculate the posterior variance outside the justice specific loop
  double theta_post_var = T0;
  for (unsigned int i = 0; i < K; ++i)
    theta_post_var += std::pow(beta(i), 2.0);
  theta_post_var = 1.0 / theta_post_var;
  const double theta_post_sd = std::sqrt(theta_post_var);

  // sample for each justice
  for (unsigned int j = 0; j < J; ++j) {
    // no equality constraints 
    if (theta_eq(j) == -999) {
      double betaTZjalpha = 0;
      for (unsigned int k = 0; k < K; ++k)
        betaTZjalpha += beta(k) * (Z(j,k) + alpha(k)); 
        const double theta_post_mean = theta_post_var 
          * (T0t0 + betaTZjalpha);

        if (theta_ineq(j) == 0) { // no inequality constraint
          theta(j) = theta_post_mean + stream.rnorm(0.0, theta_post_sd); 
        } else if (theta_ineq(j) > 0) { // theta[j] > 0
          theta(j) = stream.rtbnorm_combo(theta_post_mean, theta_post_var,
              0);  
        } else { // theta[j] < 0
          theta(j) = stream.rtanorm_combo(theta_post_mean, theta_post_var,
              0);  	  
        }
    } else { // equality constraints
      theta(j) = theta_eq(j);
    }
  } 
}






// factor analysis model with normal mean 0, precision F0 prior on 
// factor scores
// X follows a multivariate normal distribution
// Lambda is the matrix of factor loadings
// Psi_inv is the inverse of the uniqueness matrix
// N is number of observations
// D is the number of factors
// this function draws the factor scores
//
//                      IMPORTANT
// ***********Psi_inv IS ASSUMED TO DIAGONAL ***********
template <typename RNGTYPE>
void 
NormNormfactanal_phi_draw(Matrix<> &phi,  
			  const Matrix<> &F0, 
			  const Matrix<> &Lambda,
			  const Matrix<> &Psi_inv,
			  const Matrix<> &X,
			  const int& N, const int& D,
			  rng <RNGTYPE>& stream){
  // If Psi_inv is *not* diagonal then use:
  // Matrix<double> phi_post_var = invpd(F0 + t(Lambda) * Psi_inv * 
  //                                     Lambda);
  //Instead of the following 2 lines: 
  const Matrix<> AAA = scythe::sqrt(Psi_inv) * Lambda;
  const Matrix<> phi_post_var = invpd(F0 + crossprod(AAA)); 
  const Matrix<> phi_post_C = cholesky(phi_post_var);
  for (int i=0; i<N; ++i){
    const Matrix<> phi_post_mean = phi_post_var * 
      (t(Lambda) * Psi_inv * t(X(i,_)));
    const Matrix<> phi_samp = gaxpy(phi_post_C, stream.rnorm(D, 1, 0.0, 1.0), 
				    phi_post_mean); 
    for (int j=0; j<D; ++j)
      phi(i,j) = phi_samp(j);
  }    
}


// samples the Psi matrix for a Normal theory factor model with IG 
// prior on diag elements of Psi
template <typename RNGTYPE>
void 
NormIGfactanal_Psi_draw(Matrix<> &Psi, const Matrix<> &X,
			const Matrix<> &phi, 
			const Matrix<> &Lambda,
			const Matrix<> &a0,
			const Matrix<> &b0,
			const int& K, const int& N,
			rng<RNGTYPE>& stream){
  for (int i=0; i<K; ++i){
    const Matrix<double> epsilon = gaxpy(phi, -1*(t(Lambda(i,_))), X(_,i));
    const Matrix<double>  SSE = crossprod(epsilon);
    const double a1 = (a0[i] + N)*0.5;
    const double b1 = (b0[i] + SSE[0])*0.5;
    Psi(i,i) = stream.rigamma(a1, b1);
  }    
}



// Psi_inv assumed diagnonal
// this function draws the factor loading matrix
template <typename RNGTYPE>
void 
NormNormfactanal_Lambda_draw(Matrix<>& Lambda, 
			     const Matrix<bool> &Lambda_free_indic,
			     const Matrix<> &Lambda_prior_mean,
			     const Matrix<> &Lambda_prior_prec,
			     const Matrix<> &phi,
			     const Matrix<> &X,
			     const Matrix<> &Psi_inv,
			     const Matrix<> &Lambda_ineq,
			     const unsigned int D, const unsigned int K, 
           rng<RNGTYPE>& stream)
{

  for (unsigned int i = 0; i < K; ++i) {
    const Matrix<bool> free_indic = t(Lambda_free_indic(i,_));
    // end replacement

    if (sumc(free_indic)(0) > 0 && sumc(! free_indic)(0) > 0) { 
      // both constrnd & unconstrnd
      const Matrix<> phifree_i =  t(selif(t(phi), free_indic));
      const Matrix<> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
					 free_indic); // prior mean
      const Matrix<> hold = selif(t(Lambda_prior_prec(i,_)), 
				  free_indic); 
      Matrix<> sig2lamfree_inv_i 
	= eye<double>(hold.rows()); // prior prec
				
      for (unsigned int j = 0; j < (hold.rows()); ++j)
	sig2lamfree_inv_i(j,j) = hold[j];

      const Matrix<> Lambdacon_i = selif(t(Lambda(i,_)), ! free_indic); 
      const Matrix<> phicon_i  = t(selif(t(phi), ! free_indic));
      const Matrix<> newX_i = gaxpy((-1.0*phicon_i), Lambdacon_i, 
				    X(_,i));
      const Matrix<> Lam_post_var = invpd(sig2lamfree_inv_i + 
					  Psi_inv(i,i) * crossprod(phifree_i));
      const Matrix<> Lam_post_C = cholesky(Lam_post_var);
      const Matrix<> Lam_post_mean = Lam_post_var * 
	(sig2lamfree_inv_i * mulamfree_i + Psi_inv(i,i) * 
	 t(phifree_i) * newX_i);
				

      Matrix<> Lambdafree_i = 
	gaxpy(Lam_post_C, stream.rnorm(hold.rows(), 1, 0, 1),Lam_post_mean);

      
      
      // check to see if inequality constraints hold
      const Matrix<> Lambda_ineq_vec = Lambda_ineq(i,_);

      double ineq_holds = 0;
      int Lam_count = 0;
      for (unsigned int j = 0; j < D; ++j) {
	if (free_indic(j)){ 
	  ineq_holds = std::min(ineq_holds, Lambda_ineq_vec(j) * 
				Lambdafree_i(Lam_count)); 
	  ++Lam_count;
	}
      }

      while (ineq_holds < 0) {
	Lambdafree_i = gaxpy(Lam_post_C, stream.rnorm(hold.rows(), 1, 0, 1),
			     Lam_post_mean);
	Lam_count = 0;
	double test = 0;
	for (unsigned int j = 0; j < D; ++j) {
	  if (free_indic(j) == 1) {
	    Matrix<> prodcheck = Lambda_ineq_vec(j) 
	      * Lambdafree_i(Lam_count);
	    test = std::min(test, prodcheck(0));
	    ++Lam_count;
	  }
	} 
	ineq_holds = test; 
      }

      // put draw into Lambda
      Lam_count = 0;
      for (unsigned int j = 0; j < D; ++j) { 
	if (free_indic(j) == 1) {
	  Lambda(i,j) = Lambdafree_i(Lam_count);
	  ++Lam_count; 
	} 
      }
    } else  if (sumc(free_indic)(0) > 0) { // just unconstrained
      const Matrix<> phifree_i =  t(selif(t(phi), free_indic));
      const Matrix<> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
					 free_indic); // prior mean
      const Matrix<> hold = selif(t(Lambda_prior_prec(i,_)), free_indic);
      Matrix<> sig2lamfree_inv_i = eye<double>(hold.rows()); // prior prec
			
      for (unsigned int j = 0; j < hold.rows(); ++j)
	sig2lamfree_inv_i(j,j) = hold(j); 
			
      const Matrix<> Lam_post_var = invpd(sig2lamfree_inv_i + 
					  Psi_inv(i,i) * crossprod(phifree_i));
      const Matrix<> Lam_post_C = cholesky(Lam_post_var);
      const Matrix<> Lam_post_mean = Lam_post_var * (sig2lamfree_inv_i 
						     * mulamfree_i + 
						     Psi_inv(i,i) * 
						     t(phifree_i) * X(_,i));
      Matrix<> Lambdafree_i = gaxpy(Lam_post_C, 
				    stream.rnorm(hold.rows(), 1, 0, 1), 
				    Lam_post_mean);


      // check to see if inequality constraints hold
      Matrix<> Lambda_ineq_vec = Lambda_ineq(i,_); 
      double ineq_holds = 0;
      for (unsigned int j = 0; j < D; ++j) { 
	ineq_holds = std::min(ineq_holds, Lambda_ineq_vec(j)
			      * Lambdafree_i(j)); 
      } 
			
      while (ineq_holds < 0) {
	Lambdafree_i = gaxpy(Lam_post_C, stream.rnorm(hold.rows(), 1, 0, 1),
			     Lam_post_mean); 
	double test = 0; 
	for (unsigned int j = 0; j < D; ++j) { 
	  //if (free_indic[j]==1) 
	  double prodcheck = Lambda_ineq_vec[j]*Lambdafree_i[j]; 
	  test = std::min(test, prodcheck);
	} 
				
	ineq_holds = test; 
      }

      // put draw into Lambda 
      for (unsigned int j = 0; j < D; ++j) { 
	Lambda(i,j) = Lambdafree_i(j); 
      }
    }
  }      
  //    return(Lambda);
}

















#endif



