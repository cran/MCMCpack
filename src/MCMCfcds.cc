// MCMCfcds.cc contains definitions for a number of functions that 
// produce draws from full conditional distributions. 
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
// modified to meet the new developer spec KQ 6/18/2004
// update to new Scythe version ADM 7/24/2004

#ifndef MCMCFCDS_CC
#define MCMCFCDS_CC

#include "rng.h"
#include "distributions.h"
#include "stat.h"
#include "smath.h"
#include "la.h"
#include "ide.h"
#include "error.h"
#include "smath.h"

namespace SCYTHE {

  // linear regression with Gaussian errors beta draw 
  // (multivariate Normal prior)
  // regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)
  // XpX is X'X
  // XpY is X'y
  // b0 is the prior mean of beta
  // B0 is the prior precision (the inverse variance) of beta
  Matrix<double> 
  NormNormregress_beta_draw (const Matrix <double> &XpX, 
			     const Matrix <double> &XpY,
			     const Matrix <double> &b0, 
			     const Matrix <double> &B0,
			     const double &sigma2,
			     rng *stream){
    
    // this function gets the cross-product matrix X'X and the matrix X'Y
    // to minimize the amount of computation within the function
    const int k = XpX.cols ();
    const double sig2_inv = 1.0 / sigma2;
    const Matrix <double> sig_beta = invpd (B0 + XpX * sig2_inv);
    const Matrix <double> C = cholesky (sig_beta);
    const Matrix <double> betahat = sig_beta * 
      gaxpy(B0, b0, XpY*sig2_inv);

    return( gaxpy(C, stream->rnorm(k,1), betahat) );
  }
  
  // linear regression with Gaussian errors sigma2 draw 
  // (inverse-Gamma  prior)
  // regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)
  // c0/2 is the prior shape parameter for sigma2
  // d0/2 is the prior scale parameter for sigma2 
  double
  NormIGregress_sigma2_draw (const Matrix <double> &X, 
			     const Matrix <double> &Y,
			     const Matrix <double> &beta, 
			     const double& c0,
			     const double& d0,
			     rng *stream){
    
    const Matrix <double> e = gaxpy(X, (-1*beta), Y);
    const Matrix <double> SSE = crossprod (e); 
    const double c_post = (c0 + X.rows ()) * 0.5;
    const double d_post = (d0 + SSE[0]) * 0.5;

    return  stream->rigamma (c_post, d_post);   
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
  void 
  NormNormfactanal_phi_draw(Matrix<double> &phi,  
			    const Matrix<double> &F0, 
			    const Matrix<double> &Lambda,
			    const Matrix<double> &Psi_inv,
			    const Matrix<double> &X,
			    const int& N, const int& D,
			    rng *stream){
    // If Psi_inv is *not* diagonal then use:
    // Matrix<double> phi_post_var = invpd(F0 + t(Lambda) * Psi_inv * 
    //                                     Lambda);
    //Instead of the following 2 lines: 
    const Matrix<double> AAA = SCYTHE::sqrt(Psi_inv) * Lambda;
    const Matrix<double> phi_post_var = invpd(F0 + crossprod(AAA)); 
    
    const Matrix<double> phi_post_C = cholesky(phi_post_var);
    for (int i=0; i<N; ++i){
      const Matrix<double> phi_post_mean = phi_post_var * 
	(t(Lambda) * Psi_inv * t(X(i,_)));
      const Matrix<double> phi_samp = gaxpy(phi_post_C, stream->rnorm(D, 1), 
				      phi_post_mean); 
      for (int j=0; j<D; ++j)
	phi(i,j) = phi_samp[j];
    }    
  }
   
  // Psi_inv assumed diagnonal
  // this function draws the factor loading matrix
  void 
  NormNormfactanal_Lambda_draw(Matrix<double>& Lambda, 
			       const Matrix<double> &Lambda_free_indic,
			       const Matrix<double> &Lambda_prior_mean,
			       const Matrix<double> &Lambda_prior_prec,
			       const Matrix<double> &phi,
			       const Matrix<double> &X,
			       const Matrix<double> &Psi_inv,
			       const Matrix<double> &Lambda_ineq,
			       const int& D, const int& K,
			       rng *stream) {
    
    for (int i=0; i<K; ++i){
      const Matrix<double> free_indic = t(Lambda_free_indic(i,_));
      const Matrix<double> not_free_indic = (free_indic-1)*-1;
      if (sumc(free_indic)[0] > 0 && 
	  sumc(not_free_indic)[0] > 0){ // both constrnd & unconstrnd
	const Matrix<double> phifree_i =  t(selif(t(phi), free_indic));
	const Matrix<double> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
					   free_indic); // prior mean
	const Matrix<double> hold = selif(t(Lambda_prior_prec(i,_)), 
				    free_indic);
	Matrix<double> sig2lamfree_inv_i = 
	  eye<double>(hold.rows());   // prior prec
	for (int j=0; j<(hold.rows()); ++j)
	  sig2lamfree_inv_i(j,j) = hold[j];
	const Matrix<double> Lambdacon_i = 
	  selif(t(Lambda(i,_)), not_free_indic);
	const Matrix<double> phicon_i  = t(selif(t(phi), not_free_indic));
	const Matrix<double> newX_i = gaxpy((-1.0*phicon_i), Lambdacon_i,
					    X(_,i));
	const Matrix<double> Lam_post_var = invpd(sig2lamfree_inv_i + 
						  Psi_inv(i,i) * 
						  crossprod(phifree_i)); 
	const Matrix<double> Lam_post_C = cholesky(Lam_post_var);
	const Matrix<double> Lam_post_mean = Lam_post_var * 
	  (sig2lamfree_inv_i * mulamfree_i + Psi_inv(i,i) * 
	   t(phifree_i) * newX_i);
	
	Matrix<double> Lambdafree_i = 
	  gaxpy(Lam_post_C, stream->rnorm(hold.rows(), 1), Lam_post_mean);
	
	// check to see if inequality constraints hold
	const Matrix<double> Lambda_ineq_vec = Lambda_ineq(i,_);
	double ineq_holds = 0;
	int Lam_count = 0;
	for (int j=0; j<D; ++j){
	  if (free_indic[j]==1)
	    ineq_holds = std::min(ineq_holds, 
				  Lambda_ineq_vec[j] * 
				  Lambdafree_i[Lam_count]);
	  ++Lam_count;
	}
	while (ineq_holds < 0){
	  Lambdafree_i = 
	    gaxpy(Lam_post_C, stream->rnorm(hold.rows(), 1), Lam_post_mean);
	  Lam_count = 0;
	  double test = 0;
	  for (int j=0; j<D; ++j){
	    if (free_indic[j]==1){
	      Matrix<double> prodcheck = 
		Lambda_ineq_vec[j]*Lambdafree_i[Lam_count];    
	      test = std::min(test, prodcheck[0]); 	      
	      ++Lam_count;
	    }
	  }
	  ineq_holds = test;
	}
	
	// put draw into Lambda 
	Lam_count = 0;
	for (int j=0; j<D; ++j){
	  if (free_indic[j] == 1){
	    Lambda(i,j) = Lambdafree_i[Lam_count];
	    ++Lam_count;
	  }
	}
      }
      else  if (sumc(free_indic)[0] > 0){ // just unconstrained
	const Matrix<double> phifree_i =  t(selif(t(phi), free_indic));
	const Matrix<double> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
						 free_indic); // prior mean
	const Matrix<double> hold = selif(t(Lambda_prior_prec(i,_)), 
					  free_indic);
	Matrix<double> sig2lamfree_inv_i = 
	  eye<double>(hold.rows());  // prior prec
	for (int j=0; j<hold.rows(); ++j)
	  sig2lamfree_inv_i(j,j) = hold[j];
	const Matrix<double> Lam_post_var = invpd(sig2lamfree_inv_i + 
						  Psi_inv(i,i) * 
						  crossprod(phifree_i)); 
	const Matrix<double> Lam_post_C = cholesky(Lam_post_var);
	const Matrix<double> Lam_post_mean = Lam_post_var * 
	  (sig2lamfree_inv_i * mulamfree_i + Psi_inv(i,i) * 
	   t(phifree_i) * X(_,i));
	Matrix<double> Lambdafree_i = 
	  gaxpy(Lam_post_C, stream->rnorm(hold.rows(), 1), Lam_post_mean);
	
	// check to see if inequality constraints hold
	Matrix<double> Lambda_ineq_vec = Lambda_ineq(i,_);
	double ineq_holds = 0;
	for (int j=0; j<D; ++j){
	  ineq_holds = 
	    std::min(ineq_holds, Lambda_ineq_vec[j]*Lambdafree_i[j]);
	}
	while (ineq_holds < 0){
	  Lambdafree_i = 
	    gaxpy(Lam_post_C, stream->rnorm(hold.rows(), 1), Lam_post_mean);
	  double test = 0;
	  for (int j=0; j<D; ++j){
	    //if (free_indic[j]==1)
	    double prodcheck = Lambda_ineq_vec[j]*Lambdafree_i[j];
	    test = std::min(test, prodcheck); 
	  }
	  ineq_holds = test;
	}
	
	// put draw into Lambda
	for (int j=0; j<D; ++j){
	  Lambda(i,j) = Lambdafree_i[j];
	}
      }	
    }      
    
    //    return(Lambda);
  }


  // samples the Psi matrix for a Normal theory factor model with IG 
  // prior on diag elements of Psi
  void 
  NormIGfactanal_Psi_draw(Matrix<double> &Psi, const Matrix<double> &X,
			  const Matrix<double> &phi, 
			  const Matrix<double> &Lambda,
			  const Matrix<double> &a0,
			  const Matrix<double> &b0,
			  const int& K, const int&N,
			  rng *stream){
    for (int i=0; i<K; ++i){
      const Matrix<double> epsilon = gaxpy(phi, -1*(t(Lambda(i,_))), X(_,i));
      const Matrix<double>  SSE = crossprod(epsilon);
      const double a1 = (a0[i] + N)*0.5;
      const double b1 = (b0[i] + SSE[0])*0.5;
      Psi(i,i) = stream->rigamma(a1, b1);
    }    
  }


  // update latent data for standard item response models
  // only works for 1 dimensional case
  void irt_Z_update1 (Matrix<double> &Z, const Matrix<int>& X, 
		      const Matrix<double> &theta, 
		      const Matrix<double> &eta, rng *stream) {
    // define constants
    const int J = theta.rows();
    const int K = eta.rows();
    
    // perform update from truncated Normal / standard Normals
    for (int i=0; i<J; ++i) {
      for (int j=0; j<K; ++j){
	const double Z_mean = -eta(j,0) + theta[i] * eta(j,1);
	if (X(i,j) == 1){
	  Z(i,j) = stream->rtbnorm_combo(Z_mean, 1.0, 0);
	}
	else if (X(i,j) == 0){
	  Z(i,j) = stream->rtanorm_combo(Z_mean, 1.0, 0);
	}
	else {
	  Z(i,j) = stream->rnorm(Z_mean, 1.0);
	}
      }
    }
  }
  
  // update item (case, roll call)  parameters for item response model
  // note: works only for one-dimensional case
  void
  irt_eta_update1 (Matrix<double> & eta, const Matrix<double> & Z,
		   const Matrix<double> & theta, 
		   const Matrix<double> & ab0,
		   const Matrix<double> & AB0,  
		   rng *stream) {
    
    // define constants
    const int J = theta.rows();
    const int K = Z.cols();
    const Matrix<double> AB0ab0 = AB0 * ab0;

    // perform update 
    const Matrix<double> Ttheta_star = t(cbind(-1.0*ones<double>(J,1),theta)); // only needed for option 2
    const Matrix<double> tpt(2,2);
    for (int i=0; i<J; ++i){
      const double theta_i = theta[i];
      tpt(0,1) -= theta_i;
      tpt(1,1) += std::pow(theta_i, 2.0);
    }
    tpt(1,0) = tpt(0,1);
    tpt(0,0) = J;
    const Matrix<double> eta_post_var = invpd(tpt + AB0);
    const Matrix<double> eta_post_C = cholesky(eta_post_var);
    
    for (int k=0; k<K; ++k){    
      const Matrix<double> TZ(2, 1);
      for (int j=0; j<J; ++j){
	TZ[0] -= Z(j,k);
	TZ[1] += Z(j,k) * theta[j];
      }
      const Matrix<double> eta_post_mean = eta_post_var * 
	(TZ + AB0ab0);     
      const Matrix<double> new_eta = gaxpy(eta_post_C, 
					   stream->rnorm(2, 1), 
					   eta_post_mean);
	 eta(k,0) = new_eta[0];
	 eta(k,1) = new_eta[1];
    }
  }
  
  // update ability parameters (ideal points) for one dimensional 
  // item response model 
  // note: works only for one-dimensional case
  void
  irt_theta_update1 (Matrix<double>& theta, const Matrix<double> & Z, 
		     const Matrix<double> & eta, 
		     const double& t0, const double& T0, 
		     const Matrix<double>& theta_eq,
		     const Matrix<double>& theta_ineq,
		     rng *stream) {
    
    const int J = Z.rows();
    const int K = Z.cols();

    // perform update from multivariate Normal
    const double T0t0 = T0*t0;
    const Matrix<double> alpha = eta(_, 0);
    const Matrix<double> beta =  eta(_, 1);
    const Matrix<double> tbeta = t(beta);  // only neede for option 2
    const Matrix<double> talpha = t(alpha); // only needed for option 2

    // calculate the posterior variance outside the justice specific loop
    double theta_post_var = T0;
    for (int i=0; i<K; ++i)
      theta_post_var += std::pow(beta[i], 2.0);
    theta_post_var = 1.0/theta_post_var;
    const double theta_post_sd = std::sqrt(theta_post_var);

    // sample for each justice
    for (int j=0; j<J; ++j) {
      // no equality constraints 
      if (theta_eq[j] == -999){
	double betaTZjalpha = 0;
	for (int k=0; k<K; ++k)
	  betaTZjalpha += beta[k] * (Z(j,k) + alpha[k]); 
	const double theta_post_mean = theta_post_var * 
	  (T0t0 + betaTZjalpha);
	
	if (theta_ineq[j] == 0){ // no inequality constraint
	  theta[j] = theta_post_mean + stream->rnorm(0.0, theta_post_sd); 
	}
	else if (theta_ineq[j] > 0){ // theta[j] > 0
	  theta[j] = stream->rtbnorm_combo(theta_post_mean, 
					   theta_post_var, 0);  
	}
	else { // theta[j] < 0
	  theta[j] = stream->rtanorm_combo(theta_post_mean, 
					   theta_post_var, 0);  	  
	}
      }
      else { // equality constraints
	theta[j] = theta_eq[j];
      }
    } 
    
  }
  
}// end namespace SCYTHE
#endif
