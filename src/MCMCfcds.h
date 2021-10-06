//////////////////////////////////////////////////////////////////////////
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
// KQ 6/10/2004
// DBP 7/01/2007 [ported to scythe 1.0.x (partial)]
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////



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

// Bayesian quantile (including median) regression beta draw 
// (multivariate Normal prior)
// 
// b0 is the prior mean of beta(tau)
// B0 is the prior precision (the inverse variance) of beta(tau)
template <typename RNGTYPE>
Matrix<double> 
ALaplaceNormregress_beta_draw (double tau, const Matrix<>& X, const Matrix<>& Y, const Matrix<>& weights,
         const Matrix<>& b0, const Matrix<>& B0,
         rng<RNGTYPE>& stream)
{

	const unsigned int k = X.cols();
	const unsigned int n_obs = X.rows();
	Matrix<> U(Y);
if (tau!=0.5){
	U -= (1.0-2.0*tau)*weights;
} 
	Matrix<> XtwX(k,k,false);
	Matrix<> XtwU(k,1,false);
	double temp_x = 0.0;
	double temp_u = 0.0;

//Calculate XtwU where w denotes a diagonal matrix with the augmented data (weights) on the diagonal   
  for (unsigned int i=0; i<k; ++i){
    for (unsigned int m=0; m<n_obs; ++m){
       temp_u += X(m,i)*U(m)/weights(m);
    }
    XtwU(i) = temp_u;
    temp_u = 0.0;
  }

//Calculate XtwX
  for (unsigned int i=0; i<k; ++i){
    for (unsigned int j=0; j<(i+1); ++j){
      for (unsigned int m=0; m<n_obs; ++m){
	temp_x += X(m,i)*X(m,j)/weights(m);
      }
      XtwX(i,j) = temp_x;
      XtwX(j,i) = temp_x;
      temp_x = 0.0;
    }
  }

	const Matrix<> var_matrix_beta = invpd(B0+0.5*XtwX);
 	const Matrix<> C = cholesky(var_matrix_beta);
	const Matrix<> betahat = var_matrix_beta*gaxpy(B0,b0,0.5*XtwU);

	return( gaxpy(C, stream.rnorm(k,1, 0, 1), betahat) );
}

// This function draws from the full conditional distribution of the latent random variables (weights) under quantile regression (including median regression) and returns a column vector of those weights.

template <typename RNGTYPE>
Matrix<double>
ALaplaceIGaussregress_weights_draw (const Matrix <> &abse,
         rng<RNGTYPE>& stream)

{
        const Matrix<double> nu_params = pow(abse,-1.0);
	Matrix<> w(abse);
        const unsigned int n_obs = abse.rows();

	// The inverse Gaussian distribution

	for (unsigned int i=0; i<n_obs; ++i){
            double chisq = stream.rchisq(1);
            double nu = nu_params(i);
	    double smallroot = nu*(nu*chisq+1.0-std::sqrt(nu*nu*chisq*chisq+2.0*nu*chisq));
	    unsigned int q = stream.rbern(nu/(nu+smallroot));
	    if (q == 1){
		w(i) = smallroot;
            }
	    else{
		w(i) = nu*nu/smallroot;
	    }
	  }
	return(pow(w,-1.0));
}

/////////////////////////////////////////////////////

// Functions for the quantile regression stochastic search variable selection

struct COV_TRIAL{
Matrix<> Cnew;
bool newtrial;
double logdetminhalf;
};

// updating the indicator variable corresponding to whether a 
// covariate is included in the model, given that it was
// previously absent
template <typename RNGTYPE>
COV_TRIAL
QR_SSVS_covariate_trials_draw_absent(const Matrix<>& C, const Matrix<>& X_gamma, const Matrix<>& U,
				     const Matrix<>& newXcol, unsigned int row_index, const Matrix<>& weights, double pi0, double newlambda, double logolddetminhalf, rng<RNGTYPE>& stream){
  const unsigned int n_obs = U.rows();
  const unsigned int k = C.rows();
  
  //Calculate new row required to update the Cholesky decomposition
  Matrix<> XUXnewtXnew(k+1,1,false);
  double temp_xux1 = 0.0;
  double temp_xux2 = 0.0;
  double temp_xux3 = 0.0;
  
  //Calculate XUXnewtXnew
  for (unsigned int i=0; i<k-1; ++i){
    for (unsigned int m=0; m<n_obs; ++m){
      temp_xux1 += X_gamma(m,i)*newXcol(m)/weights(m);
    }
    XUXnewtXnew(i) = 0.5*temp_xux1;
    temp_xux1 = 0.0;
  }
  for (unsigned int m=0; m<n_obs; ++m){
    temp_xux2 += U(m)*newXcol(m)/weights(m);
    temp_xux3 += newXcol(m)*newXcol(m)/weights(m);
  }
  XUXnewtXnew(k-1) = 0.5*temp_xux2;
  XUXnewtXnew(k) = 0.5*temp_xux3+newlambda;
  
  //Obtain the Cholesky Decomposition of the new matrix formed by adding newXcol onto the right hand side
  
  Matrix<> z(k,1,false);
  for (unsigned int i = 0; i < k; ++i) {
    double sum = 0;
    for (unsigned int j = 0; j < i; ++j) {
      sum += C(i,j) * z(j);
    }
    z(i) = (XUXnewtXnew(i) - sum) / C(i, i);
  }
  double rho = std::sqrt(XUXnewtXnew(k)-crossprod(z)(0));
  
  Matrix<> Cnew(k+1, k+1, true, 0.0);
  Cnew(0,0,k-1,k-1) = C;
  Cnew(k,0,k,k-1) = z;
  Cnew(k,k) = rho;
  
  // Permuting the Cholesky decomposition so that it corresponds to what would be obtained if the X matrix included the covariate
  
  Matrix<> temp(Cnew);
  if (row_index != 0){
    temp(0,0,row_index-1,k) = Cnew(0,0,row_index-1,k);
  }
  temp(row_index,_) = Cnew(k,_);
  temp(row_index+1,0,k,k) = Cnew(row_index,0,k-1,k);
  
  // Givens rotations
  
  Matrix<> Q(2,2,false);
  for (unsigned int i=k; i>row_index; --i)
    {
      double two_norm = std::sqrt(temp(row_index,i)*temp(row_index,i)
				  +temp(row_index,i-1)*temp(row_index,i-1));
      Q(0,0) = temp(row_index,i-1)/two_norm;
      Q(1,0) = temp(row_index,i)/two_norm;
      Q(1,1) = Q(0,0);
      Q(0,1) = -1.0*Q(1,0);
      if (i!=k){
	temp(i+1,i-1,k,i) = temp(i+1,i-1,k,i) * Q;
      }
      double temp2 = temp(i,i-1);
      temp(i,i-1) = Q(0,0)*temp2;
      temp(i,i) = Q(0,1)*temp2;
      if (temp(i,i) < 0){
	temp(i,i,k,i) = -1.0*temp(i,i,k,i);
      }
      temp(row_index,i-1) = two_norm;
      temp(row_index,i) = 0.0;
    }
  Cnew=temp;
  
  //Work out -0.5*log(det(Cnew'Cnew))
  double lognewdetminhalf = 0.0;
  for (unsigned int i=0; i<k; ++i){
    lognewdetminhalf -= std::log(Cnew(i,i));
  }
  
  
  
  double log_g0 = logolddetminhalf-0.5*C(k-1,k-1)*C(k-1,k-1);
  double log_g1 = 0.5*std::log(newlambda)+lognewdetminhalf-0.5*Cnew(k,k)*Cnew(k,k);
  
  double log_ratio = log_g0+std::log(1.0-pi0)-log_g1-std::log(pi0);
  double success_prob = 1.0/(1.0+std::exp(log_ratio));
  bool new_covariate_trial = stream.rbern(success_prob);
  COV_TRIAL result;
  result.newtrial = new_covariate_trial;
  if (new_covariate_trial == false){
    result.Cnew = C;
    result.logdetminhalf = logolddetminhalf;
  }
  else {
    result.Cnew = Cnew;
    result.logdetminhalf = lognewdetminhalf;
  }  
  return result;
}

// updating the indicator variable corresponding to whether a 
// covariate is included in the model, given that it was
// previously present in the model

template <typename RNGTYPE>
COV_TRIAL
QR_SSVS_covariate_trials_draw_present(const Matrix<>& C, unsigned int row_index, unsigned int n_obs, double pi0, double oldlambda, double logolddetminhalf, rng<RNGTYPE>& stream){
  unsigned int k = C.rows();
  
  // Permuting the Cholesky decomposition so that it corresponds to what would be obtained if the X matrix had the covariate in the final column
  
  Matrix<> temp(C);
  if (row_index != 0){
    temp(0,0,row_index-1,k-1) = C(0,0,row_index-1,k-1);
  }
  temp(k-1,_) = C(row_index,_);
  temp(row_index,0,k-2,k-1) = C(row_index+1,0,k-1,k-1);
  
  // Givens rotations
  
  Matrix<> Q(2,2,false);
  for (unsigned int i=row_index; i<k-1; ++i)
    {
      double two_norm = std::sqrt(temp(i,i)*temp(i,i)
				  +temp(i,i+1)*temp(i,i+1));
      Q(0,0) = temp(i,i)/two_norm;
      Q(1,0) = temp(i,i+1)/two_norm;
      Q(1,1) = Q(0,0);
      Q(0,1) = -1.0*Q(1,0);
      if (i!=k-2){
	temp(i+1,i,k-2,i+1) = temp(i+1,i,k-2,i+1) * Q;
      }
      double temp2 = temp(k-1,i);
      temp(k-1,i) = Q(0,0)*temp2;
      temp(k-1,i+1) = Q(0,1)*temp2;
      temp(i,i) = two_norm;
      temp(i,i+1) = 0.0;
    }
  if (temp(k-1,k-1) < 0){
    temp(k-1,k-1) = -1.0*temp(k-1,k-1);
  }
  
  Matrix<> Cnew = temp(0,0,k-2,k-2);
  
  // Work out -1/2*log(det(Cnew'Cnew))
  double lognewdetminhalf = 0.0;
  for (unsigned int i=0; i<k-2; ++i){
    lognewdetminhalf -= std::log(Cnew(i,i));
  }
  
  
  double log_g0 = lognewdetminhalf-0.5*Cnew(k-2,k-2)*Cnew(k-2,k-2);
  double log_g1 = 0.5*std::log(oldlambda)+logolddetminhalf-0.5*C(k-1,k-1)*C(k-1,k-1);
  
  
  double log_ratio = log_g0+std::log(1.0-pi0)-log_g1-std::log(pi0);
  double success_prob = 1.0/(1.0+std::exp(log_ratio));
  bool new_covariate_trial = stream.rbern(success_prob);
  COV_TRIAL result;
  result.newtrial = new_covariate_trial;
  if (new_covariate_trial == false){
    result.Cnew = Cnew;
    result.logdetminhalf = lognewdetminhalf;
  }
  else {
    result.Cnew = C;
    result.logdetminhalf = logolddetminhalf;
  }  
  return result;
}

// update betas using Cholesky decomposition
template <typename RNGTYPE>
Matrix<>
QR_SSVS_beta_draw(const Matrix<>& C, rng<RNGTYPE>& stream){
  unsigned int k = C.rows();
  Matrix<> standnorm = stream.rnorm(k-1,1,0,1);
  Matrix<> z(k-1,1,false);
  z = t(C(k-1,0,k-1,k-2));
  Matrix<> Q = z+standnorm*std::sqrt(2.0);
  Matrix<> result(k-1,1,false);
  for (int i = k-2; i >= 0; --i) {
    double sum = 0;
    for (unsigned int j = i+1; j < k-1; ++j) {
      sum += C(j,i) * result(j);
    }
    result(i) = (Q(i) - sum) / C(i, i);
  }
  return result; 
}

//hyperparameter pi0 updating

template <typename RNGTYPE>
double
QR_SSVS_pi0_draw(unsigned int n_uncert_cov, unsigned int tot_n_uncert_cov, 
double pi0a0, double pi0b0, rng<RNGTYPE>& stream){
    double pi0a1 = pi0a0 + n_uncert_cov;
    double pi0b1 = pi0b0 + tot_n_uncert_cov - n_uncert_cov;
	return(stream.rbeta(pi0a1,pi0b1));	
}

//update latent lambdas

template<typename RNGTYPE>
Matrix<double>
QR_SSVS_lambda_draw(const Matrix<>& beta_red, const Matrix<>& gamma, unsigned int tot_n_cov, unsigned int n_cert_cov, rng<RNGTYPE>& stream)
    {
unsigned int n_uncert_cov = tot_n_cov - n_cert_cov;
Matrix<> newlambda(n_uncert_cov,1,false);

for (unsigned int i=n_cert_cov; i<tot_n_cov; ++i){

  unsigned int j = i-n_cert_cov;
  
  if (gamma(i) == true){
    
    unsigned int col_index = n_cert_cov;
    //obtain column index of betas
    for (unsigned int m=n_cert_cov; m<i; ++m){
      if (gamma(m) == true){
	++col_index;
      }
    }
    
    newlambda(j) = stream.rexp(0.5*(1.0+beta_red(col_index)*beta_red(col_index)));
  }
  
else {
newlambda(j) = stream.rexp(0.5);
}

}
    return newlambda;
    }


/////////////////////////////////////////////////////

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
  for (unsigned int i = 0; i < K; ++i){
    theta_post_var += std::pow(beta(i), 2.0);
  }
  theta_post_var = 1.0 / theta_post_var;
  const double theta_post_sd = std::sqrt(theta_post_var);
  
  // sample for each justice
  for (unsigned int j = 0; j < J; ++j) {
    // no equality constraints 
    if (theta_eq(j) == -999) {
      double betaTZjalpha = 0;
      for (unsigned int k = 0; k < K; ++k){
        betaTZjalpha += beta(k) * (Z(j,k) + alpha(k));
      }
      const double theta_post_mean = theta_post_var*(T0t0 + betaTZjalpha);
      if (theta_ineq(j) == 0) { // no inequality constraint
	theta(j) = theta_post_mean + stream.rnorm(0.0, theta_post_sd); 
      } else if (theta_ineq(j) > 0) { // theta[j] > 0
	theta(j) = stream.rtbnorm_combo(theta_post_mean, theta_post_var, 0);  
      } else { // theta[j] < 0
	theta(j) = stream.rtanorm_combo(theta_post_mean, theta_post_var, 0);    
      }
    } else { // equality constraints
      theta(j) = theta_eq(j);
    }
  } 
}







// update latent Ystar values in paired comparisons model
template <typename RNGTYPE>
void paircompare_Ystar_update (Matrix<>& Ystar,
			       const Matrix<unsigned int>& MD, 
			       const Matrix<>& theta, const Matrix<>& alpha,
			       rng<RNGTYPE>& stream){

  const unsigned int N = MD.rows();
  for (unsigned int i = 0; i < N; ++i){
    const double alpha_i = alpha(MD(i,0));
    const double theta_1 = theta(MD(i,1));
    const double theta_2 = theta(MD(i,2));
    const double mean = alpha_i * (theta_1 - theta_2);
    const bool above = MD(i,1) == MD(i,3); // cand 1 chosen
    const bool below = MD(i,2) == MD(i,3); // cand 2 chosen
    if (above){
      Ystar(i) = stream.rtbnorm_combo(mean, 1.0, 0.0);
    }
    else if (below){
      Ystar(i) = stream.rtanorm_combo(mean, 1.0, 0.0);
    }
    else{
      Ystar(i) = stream.rnorm(mean, 1.0);
    }
  }
}





// update theta values in paired comparisons model
template <typename RNGTYPE>
void paircompare_theta_update (Matrix<>& theta, const Matrix<>& Ystar, 
			       const Matrix<unsigned int>& MD,
			       const Matrix<>& alpha,
			       const Matrix<unsigned int>& theta_n,
			       const Matrix<>& theta_eq,
			       const Matrix<>& theta_ineq,
			       const vector< vector < double* > >& theta_Ystar_ptr,
			       const vector< vector < double* > >& theta_alpha_ptr,
			       const vector< vector < double* > >& theta_theta_ptr,
			       const vector< vector < double > >& theta_sign,
			       rng<RNGTYPE>& stream){
  
  const unsigned int J = theta.rows();  
  for (unsigned int j = 0; j < J; ++j){
    double xx = 0.0;
    double xz = 0.0;
    for (unsigned int i = 0; i < theta_n[j]; ++i){
      const double x_ji = theta_sign[j][i] *  *theta_alpha_ptr[j][i];
      const double z_ji = *theta_Ystar_ptr[j][i] + theta_sign[j][i] *
	*theta_alpha_ptr[j][i] * *theta_theta_ptr[j][i]; 
	 /*
	 The reason for the above plus sign: -(-theta_sign[j][i])=+theta_sign[j][i]. 
	 The other theta's sign in a comparison is the opposite to the theta in scrutiny.
	 */
      xx += x_ji * x_ji;
      xz += x_ji * z_ji;
    }
	// prior for theta is N(0,1), therefore the posterior mean and variance are as follow
	
	const double v = 1.0/(xx + 1.0); 
    const double m = v * (xz);
    // no equality constraints 
    if (theta_eq(j) == -999) {      
      if (theta_ineq(j) == 0) { // no inequality constraint
	theta(j) =  stream.rnorm(m, std::sqrt(v)); 
      } else if (theta_ineq(j) > 0) { // theta[j] > 0
	theta(j) = stream.rtbnorm_combo(m, v, 0);  
      } else { // theta[j] < 0
	theta(j) = stream.rtanorm_combo(m, v, 0);  	  
      }
    } else { // equality constraints
      theta(j) = theta_eq(j);
    }    
  }
  
}









// update alpha values in paired comparisons model
template <typename RNGTYPE>
void paircompare_alpha_update (Matrix<>& alpha,
			       const Matrix<>& Ystar, 
			       const Matrix<unsigned int>& MD,
			       const Matrix<>& theta,
			       const double& A0,
			       const double& A0a0,
			       const Matrix<unsigned int>& alpha_n,
			       const vector< vector < double* > >& alpha_Ystar_ptr,
			       const vector< vector < double* > >& alpha_theta1_ptr,
			       const vector< vector < double* > >& alpha_theta2_ptr,
			       rng<RNGTYPE>& stream){
  
  const unsigned int I = alpha.rows();  
  for (unsigned int i = 0; i < I; ++i){
    double xx = 0.0;
    double xz = 0.0;
    for (unsigned int j = 0; j < alpha_n[i]; ++j){
      const double x_ji = *alpha_theta1_ptr[i][j] -  *alpha_theta2_ptr[i][j];
      const double z_ji = *alpha_Ystar_ptr[i][j];
      xx += x_ji * x_ji;
      xz += x_ji * z_ji;
    }
	// prior for alpha is N(0,1), therefore the posterior mean and variance are as follow. Same as Bayesian regression
    const double v = 1.0/(xx + A0);
    const double m = v * (xz + A0a0);
    alpha(i) = stream.rnorm(m, std::sqrt(v));    
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


// update ability parameters (ideal points) for one dimensional 
// Hierarchical item response model.
// 2008-11-18 now deals with PX alpha and holds on to thetahat
template <typename RNGTYPE>
void hirt_theta_update1 (Matrix<>& theta, Matrix<>& thetahat, 
			 const Matrix<>& Z, 
			 const Matrix<>& eta, 
			 const Matrix<>& beta, const Matrix<>& Xj,
			 const double& sigma2,
			 const double& alpha,
			 rng<RNGTYPE>& stream)
{
  
  const unsigned int J = Z.rows();
  const unsigned int K = Z.cols();
  // Get level1 prior mean
  const Matrix<double> Xbeta = (Xj * beta);
  
  // a and b are backwards here, different parameterizations
  // common for edu people and us.
  const Matrix<> b = eta(_, 0); // location
  const Matrix<> a = eta(_, 1); // relevance

  // calculate the posterior variance outside the justice specific loop
  const double sig2_inv = 1.0 / sigma2;
  const Matrix<double> apa = crossprod(a);
  const Matrix<double> theta_post_var = scythe::invpd(apa + sig2_inv);
  const double theta_post_sd = std::sqrt(theta_post_var[0]);
  // sample for each justice
  for (unsigned int j = 0; j < J; ++j) {
    thetahat(j) = 0.0;
    for (unsigned int k = 0; k < K; ++k) {
      thetahat(j) += a[k] * ( Z(j,k) + b[k]); // bill contribution
    }
    thetahat(j) += (Xbeta[j] / sigma2); // j prior level1 contribution
    
    thetahat(j) *= theta_post_var[0];
    const double t = thetahat(j) / alpha;
    theta(j) = stream.rnorm( t , theta_post_sd );
  }
}


// update item (case, roll call)  parameters for item response model
// note: works only for one-dimensional case
// updated for PX (alpha) and hold on to etahat 2008-11-18
template <typename RNGTYPE>
void hirt_eta_update1 (Matrix<>& eta, Matrix<>& etahat, const Matrix<>& Z,
			  const Matrix<>& theta, const Matrix<>& AB0, 
			  const Matrix<>& AB0ab0, 
			  const double& alpha,
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
    Matrix<> eta_post_mean = eta_post_var * (TZ + AB0ab0);
    etahat(k,0) = eta_post_mean(0);
    etahat(k,1) = eta_post_mean(1);
    eta_post_mean /= alpha;    
    const Matrix<> new_eta = gaxpy(eta_post_C, stream.rnorm(2, 1, 0, 1), 
                                   (eta_post_mean) );
     eta(k,0) = new_eta(0);
     eta(k,1) = new_eta(1);

     //Rprintf("\n\a: %3.1f,b:%3.1f ",eta(k,0),eta(k,1)); 

  }
}







#endif



