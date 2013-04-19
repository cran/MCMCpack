////////////////////////////////////////////////////////////////////
// MCMCintervention.cc is a C++ code to estimate 
// linear regression changepoint model
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr
//
// Written 03/03/2009
////////////////////////////////////////////////////////////////////

#ifndef MCMCREGRESSCHANGE_CC
#define MCMCREGRESSCHANGE_CC

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


double lndinvgamma_pjh (const double x, const double shape, const double scale){
  double log_density = shape *::log(scale) - lngammafn(shape) - (shape + 1) * ::log(x) - (scale/x);
  return (log_density);
}

template <typename RNGTYPE>
Matrix<double> gaussian_state_fixedBeta_sampler(rng<RNGTYPE>& stream, 
						 const int m, 
						 const Matrix<double>& Y,
						 const Matrix<double>& X,
						 const Matrix<double>& beta,
						 const Matrix<double>& Sigma,
						 const Matrix<double>& P){
  
  const int ns = m + 1;
  const int n = Y.rows();
  Matrix<double> F(n, ns);
  Matrix<double> pr1(ns, 1);
  pr1[0] = 1;
  Matrix<double> py(ns, 1);
  Matrix<double> pstyt1(ns, 1);
  Matrix<int> s(n, 1);                        // holder for state variables
  Matrix<double> ps = Matrix<double>(n, ns);  // holder for state probabilities
  for (int tt=0; tt<n ; ++tt){
    Matrix<double> mu = X(tt,_)*beta; //k by 1 vector
    for (int j = 0; j< ns; ++j){
      py[j]  =  dnorm(Y[tt], mu[0], sqrt(Sigma[j]));
    }
    if (tt==0) pstyt1 = pr1;
    else {
      pstyt1 =  ::t(F(tt-1,_)*P); // make it an ns by 1 matrix
    }
    Matrix<double> unnorm_pstyt = pstyt1%py;
    const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
    for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);
    
  }// end of F matrix filtering
  ps(n-1,_) = F(n-1,_);                     
  s(n-1) = ns;                                
  Matrix<double> pstyn = Matrix<double>(ns, 1);
  double pone = 0.0;
  int tt = n-2;
  while (tt >= 0){
    int st = s(tt+1);
    Matrix<double> Pst_1 = ::t(P(_,st-1)); 
    Matrix<double> unnorm_pstyn = F(tt,_)%Pst_1;
    pstyn = unnorm_pstyn/sum(unnorm_pstyn); 
    if (st==1)   s(tt) = 1;                  
    else{
      pone = pstyn(st-2);
      if(stream.runif() < pone) s(tt) = st-1;
      else s(tt) = st;
    }
    ps(tt,_) = pstyn;
    --tt;
  }// end of while loop
  Matrix<double> Sout(n, ns+1); 
  Sout(_, 0) = s(_,0);
  for (int j = 0; j<ns; ++j){
    Sout(_,j+1) = ps(_, j);
  }
  return Sout;
} // end of state sampler


template <typename RNGTYPE>
Matrix<double> gaussian_state_sampler(rng<RNGTYPE>& stream, 
				      const int m, 
				      const Matrix<double>& Y,
				      const Matrix<double>& X,
				      const Matrix<double>& beta,
				      const Matrix<double>& Sigma,
				      const Matrix<double>& P){
  
  const int ns = m + 1;
  const int n = Y.rows();
  
  // P is a (m+1 by m+1) transition matrix
  // F matrix contains all the information of Pr(st|Yt)
  Matrix<double> F(n, ns);
  Matrix<double> pr1(ns, 1);
  pr1[0] = 1;
  Matrix<double> py(ns, 1);
  Matrix<double> pstyt1(ns, 1);
  Matrix<int> s(n, 1);                        // holder for state variables
  Matrix<double> ps = Matrix<double>(n, ns);  // holder for state probabilities
  
  //
  // Forward sampling: update F matrix
  //
  for (int tt=0; tt<n ; ++tt){
    Matrix<double> mu = X(tt,_)*::t(beta); //k by 1 vector
    for (int j = 0; j< ns; ++j){
      py[j]  =  dnorm(Y[tt], mu[j], sqrt(Sigma[j]));
    }
    if (tt==0) pstyt1 = pr1;
    else {
      pstyt1 =  ::t(F(tt-1,_)*P); // make it an ns by 1 matrix
    }
    Matrix<double> unnorm_pstyt = pstyt1%py;
    const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
    for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);
    
  }// end of F matrix filtering
  
  //
  // Backward recursions: sample s using F and a transition matrix P
  //
  ps(n-1,_) = F(n-1,_);                       // we know last elements of ps and s
  s(n-1) = ns;                                // the last period is s(*n-1) but one down in C location
  
  Matrix<double> pstyn = Matrix<double>(ns, 1);
  double pone = 0.0;
  int tt = n-2;
  while (tt >= 0){
    int st = s(tt+1);
    Matrix<double> Pst_1 = ::t(P(_,st-1)); // prob of being at a previous state
    Matrix<double> unnorm_pstyn = F(tt,_)%Pst_1;
    pstyn = unnorm_pstyn/sum(unnorm_pstyn); // normalize into a prob. density
    if (st==1)   s(tt) = 1;                  // If this is the first period, state should be 1.
    // Otherwise, draw a state from a discrete prob. distribution("pstyn")
    // using the inverse CDF method.
    else{
      pone = pstyn(st-2);
      if(stream.runif() < pone) s(tt) = st-1;
      else s(tt) = st;
    }
    ps(tt,_) = pstyn;
    --tt;
  }// end of while loop
  
  // name and report outputs
  Matrix<double> Sout(n, ns+1); 
  Sout(_, 0) = s(_,0);
  for (int j = 0; j<ns; ++j){
    Sout(_,j+1) = ps(_, j);
  }
  
  return Sout;
  
} // end of state sampler

//////////////////////////////////////////// 
// MCMCintervention random effect changes only  
//////////////////////////////////////////// 
template <typename RNGTYPE>
void MCMCintervention_random_impl(rng<RNGTYPE>& stream, 
				  const double m, const int intervention, 				 
				  const Matrix<>& Y, const Matrix<>& X,
				  Matrix<>& beta, Matrix<>& Sigma, Matrix<>& P, Matrix<int>& s, 
				  Matrix<>& b0, Matrix<>& B0,   
				  const double c0, const double d0,
				  const Matrix<>& A0, 
				  unsigned int burnin, unsigned int mcmc, unsigned int thin,
				  unsigned int verbose, bool chib, bool ar, 
				  Matrix<>& beta_store, Matrix<>& Sigma_store, 
				  Matrix<>& P_store, Matrix<>& ps_store, Matrix<int>& s_store, 
				  double& logmarglike, 
				  Matrix<>& yhat_mat, 
				  Matrix<>& yerror_mat, 
				  Matrix<>& yfore_pred_mat, 
				  Matrix<>& yback_pred_mat)
{
  // define constants and form cross-product matrices
  const int tot_iter = burnin + mcmc; //total iterations
  const int nstore = mcmc / thin; // number of draws to store
  const int n = Y.rows();
  const int ns = m + 1;                 // number of states
  const int k = X.cols();
  const Matrix<> B0inv = invpd(B0);
  Matrix<> sigma(ns, 1);
  
  //MCMC loop
  unsigned int count = 0;  
  unsigned int reject = 0;  
  double rejectionrate = 0; 
  Matrix <> sum_e(ns, 1);

  for (int iter = 0; iter < tot_iter; ++iter){
    //////////////////////
    // 1. Sample beta and Sigma
    //////////////////////
    int beta_count = 0;
    Matrix<int> nstate(ns, 1); 
    Matrix<double> XtX(k, k);   
    Matrix<double> XtY(k, 1);   
    for (int j = 0; j<ns; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == j + 1) { // since j starts from 0
	  nstate[j] = nstate[j] + 1;
	}// end of if
      }// end of int i<n
      beta_count = beta_count + nstate[j];        
 
      Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);

      // SIGMA UPDATE
      double shape = (c0 + (double)nstate[j])/2;
      Matrix<> yhat_j = Xj*beta; 
      Matrix<double> ej = yj - yhat_j;
      Matrix<double> sum_ej = t(ej)*ej;
      double scale =(d0 + sum_ej[0])/2;
      
      Sigma[j] = 1/stream.rgamma(shape, scale);
      sigma[j] = sqrt(Sigma[j]);
      
      if (iter >= burnin && ((iter % thin)==0)){
	yhat_mat(count, (beta_count - nstate[j]), count, (beta_count - 1)) = yhat_j(_,0);
	yerror_mat(count, (beta_count - nstate[j]), count, (beta_count - 1)) = ej(_,0);
      }    
      
      // CONSTANT BETA UPDATE
      XtX = XtX + (crossprod(Xj))/Sigma[j];
      XtY = XtY + ::t(Xj)*yj/Sigma[j];
    }
    Matrix<double> Bn = invpd(B0 + XtX);
    Matrix<double> bn = Bn*(B0*b0 + XtY);
    if (ar == 1){
      Matrix<double> beta_can = stream.rmvnorm(bn, Bn);
      if (beta_can(1) > 1 | beta_can(1) < -1){
	// Rprintf("\n AR coefficient %10.5f is outside the stationary region! \n", beta_can(1));  
	++reject;
      }
      else{
	beta = beta_can;
      }
    }
    else{
      beta = stream.rmvnorm(bn, Bn);
    }
    //////////////////////
    // 2. Sample P
    //////////////////////        
    double shape1 = 0;
    double shape2 = 0;    
    P(ns-1, ns-1) = 1; //no jump at the last state
    
    for (int j =0; j<(ns-1); ++j){
      shape1 =  A0(j,j) + (double)nstate[j] - 1;  
      shape2 =  A0(j,j+1) + 1; // SS(j,j+1);        
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }
    
    //////////////////////
    // 3. Sample s
    //////////////////////
    Matrix<double> F(n, ns);
    Matrix<double> pr1(ns, 1);
    pr1[0] = 1;
    Matrix<double> py(ns, 1);
    Matrix<double> pstyt1(ns, 1);
    Matrix<double> ps = Matrix<double>(n, ns);  // holder for state probabilities
    
    //
    // Forward sampling: update F matrix
    //
    for (int tt=0; tt<n ; ++tt){
      Matrix<double> mu = X(tt,_)*beta; //k by 1 vector
      for (int j = 0; j< ns; ++j){
	py[j]  =  dnorm(Y[tt], mu[0], sigma[j]);
      }
      if (tt==0) pstyt1 = pr1;
      else {
	pstyt1 =  ::t(F(tt-1,_)*P); // make it an ns by 1 matrix
      }
      Matrix<double> unnorm_pstyt = pstyt1%py;

      /////////////////////////////////////////////////////////////////////
      // Prediction of future outcomes based on pre-intervention state
      /////////////////////////////////////////////////////////////////////
      if (tt==(intervention - 1)&&iter >= burnin && ((iter % thin)==0)){
	// Forward prediction
	Matrix <> yfore_pred(1, n);
	for (int ttt=tt; ttt<n ; ++ttt){
	  int ss = s(tt-1); // 1 unit before the intervention
	  mu = X(ttt,_)*beta; //k by 1 vector
	  yfore_pred(ttt) = stream.rnorm(mu(0), sigma[ss-1]);
	}
	yfore_pred_mat(count, _) = yfore_pred(0, _);
	// backward prediction
	Matrix <> yback_pred(1, n);
	for (int ttt=tt; ttt>=0 ; --ttt){
	  int ss = s(tt+1);
	  mu = X(ttt,_)*beta; //k by 1 vector
	  yback_pred(ttt) = stream.rnorm(mu(0), sigma[ss-1]);
	}
	yback_pred_mat(count, _) = yback_pred(0, _);
      }	
    
      const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
      for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);
      
    }// end of F matrix filtering
    
    //
    // Backward recursions: sample s using F and a transition matrix P
    //
    ps(n-1,_) = F(n-1,_);                       // we know last elements of ps and s
    s(n-1) = ns;                                // the last period is s(*n-1) but one down in C location
    
    Matrix<double> pstyn = Matrix<double>(ns, 1);
    double pone = 0.0;
    int tt = n-2;
    while (tt >= 0){
      int st = s(tt+1);
      Matrix<double> Pst_1 = ::t(P(_,st-1)); // prob of being at a previous state
      Matrix<double> unnorm_pstyn = F(tt,_)%Pst_1;
      pstyn = unnorm_pstyn/sum(unnorm_pstyn); // normalize into a prob. density
      if (st==1)   s(tt) = 1;                  // If this is the first period, state should be 1.
      // Otherwise, draw a state from a discrete prob. distribution("pstyn")
      // using the inverse CDF method.
      else{
	pone = pstyn(st-2);
	if(stream.runif() < pone) s(tt) = st-1;// jump from tt-1 to tt 
	else s(tt) = st;// stay 
      }
      ps(tt,_) = pstyn;
      --tt;
    }// end of while loop
    

    // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){
      for (int i=0; i<k; ++i)
	beta_store(count,i) = beta[i];// stored by the order of (11, 12, 13, 21, 22, 23)
      for (int i=0; i<ns; ++i)
	Sigma_store(count, i) = Sigma[i]; 
      for (int j=0; j<ns*ns; ++j)    
	P_store(count,j)= P[j];
      for (int l=0; l<n ; ++l)           
	ps_store(l,_) = ps_store(l,_) + ps(l,_);           // add two matrices
      s_store(count,_) = s;
      
      ++count; // count when (iter % *thin)==0
      
    }   // end of if(iter...)
    
    
    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\nMCMCintervention iteration %i of %i \n", (iter+1), tot_iter);
      if (ar == 1 ){
	double rejectionrate = (double)reject/(double)(iter+1);
	Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	Rprintf("The acceptance rate was %3.5f", 1 - rejectionrate);
	Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       }
      for (int j = 0;j<ns; ++j){
	Rprintf("\n nstate[j] %10.5f", static_cast<double>(nstate[j]));     
      }
      Rprintf("\n beta \n");
      for (int i = 0; i<k; ++i){
	Rprintf("%10.5f\n", beta(i)); 
      }
      Rprintf("\n sigma^2 \n");
      for (int i = 0; i<ns; ++i){
	Rprintf("%10.5f\n", Sigma(i)); 
      }
    }
    
  }// end MCMC loop 
  

  if(chib ==1){\
    double rejectionrate = (double)reject/(double)nstore;    
    if(rejectionrate > .05){
      Rprintf("The acceptance rate is too low.\n"); 
    }
    else {
      Matrix<double> beta_st(k, 1);
      Matrix<double> betast = meanc(beta_store); //meanc(beta_store)=(11, 12, 13, 21, 22, 23) 
      for (int i = 0; i<k ; ++i){
	beta_st(i) = betast(i);
	//Rprintf("beta_st(i) %10.5f\n", beta_st(i)); 
      }
      Matrix<double> Sigma_st = meanc(Sigma_store);
      Matrix<double> sigma_st(ns, 1);
      for (int j = 0; j<ns ; ++j){
	sigma_st(j) = sqrt(Sigma_st(j));
  	//Rprintf("sigma_st(j) %10.5f\n", sigma_st(j)); 
      }
      
      Matrix<double> P_vec_st = meanc(P_store);
      const Matrix<double> P_st(ns, ns);
      for (int j = 0; j< ns*ns; ++j){  
	P_st[j] = P_vec_st[j]; 
      }    
      
      //////////////////////
      // 1. pdf.beta
      //////////////////////      
      Matrix<double> density_beta(nstore, 1);      
      for (int iter = 0; iter<nstore; ++iter){  
	Matrix<double> XtX(k, k);   
	Matrix<double> XtY(k, 1);   
	Matrix<int> nstate(ns, 1); // contains total numbers of each state
	int beta_count = 0;
	
	for (int j = 0; j<ns ; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s_store(iter, i) == (j+1)) { // since j starts from 0
	      nstate[j] = nstate[j] + 1;
	    }// end of if
	  }// end of int i<n
	  beta_count = beta_count + nstate[j];     
	  
	  const Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  const Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	  const double precision = 1.0/Sigma_store(iter, j);
	  XtX = XtX + (crossprod(Xj))*precision;
	  XtY = XtY + ::t(Xj)*yj*precision;
	}
	Matrix<double> Bn = invpd(B0 + XtX);
	Matrix<double> bn = Bn*(B0*b0 + XtY);
	if (k == 1){
	  density_beta(iter) = dnorm(beta_st(0), bn(0), sqrt(Bn(0)));
	}
	else{
	  density_beta(iter) = ::exp(lndmvn(beta_st, bn, Bn));
	}
      } 
      
      double pdf_beta = log(prod(meanc(density_beta)));     
      
      ////////////////////// ////////////////////// //////////////////////
      // 2. pdf.Sigma|beta_st, S, P
      //////////////////////   ////////////////////// //////////////////////    
      Matrix<double> density_Sigma(nstore, ns);
      for (int iter = 0; iter<nstore; ++iter){   
	// STEP 2.1 S|y, beta.st, Sigma, P     
	Matrix <double> Sout = gaussian_state_fixedBeta_sampler(stream, m, Y, X, beta_st, Sigma, P);
	Matrix <double> s = Sout(_, 0);
	
	// STEP 2.2 Sigma|y, beta.st, S, P
	int beta_count = 0;
	Matrix<int> nstate(ns, 1); // contains total numbers of each state
	
	for (int j = 0; j <ns ; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s[i] == (j+1)) { // since j starts from 0
	      nstate[j] = nstate[j] + 1;
	    }// end of if
	  }// end of int i<n
	  beta_count = beta_count + nstate[j];        

	  Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);  
	  Matrix<double> ej = yj - Xj*beta_st;
	  Matrix<double> sum_ej = ::t(ej)*ej;
	  double shape = (c0 + (double)nstate[j])/2;
	  double scale =(d0 + sum_ej[0])/2;
	  
	  Sigma[j] = stream.rigamma(shape, scale);
	  density_Sigma(iter, j) = ::exp(lndinvgamma_pjh(Sigma_st[j], shape, scale));
	}// end of sampling beta and Sigma

	// STEP 2.3 P|S
	double shape1 = 0;
	double shape2 = 0;    
	P(ns-1, ns-1) = 1; //no jump at the last state
	
	for (int j =0; j< (ns-1); ++j){
	  shape1 =  A0(j,j) + (double)nstate[j] - 1;  
	  shape2 =  A0(j,j+1) + 1; //       
	  P(j,j) =  stream.rbeta(shape1, shape2);
	  P(j,j+1) = 1 - P(j,j);
	}  	
      }// end of pdf.Sigma  
      
      double pdf_Sigma = log(prod(meanc(density_Sigma)));   
      
      //////////////////////
      // 3. pdf.P|beta_st, Sigma_st, S
      //////////////////////
      Matrix<double> density_P(nstore, ns);
      for (int iter = 0; iter < nstore; ++iter){
	
	// STEP 2.1 S|y, beta.st, Sigma, P     
	Matrix <double> Sout = gaussian_state_fixedBeta_sampler(stream, m, Y, X, beta_st, Sigma_st, P);
	Matrix <double> s = Sout(_, 0);
	
	double shape1 = 0;
	double shape2 = 0;    
	P(ns-1, ns-1) = 1; //no jump at the last state
	
	// compute addN
	Matrix <double> P_addN(ns, 1); 
	for (int j = 0; j<ns; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s(i) == (j+1)) { // since j starts from 0
	      P_addN[j] = P_addN[j] + 1;            
	    }// end of if
	  } // end of for i
	} // end of for j        
	
	// evaluate P
	for (int j =0; j< (ns-1); ++j){
	  shape1 =  A0(j,j) + P_addN[j] - 1;  
	  shape2 =  A0(j,j+1) + 1; //         
	  P(j,j) = stream.rbeta(shape1, shape2);
	  P(j,j+1) = 1 - P(j,j);
	  density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2); 
	} // end of for j
	density_P(iter, ns-1) = 1; //
	
      }// end of pdf.P MCMC run            
      double pdf_P = log(prod(meanc(density_P)));
      
      //////////////////////
      // likelihood
      //////////////////////       
      Matrix<double> F = Matrix<double>(n, ns);
      Matrix<double> like(n, 1);
      Matrix<double> pr1 = Matrix<double>(ns, 1);
      pr1[0] = 1;
      Matrix<double> py(ns, 1);
      Matrix<double> pstyt1(ns, 1);
      
      for (int t=0; t<n ; ++t){
	int yt = (int) Y[t];
	Matrix<double> mu = X(t,_)*beta_st; 
	for (int j = 0; j< ns; ++j){
	  py[j]  =  dnorm(Y[t], mu[0], sigma_st[j]);
	}
	if (t==0) pstyt1 = pr1;
	else {
	  pstyt1 =  ::t(F(t-1,_)*P_st); 
	}
	Matrix<double> unnorm_pstyt = pstyt1%py;
	Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
	for (int j=0; j<ns ; ++j){
	  F(t,j) = pstyt(j);
	}
	like[t] = sum(unnorm_pstyt);
      }// end of likelihood computation
      
      const double loglike = sum(log(like));
      
      //////////////////////
      // log prior ordinate 
      //////////////////////
      double density_beta_prior = 0.0;
      Matrix<double> density_Sigma_prior(ns, 1);
      Matrix<double> density_P_prior(ns, 1);
      density_P[ns-1] = 1; //
      if (k == 1){
	density_beta_prior = log(dnorm(beta_st(0), b0(0), sqrt(B0inv(0)))); 
      }
      else{
	density_beta_prior = lndmvn(beta_st, b0, B0inv); 
      }
      for (int j=0; j<ns ; ++j){
	density_Sigma_prior[j] = lndinvgamma_pjh(Sigma_st[j], c0/2, d0/2); 
      }   
      
      for (int j =0; j< (ns-1); ++j){
	density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1))); 
      }        
      
      // compute marginal likelihood
      double logprior = density_beta_prior + sum(density_Sigma_prior) + sum(density_P_prior);
      logmarglike = (loglike + logprior) - (pdf_beta + pdf_Sigma + pdf_P);
      if(verbose > 0){
	Rprintf("logmarglike = %10.5f\n", logmarglike);
	Rprintf("loglike = %10.5f\n", loglike);
	Rprintf("logprior = %10.5f\n", logprior);
	Rprintf("pdf_beta = %10.5f\n", pdf_beta);
	Rprintf("pdf_Sigma = %10.5f\n", pdf_Sigma);
	Rprintf("pdf_P = %10.5f\n", pdf_P);
      }
    }
  }// end of marginal likelihood computation
}


//////////////////////////////////////////// 
// MCMCintervention fixed effect changes only  
//////////////////////////////////////////// 
template <typename RNGTYPE>
void MCMCintervention_fixed_impl(rng<RNGTYPE>& stream, 
				 const double m, const int intervention, 				 
				 const Matrix<>& Y, const Matrix<>& X,
				 Matrix<>& beta, Matrix<>& Sigma, Matrix<>& P, Matrix<int>& s, 
				 Matrix<>& b0, Matrix<>& B0,   
				 const double c0, const double d0,
				 const Matrix<>& A0, 
				 unsigned int burnin, unsigned int mcmc, unsigned int thin,
				 unsigned int verbose, bool chib, bool ar, 
				 Matrix<>& beta_store, Matrix<>& Sigma_store, 
				 Matrix<>& P_store, Matrix<>& ps_store, Matrix<int>& s_store, 
				 double& logmarglike, 
				 Matrix<>& yhat_mat, 
				 Matrix<>& yerror_mat, 
				 Matrix<>& yfore_pred_mat, 
				 Matrix<>& yback_pred_mat, 
				 double acceptance)
{
  // define constants and form cross-product matrices
  const int tot_iter = burnin + mcmc; //total iterations
  const int nstore = mcmc / thin; // number of draws to store
  const int n = Y.rows();
  const int ns = m + 1;                 // number of states
  const int k = X.cols();
  const Matrix<> B0inv = invpd(B0);
  double sigma2 = Sigma(0);
  double sigma = sqrt(sigma2);
  
  //MCMC loop
  unsigned int count = 0;  
  unsigned int reject = 0;  
  Matrix <> sum_e(ns, 1);

  for (int iter = 0; iter < tot_iter; ++iter){
    //////////////////////
    // 1. Sample beta and Sigma
    //////////////////////
    int beta_count = 0;
    Matrix<int> nstate(ns, 1); // contains total numbers of each state    
    for (int j = 0; j<ns; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == j + 1) { // since j starts from 0
	  nstate[j] = nstate[j] + 1;
	}// end of if
      }// end of int i<n
      beta_count = beta_count + nstate[j];        
 
      // BETA UPDATE
      Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
      Matrix<double> tXj = ::t(Xj);
      Matrix<double> Bn = invpd(B0 + crossprod(Xj)/sigma2);
      Matrix<double> bn = Bn*(B0*b0 + tXj*yj/sigma2);
      if (ar == 1){
	Matrix<double> beta_can = stream.rmvnorm(bn, Bn);
	if (beta_can(1) > 1 | beta_can(1) < -1){
	  // Rprintf("\n AR coefficient %10.5f is outside the stationary region! \n", beta_can(1));  
	  ++reject;
	}
	else{
	  for (int kk = 0; kk<k; ++kk){
	    beta(j,kk) = beta_can(kk);
	  }
	}
      }
      else{
	beta(j,_) = stream.rmvnorm(bn, Bn);
      }
      Matrix<> yhat_j = Xj*::t(beta(j,_)); 
      Matrix<double> ej = yj - yhat_j;
      Matrix<double> sum_ej = t(ej)*ej;
      sum_e(j) = sum_ej[0];
      if (iter >= burnin && ((iter % thin)==0)){
	yhat_mat(count, (beta_count - nstate[j]), count, (beta_count - 1)) = yhat_j(_,0);
	yerror_mat(count, (beta_count - nstate[j]), count, (beta_count - 1)) = ej(_,0);
      }
      
    }// end of sampling beta and Sigma
    
    // SIGMA UPDATE
    double shape = (c0 + (double)n)/2;
    double scale =(d0 + sum(sum_e))/2;
    
    sigma2 = 1/stream.rgamma(shape, scale);
    sigma = sqrt(sigma2);
     
    //////////////////////
    // 2. Sample P
    //////////////////////        
    double shape1 = 0;
    double shape2 = 0;    
    P(ns-1, ns-1) = 1; //no jump at the last state
    
    for (int j =0; j<(ns-1); ++j){
      shape1 =  A0(j,j) + (double)nstate[j] - 1;  
      shape2 =  A0(j,j+1) + 1; // SS(j,j+1);        
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }
    
    //////////////////////
    // 3. Sample s
    //////////////////////
    Matrix<double> F(n, ns);
    Matrix<double> pr1(ns, 1);
    pr1[0] = 1;
    Matrix<double> py(ns, 1);
    Matrix<double> pstyt1(ns, 1);
    Matrix<double> ps = Matrix<double>(n, ns);  // holder for state probabilities
    
    //
    // Forward sampling: update F matrix
    //
    for (int tt=0; tt<n ; ++tt){
      Matrix<double> mu = X(tt,_)*::t(beta); //k by 1 vector
      for (int j = 0; j< ns; ++j){
	py[j]  =  dnorm(Y[tt], mu[j], sigma);
      }
      if (tt==0) pstyt1 = pr1;
      else {
	pstyt1 =  ::t(F(tt-1,_)*P); // make it an ns by 1 matrix
      }
      Matrix<double> unnorm_pstyt = pstyt1%py;

      /////////////////////////////////////////////////////////////////////
      // Prediction of future outcomes based on pre-intervention state
      /////////////////////////////////////////////////////////////////////
      if (tt==(intervention - 1)&&iter >= burnin && ((iter % thin)==0)){
	// Forward prediction
	Matrix <> yfore_pred(1, n);
	for (int ttt=tt; ttt<n ; ++ttt){
	  int ss = s(tt-1); // 1 unit before the intervention
	  mu = X(ttt,_)*::t(beta); //k by 1 vector
	  yfore_pred(ttt) = stream.rnorm(mu[ss-1], sigma);
	}
	yfore_pred_mat(count, _) = yfore_pred(0, _);
	// backward prediction
	Matrix <> yback_pred(1, n);
	for (int ttt=tt; ttt>=0 ; --ttt){
	  int ss = s(tt+1);
	  mu = X(ttt,_)*::t(beta); //k by 1 vector
	  yback_pred(ttt) = stream.rnorm(mu[ss-1], sigma);
	}
	yback_pred_mat(count, _) = yback_pred(0, _);
      }	
    
      const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
      for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);
      
    }// end of F matrix filtering
    
    //
    // Backward recursions: sample s using F and a transition matrix P
    //
    ps(n-1,_) = F(n-1,_);                       // we know last elements of ps and s
    s(n-1) = ns;                                // the last period is s(*n-1) but one down in C location
    
    Matrix<double> pstyn = Matrix<double>(ns, 1);
    double pone = 0.0;
    int tt = n-2;
    while (tt >= 0){
      int st = s(tt+1);
      Matrix<double> Pst_1 = ::t(P(_,st-1)); // prob of being at a previous state
      Matrix<double> unnorm_pstyn = F(tt,_)%Pst_1;
      pstyn = unnorm_pstyn/sum(unnorm_pstyn); // normalize into a prob. density
      if (st==1)   s(tt) = 1;                  // If this is the first period, state should be 1.
      // Otherwise, draw a state from a discrete prob. distribution("pstyn")
      // using the inverse CDF method.
      else{
	pone = pstyn(st-2);
	if(stream.runif() < pone) s(tt) = st-1;// jump from tt-1 to tt 
	else s(tt) = st;// stay 
      }
      ps(tt,_) = pstyn;
      --tt;
    }// end of while loop
    

    // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){
      Matrix<double> tbeta = ::t(beta); //transpose beta for R output 
      for (int i=0; i<(ns*k); ++i)
	beta_store(count,i) = tbeta[i];// stored by the order of (11, 12, 13, 21, 22, 23)
      Sigma_store(count) = sigma2; 
      for (int j=0; j<ns*ns; ++j)    
	P_store(count,j)= P[j];
      for (int l=0; l<n ; ++l)           
	ps_store(l,_) = ps_store(l,_) + ps(l,_);           // add two matrices
      s_store(count,_) = s;
      
      ++count; // count when (iter % *thin)==0
      
    }   // end of if(iter...)
    
    
    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\nMCMCintervention iteration %i of %i \n", (iter+1), tot_iter);
      if (ar == 1 ){
	double rejectionrate = (double)reject/(double)(iter+1);
	Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	Rprintf("The acceptance rate was %3.5f", 1 - rejectionrate);
	Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      }
      for (int j = 0;j<ns; ++j){
	Rprintf("\n nstate[j] %10.5f", static_cast<double>(nstate[j]));     
      }
      Rprintf("\n beta \n");
      for (int i = 0; i<ns; ++i){
	for (int j = 0; j<k; ++j){
	  Rprintf("%10.5f\n", beta(i, j)); 
	}
      }
      Rprintf("\n sigma^2 \n");
      Rprintf("%10.5f\n", sigma2); 
    }
    
  }// end MCMC loop 
  
  acceptance = 1 - (double)reject/(double)nstore;    
  
  if(chib ==1){
    // 0.05 is used for an arbitrary threshold to avoid the computation of marg like
    // Later, the threshold should be set by users
    if(acceptance < .95){
      Rprintf("The acceptance rate is too low.\n"); 
    }
    else{
      Matrix<double> betast = meanc(beta_store); //meanc(beta_store)=(11, 12, 13, 21, 22, 23) 
      Matrix<double, Row> beta_st(ns, k);
      for (int j = 0; j<ns*k; ++j){
	beta_st[j] = betast[j];
      }
      
      Matrix<double> Sigma_st(ns, 1);
      for (int j = 0; j<ns ; ++j){
	Sigma_st(j) = mean(Sigma_store);
      }
      Matrix<double> P_vec_st = meanc(P_store);
      const Matrix<double> P_st(ns, ns);
      for (int j = 0; j< ns*ns; ++j){  
	P_st[j] = P_vec_st[j]; 
      }    
      
      //////////////////////
      // 1. pdf.beta
      //////////////////////      
      Matrix<double> density_beta(nstore, ns);      
      for (int iter = 0; iter<nstore; ++iter){  
	
	Matrix<int> nstate(ns, 1); // contains total numbers of each state
	int beta_count = 0;
	
	for (int j = 0; j<ns ; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s_store(iter, i) == (j+1)) { // since j starts from 0
	      nstate[j] = nstate[j] + 1;
	    }// end of if
	  }// end of int i<n
	  beta_count = beta_count + nstate[j];     
	  
	  // BETA UPDATE
	  const Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  const Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	  const double precision = 1.0/Sigma_store(iter);
	  const Matrix<double> XpX = crossprod(Xj);
	  const Matrix<double> XpY = (::t(Xj)*yj);
	  const Matrix<double> Bn = invpd(B0 + XpX*precision);
	  const Matrix<double> bn = Bn*gaxpy(B0, b0, XpY*precision);
	  if (k == 1){
	    density_beta(iter, j) = dnorm(beta_st(j), bn(0), sqrt(Bn(0)));
	  }
	  else{
	    density_beta(iter, j) = ::exp(lndmvn(::t(beta_st(j,_)), bn, Bn));
	  }
	}// end of sampling beta and Sigma
	
      }// end of pdf.beta   
      
      double pdf_beta = log(prod(meanc(density_beta)));     
    
      //////////////////////
      // 2. pdf.Sigma|beta_st, S, P
      //////////////////////      
      Matrix<double> density_Sigma(nstore, 1);
      for (int iter = 0; iter<nstore; ++iter){   
	
	// STEP 2.1 S|y, beta.st, Sigma, P     
	Matrix <double> Sout = gaussian_state_sampler(stream, m, Y, X, beta_st, Sigma, P);
	Matrix <double> s = Sout(_, 0);
	
	// STEP 2.2 Sigma|y, beta.st, S, P
	int beta_count = 0;
	Matrix<int> nstate(ns, 1); 
	Matrix<> sum_e(ns, 1);
	
	for (int j = 0; j <ns ; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s[i] == (j+1)) { // since j starts from 0
	      nstate[j] = nstate[j] + 1;
	    }// end of if
	  }// end of int i<n
	  beta_count = beta_count + nstate[j];        
	  
	  Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);  
	  Matrix<double> ej = yj - Xj*::t(beta_st(j,_));
	  Matrix<double> sum_ej = ::t(ej)*ej;
	  sum_e(j) = sum_ej(0);
	}// end of sampling beta and Sigma
	double shape = (c0 + (double)n)/2;
	double scale =(d0 + sum(sum_e))/2;
	
	sigma2 = stream.rigamma(shape, scale);
	density_Sigma(iter) = ::exp(lndinvgamma_pjh(Sigma_st(0), shape, scale));

	for (int j = 0; j <ns ; ++j){
	  Sigma[j] = sigma2;// to use gaussian_state_sampler
	}

	// STEP 2.3 P|S
	double shape1 = 0;
	double shape2 = 0;    
	P(ns-1, ns-1) = 1; //no jump at the last state
	
	for (int j =0; j< (ns-1); ++j){
	  shape1 =  A0(j,j) + (double)nstate[j] - 1;  
	  shape2 =  A0(j,j+1) + 1; //       
	  P(j,j) =  stream.rbeta(shape1, shape2);
	  P(j,j+1) = 1 - P(j,j);

	  // not to change gaussian_state_sampler, I use Sigma matrix
	  Sigma[j] = sigma2;
	}  
	
      }// end of pdf.Sigma  
      
      double pdf_Sigma = log(prod(meanc(density_Sigma)));   
      
      //////////////////////
      // 3. pdf.P|beta_st, Sigma_st, S
      //////////////////////
      Matrix<double> density_P(nstore, ns);
      for (int iter = 0; iter < nstore; ++iter){
	
	// STEP 2.1 S|y, beta.st, Sigma, P     
	Matrix <double> Sout = gaussian_state_sampler(stream, m, Y, X, beta_st, Sigma_st, P);
	Matrix <double> s = Sout(_, 0);
	
	double shape1 = 0;
	double shape2 = 0;    
	P(ns-1, ns-1) = 1; //no jump at the last state
	
	// compute addN
	Matrix <double> P_addN(ns, 1); 
	for (int j = 0; j<ns; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s(i) == (j+1)) { // since j starts from 0
	      P_addN[j] = P_addN[j] + 1;            
	    }// end of if
	  } // end of for i
	} // end of for j        
	
	// evaluate P
	for (int j =0; j< (ns-1); ++j){
	  shape1 =  A0(j,j) + P_addN[j] - 1;  
	  shape2 =  A0(j,j+1) + 1; //         
	  P(j,j) = stream.rbeta(shape1, shape2);
	  P(j,j+1) = 1 - P(j,j);
	  density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2); 
	} // end of for j
	density_P(iter, ns-1) = 1; //
	
      }// end of pdf.P MCMC run            
      double pdf_P = log(prod(meanc(density_P)));
      
      //////////////////////
      // likelihood
      //////////////////////       
      Matrix<double> F = Matrix<double>(n, ns);
      Matrix<double> like(n, 1);
      Matrix<double> pr1 = Matrix<double>(ns, 1);
      pr1[0] = 1;
      Matrix<double> py(ns, 1);
      Matrix<double> pstyt1(ns, 1);
      
      for (int t=0; t<n ; ++t){
	int yt = (int) Y[t];
	Matrix<double> mu = X(t,_)*::t(beta_st); //k by 1 vector
	for (int j = 0; j< ns; ++j){
	  py[j]  =  dnorm(Y[t], mu[j], sqrt(Sigma_st[0]));
	}
	if (t==0) pstyt1 = pr1;
	else {
	  pstyt1 =  ::t(F(t-1,_)*P_st); // make it an ns by 1 matrix
	}
	Matrix<double> unnorm_pstyt = pstyt1%py;
	Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
	for (int j=0; j<ns ; ++j){
	  F(t,j) = pstyt(j);
	}
	like[t] = sum(unnorm_pstyt);
      }// end of likelihood computation
      
      const double loglike = sum(log(like));
      
      //////////////////////
      // log prior ordinate 
      //////////////////////
      Matrix<double> density_beta_prior(ns, 1);
      double density_Sigma_prior = 0.0;
      Matrix<double> density_P_prior(ns, 1);
      density_P[ns-1] = 1; //
      if (k == 1){
	for (int j=0; j<ns ; ++j){
	  density_beta_prior[j] = log(dnorm(beta_st(j), b0(0), sqrt(B0inv(0)))); 
	}   
      }
      else{
	for (int j=0; j<ns ; ++j){
	  density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv); 
	}   
      }
      
      density_Sigma_prior = lndinvgamma_pjh(Sigma_st(0), c0/2, d0/2); 
      
      for (int j =0; j< (ns-1); ++j){
	density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1))); 
      }        
      
      // compute marginal likelihood
      double logprior = sum(density_beta_prior) + density_Sigma_prior + sum(density_P_prior);
      logmarglike = (loglike + logprior) - (pdf_beta + pdf_Sigma + pdf_P);
      if(verbose > 0){
	Rprintf("logmarglike = %10.5f\n", logmarglike);
	Rprintf("loglike = %10.5f\n", loglike);
	Rprintf("logprior = %10.5f\n", logprior);
	Rprintf("pdf_beta = %10.5f\n", pdf_beta);
	Rprintf("pdf_Sigma = %10.5f\n", pdf_Sigma);
	Rprintf("pdf_P = %10.5f\n", pdf_P);
      }
    }
  }// end of marginal likelihood computation
}

template <typename RNGTYPE>
void MCMCintervention_impl(rng<RNGTYPE>& stream, 
			   const double m, const int intervention, 				 
			   const Matrix<>& Y, const Matrix<>& X,
			   Matrix<>& beta, Matrix<>& Sigma, Matrix<>& P, Matrix<int>& s, 
			   Matrix<>& b0, Matrix<>& B0,   
			   const double c0, const double d0,
			   const Matrix<>& A0, 
			   unsigned int burnin, unsigned int mcmc, unsigned int thin,
			   unsigned int verbose, bool chib, bool ar, 
			   Matrix<>& beta_store, Matrix<>& Sigma_store, 
			   Matrix<>& P_store, Matrix<>& ps_store, Matrix<int>& s_store, 
			   double& logmarglike, 
			   Matrix<>& yhat_mat, 
			   Matrix<>& yerror_mat, 
			   Matrix<>& yfore_pred_mat, 
			   Matrix<>& yback_pred_mat, 
			   double acceptance)
{
  // define constants and form cross-product matrices
  const int tot_iter = burnin + mcmc; //total iterations
  const int nstore = mcmc / thin; // number of draws to store
  const int n = Y.rows();
  const int ns = m + 1;                 // number of states
  const int k = X.cols();
  const Matrix<> B0inv = invpd(B0);
  Matrix<double> sigma(ns, 1);

  //MCMC loop
  unsigned int count = 0;  
  int reject = 0;  
  for (int iter = 0; iter < tot_iter; ++iter){
    //////////////////////
    // 1. Sample beta and Sigma
    //////////////////////
    int beta_count = 0;
    Matrix<int> nstate(ns, 1); // contains total numbers of each state    
    for (int j = 0; j<ns; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == j + 1) { // since j starts from 0
	  nstate[j] = nstate[j] + 1;
	}// end of if
      }// end of int i<n
      beta_count = beta_count + nstate[j];        
 
      // BETA UPDATE
      Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
      Matrix<double> tXj = ::t(Xj);
      Matrix<double> Bn = invpd(B0 + crossprod(Xj)/Sigma[j]);
      Matrix<double> bn = Bn*(B0*b0 + tXj*yj/Sigma[j]);
      if (ar == 1){
	Matrix<double> beta_can = stream.rmvnorm(bn, Bn);
	if (beta_can(1) > 1 | beta_can(1) < -1){
	  // Rprintf("\n AR coefficient %10.5f is outside the stationary region! \n", beta_can(1));  
	  ++reject;
	}
	else{
	  for (int kk = 0; kk<k; ++kk){
	    beta(j,kk) = beta_can(kk);
	  }
	}
      }
      else{
	beta(j,_) = stream.rmvnorm(bn, Bn);
      }
      
      // SIGMA UPDATE
      double shape = (c0 + (double)nstate[j])/2;
      Matrix<> yhat_j = Xj*::t(beta(j,_)); 
      Matrix<double> ej = yj - yhat_j;
      Matrix<double> sum_ej = t(ej)*ej;
      double scale =(d0 + sum_ej[0])/2;
      
      Sigma[j] = 1/stream.rgamma(shape, scale);
      sigma(j) = sqrt(Sigma(j));
      

      if (iter >= burnin && ((iter % thin)==0)){
	yhat_mat(count, (beta_count - nstate[j]), count, (beta_count - 1)) = yhat_j(_,0);
	yerror_mat(count, (beta_count - nstate[j]), count, (beta_count - 1)) = ej(_,0);
      }
      
    }// end of sampling beta and Sigma
    
    //////////////////////
    // 2. Sample P
    //////////////////////        
    double shape1 = 0;
    double shape2 = 0;    
    P(ns-1, ns-1) = 1; //no jump at the last state
    
    for (int j =0; j<(ns-1); ++j){
      shape1 =  A0(j,j) + (double)nstate[j] - 1;  
      shape2 =  A0(j,j+1) + 1; // SS(j,j+1);        
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }
    
    //////////////////////
    // 3. Sample s
    //////////////////////
    Matrix<double> F(n, ns);
    Matrix<double> pr1(ns, 1);
    pr1[0] = 1;
    Matrix<double> py(ns, 1);
    Matrix<double> pstyt1(ns, 1);
    Matrix<double> ps = Matrix<double>(n, ns);  // holder for state probabilities
    
    //
    // Forward sampling: update F matrix
    //
    for (int tt=0; tt<n ; ++tt){
      Matrix<double> mu = X(tt,_)*::t(beta); //k by 1 vector
      for (int j = 0; j< ns; ++j){
	py[j]  =  dnorm(Y[tt], mu[j], sigma[j]);
      }
      if (tt==0) pstyt1 = pr1;
      else {
	pstyt1 =  ::t(F(tt-1,_)*P); // make it an ns by 1 matrix
      }
      Matrix<double> unnorm_pstyt = pstyt1%py;

      /////////////////////////////////////////////////////////////////////
      // Prediction of future outcomes based on pre-intervention state
      /////////////////////////////////////////////////////////////////////
      if (tt==(intervention - 1)&&iter >= burnin && ((iter % thin)==0)){
	// Forward prediction
	Matrix <> yfore_pred(1, n);
	for (int ttt=tt; ttt<n ; ++ttt){
	  int ss = s(tt-1); // 1 unit before the intervention
	  mu = X(ttt,_)*::t(beta); //k by 1 vector
	  yfore_pred(ttt) = stream.rnorm(mu[ss-1], sigma[ss-1]);
	}
	yfore_pred_mat(count, _) = yfore_pred(0, _);
	// backward prediction
	Matrix <> yback_pred(1, n);
	for (int ttt=tt; ttt>=0 ; --ttt){
	  int ss = s(tt+1);
	  mu = X(ttt,_)*::t(beta); //k by 1 vector
	  yback_pred(ttt) = stream.rnorm(mu[ss-1], sigma[ss-1]);
	}
	yback_pred_mat(count, _) = yback_pred(0, _);
      }	
    
      const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
      for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);
      
    }// end of F matrix filtering
    
    //
    // Backward recursions: sample s using F and a transition matrix P
    //
    ps(n-1,_) = F(n-1,_);                       // we know last elements of ps and s
    s(n-1) = ns;                                // the last period is s(*n-1) but one down in C location
    
    Matrix<double> pstyn = Matrix<double>(ns, 1);
    double pone = 0.0;
    int tt = n-2;
    while (tt >= 0){
      int st = s(tt+1);
      Matrix<double> Pst_1 = ::t(P(_,st-1)); // prob of being at a previous state
      Matrix<double> unnorm_pstyn = F(tt,_)%Pst_1;
      pstyn = unnorm_pstyn/sum(unnorm_pstyn); // normalize into a prob. density
      if (st==1)   s(tt) = 1;                  // If this is the first period, state should be 1.
      // Otherwise, draw a state from a discrete prob. distribution("pstyn")
      // using the inverse CDF method.
      else{
	pone = pstyn(st-2);
	if(stream.runif() < pone) s(tt) = st-1;// jump from tt-1 to tt 
	else s(tt) = st;// stay 
      }
      ps(tt,_) = pstyn;
      --tt;
    }// end of while loop
    

    // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){
      Matrix<double> tbeta = ::t(beta); //transpose beta for R output 
      for (int i=0; i<(ns*k); ++i)
	beta_store(count,i) = tbeta[i];// stored by the order of (11, 12, 13, 21, 22, 23)
      for (int i=0; i<ns; ++i)
	Sigma_store(count,i) = Sigma[i]; 
      for (int j=0; j<ns*ns; ++j)    
	P_store(count,j)= P[j];
      for (int l=0; l<n ; ++l)           
	ps_store(l,_) = ps_store(l,_) + ps(l,_);           // add two matrices
      s_store(count,_) = s;
      
      ++count; // count when (iter % *thin)==0
      
    }   // end of if(iter...)
    
    
    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\nMCMCintervention iteration %i of %i \n", (iter+1), tot_iter);
      if (ar == 1 ){
	double rejectionrate = (double)reject/(double)(iter+1);
	Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	Rprintf("The acceptance rate was %3.5f", 1 - rejectionrate);
	Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      }
      for (int j = 0;j<ns; ++j){
	Rprintf("\n nstate[j] %10.5f", static_cast<double>(nstate[j]));     
      }
      Rprintf("\n beta \n");
      for (int i = 0; i<ns; ++i){
	for (int j = 0; j<k; ++j){
	  Rprintf("%10.5f\n", beta(i, j)); 
	}
      }
      Rprintf("\n sigma^2 \n");
      for (int i = 0; i<ns; ++i){
	Rprintf("%10.5f\n", Sigma(i)); 
      }
    }
    
  }// end MCMC loop 
 
  acceptance = 1 - (double)reject/(double)nstore;    
  
  if(chib ==1){
    // 0.05 is used for an arbitrary threshold to avoid the computation of marg like
    // Later, the threshold should be set by users
    if(acceptance < .95){
      Rprintf("The acceptance rate is too low.\n"); 
    }
    else{
      Matrix<double> betast = meanc(beta_store); 
      Matrix<double, Row> beta_st(ns, k);
      for (int j = 0; j<ns*k; ++j){
	beta_st[j] = betast[j];
      }
      
      Matrix<double> Sigma_st = meanc(Sigma_store);
      Matrix<double> sigma_st(ns, 1);
      for (int j = 0; j<ns ; ++j){
	sigma_st(j) = sqrt(Sigma_st(j));
      }
      
      Matrix<double> P_vec_st = meanc(P_store);
      const Matrix<double> P_st(ns, ns);
      for (int j = 0; j< ns*ns; ++j){  
	P_st[j] = P_vec_st[j]; 
      }    
      
      //////////////////////
      // 1. pdf.beta
      //////////////////////      
      Matrix<double> density_beta(nstore, ns);      
      for (int iter = 0; iter<nstore; ++iter){  
	
	Matrix<int> nstate(ns, 1); 
	int beta_count = 0;
	
	for (int j = 0; j<ns ; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s_store(iter, i) == (j+1)) { 
	      nstate[j] = nstate[j] + 1;
	    }// end of if
	  }// end of int i<n
	  beta_count = beta_count + nstate[j];     
	  
	  // BETA UPDATE
	  const Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  const Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	  const double precision = 1.0/Sigma_store(iter, j);
	  const Matrix<double> XpX = crossprod(Xj);
	  const Matrix<double> XpY = (::t(Xj)*yj);
	  const Matrix<double> Bn = invpd(B0 + XpX*precision);
	  const Matrix<double> bn = Bn*gaxpy(B0, b0, XpY*precision);
	  if (k == 1){
	    density_beta(iter, j) = log(dnorm(beta_st(j), bn(0), sqrt(Bn(0))));
	  }
	  else{
	    density_beta(iter, j) = lndmvn(::t(beta_st(j,_)), bn, Bn);
	  }
	}// end of sampling beta and Sigma
	
      }// end of pdf.beta   
      
      double pdf_beta = sum(meanc(density_beta));     
    
      //////////////////////
      // 2. pdf.Sigma|beta_st, S, P
      //////////////////////      
      Matrix<double> density_Sigma(nstore, ns);
      for (int iter = 0; iter<nstore; ++iter){   
	
	// STEP 2.1 S|y, beta.st, Sigma, P     
	Matrix <double> Sout = gaussian_state_sampler(stream, m, Y, X, beta_st, Sigma, P);
	Matrix <double> s = Sout(_, 0);
	
	// STEP 2.2 Sigma|y, beta.st, S, P
	int beta_count = 0;
	Matrix<int> nstate(ns, 1); // contains total numbers of each state
	
	for (int j = 0; j <ns ; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s[i] == (j+1)) { // since j starts from 0
	      nstate[j] = nstate[j] + 1;
	    }// end of if
	  }// end of int i<n
	  beta_count = beta_count + nstate[j];        
	  
	  Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);  
	  Matrix<double> ej = yj - Xj*::t(beta_st(j,_));
	  Matrix<double> sum_ej = ::t(ej)*ej;
	  double shape = (c0 + (double)nstate[j])/2;
	  double scale =(d0 + sum_ej[0])/2;
	  
	  Sigma[j] = stream.rigamma(shape, scale);
	  density_Sigma(iter, j) = ::exp(lndinvgamma_pjh(Sigma_st[j], shape, scale));
	  
	}// end of sampling beta and Sigma
	
	// STEP 2.3 P|S
	double shape1 = 0;
	double shape2 = 0;    
	P(ns-1, ns-1) = 1; //no jump at the last state
	
	for (int j =0; j< (ns-1); ++j){
	  shape1 =  A0(j,j) + (double)nstate[j] - 1;  
	  shape2 =  A0(j,j+1) + 1; //       
	  P(j,j) =  stream.rbeta(shape1, shape2);
	  P(j,j+1) = 1 - P(j,j);
	}  
	
      }// end of pdf.Sigma  
      
      double pdf_Sigma = log(prod(meanc(density_Sigma)));   
      
      //////////////////////
      // 3. pdf.P|beta_st, Sigma_st, S
      //////////////////////
      Matrix<double> density_P(nstore, ns);
      for (int iter = 0; iter < nstore; ++iter){
	
	// STEP 2.1 S|y, beta.st, Sigma, P     
	Matrix <double> Sout = gaussian_state_sampler(stream, m, Y, X, beta_st, Sigma_st, P);
	Matrix <double> s = Sout(_, 0);
	
	double shape1 = 0;
	double shape2 = 0;    
	P(ns-1, ns-1) = 1;
	
	// compute addN
	Matrix <double> P_addN(ns, 1); 
	for (int j = 0; j<ns; ++j){
	  for (int i = 0; i<n; ++i){
	    if (s(i) == (j+1)) { 
	      P_addN[j] = P_addN[j] + 1;            
	    }// end of if
	  } // end of for i
	} // end of for j        
	
	// evaluate P
	for (int j =0; j< (ns-1); ++j){
	  shape1 =  A0(j,j) + P_addN[j] - 1;  
	  shape2 =  A0(j,j+1) + 1; //         
	  P(j,j) = stream.rbeta(shape1, shape2);
	  P(j,j+1) = 1 - P(j,j);
	  density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2); 
	} // end of for j
	density_P(iter, ns-1) = 1; //
	
      }// end of pdf.P MCMC run            
      double pdf_P = log(prod(meanc(density_P)));
      
      //////////////////////
      // likelihood
      //////////////////////       
      Matrix<double> F = Matrix<double>(n, ns);
      Matrix<double> like(n, 1);
      Matrix<double> pr1 = Matrix<double>(ns, 1);
      pr1[0] = 1;
      Matrix<double> py(ns, 1);
      Matrix<double> pstyt1(ns, 1);
      
      for (int t=0; t<n ; ++t){
	int yt = (int) Y[t];
	Matrix<double> mu = X(t,_)*::t(beta_st); 
	for (int j = 0; j< ns; ++j){
	  py[j]  =  dnorm(Y[t], mu[j], sigma_st[j]);
	}
	if (t==0) pstyt1 = pr1;
	else {
	  pstyt1 =  ::t(F(t-1,_)*P_st); 
	}
	Matrix<double> unnorm_pstyt = pstyt1%py;
	Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); 
	for (int j=0; j<ns ; ++j){
	  F(t,j) = pstyt(j);
	}
	like[t] = sum(unnorm_pstyt);
      }// end of likelihood computation
      
      const double loglike = sum(log(like));
      
      //////////////////////
      // log prior ordinate 
      //////////////////////
      Matrix<double> density_beta_prior(ns, 1);
      Matrix<double> density_Sigma_prior(ns, 1);
      Matrix<double> density_P_prior(ns, 1);
      density_P[ns-1] = 1; //
      if (k == 1){
	for (int j=0; j<ns ; ++j){
	  density_beta_prior[j] = log(dnorm(beta_st(j), b0(0), sqrt(B0inv(0)))); 
	}   
      }
      else{
	for (int j=0; j<ns ; ++j){
	  density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv); 
	}   
      }
      
      for (int j=0; j<ns ; ++j){
	density_Sigma_prior[j] = lndinvgamma_pjh(Sigma_st[j], c0/2, d0/2); 
      }   
      
      for (int j =0; j< (ns-1); ++j){
	density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1))); 
      }        
      
      // compute marginal likelihood
      double logprior = sum(density_beta_prior) + sum(density_Sigma_prior) + sum(density_P_prior);
      logmarglike = (loglike + logprior) - (pdf_beta + pdf_Sigma + pdf_P);
      if(verbose > 0){
	Rprintf("logmarglike = %10.5f\n", logmarglike);
	Rprintf("loglike = %10.5f\n", loglike);
	Rprintf("logprior = %10.5f\n", logprior);
	Rprintf("pdf_beta = %10.5f\n", pdf_beta);
	Rprintf("pdf_Sigma = %10.5f\n", pdf_Sigma);
	Rprintf("pdf_P = %10.5f\n", pdf_P);
      }
    }// end of marginal likelihood computation
  }
}

//////////////////////////////////////////// 
// Start MCMCinterventionpoint function
///////////////////////////////////////////
extern "C"{
  void MCMCintervention(double *accept, double *betaout, double *Sigmaout, double *Pout, 
			double *psout, double *sout,  
			double *yhatout, double *yerrorout, double *yforepredout, double *ybackpredout,
			const double *Ydata, const int *Yrow, const int *Ycol, 
			const double *Xdata, const int *Xrow, const int *Xcol, 
			const int *m, const int *intervention, 
			const int *burnin, const int *mcmc, const int *thin, const int *verbose, 
			const int *uselecuyer, const int *seedarray, const int *lecuyerstream, 
			const double *betastart, const double *Sigmastart, const double *Pstart, 
			const int *statestart, 
			const double *a, const double *b, 
			const double *b0data, const double *B0data, 
			const double *c0, const double *d0, 
			const double *A0data, 
			double *logmarglikeholder, double *loglikeholder, 
			const int *ar, const int *change, const int *chib){
    
    // pull together Matrix objects
    const Matrix <double> Y(*Yrow, *Ycol, Ydata); 
    const Matrix <double> X(*Xrow, *Xcol, Xdata); 
    const unsigned int tot_iter = *burnin + *mcmc; //total iterations
    const unsigned int nstore = *mcmc / *thin; // number of draws to store
    const int n = Y.rows();
    const int k = X.cols();
    const int ns = *m + 1;                 // number of states
    
    // generate starting values
    Matrix <> Sigma(ns, 1, Sigmastart);
    Matrix <> P(ns, ns, Pstart);    
    Matrix <int> s(n, 1, statestart); 
    Matrix <> b0(k, 1, b0data);
    Matrix <> B0(k, k, B0data);
    const Matrix <> A0(ns, ns, A0data);
    double logmarglike;
    double acceptance;
    
    // storage matrices
    Matrix<> P_store(nstore, ns*ns);
    Matrix<> ps_store(n, ns);
    Matrix<int> s_store(nstore, n);
    
    Matrix<> yhat_mat(nstore, n);
    Matrix<> yerror_mat(nstore, n);
    Matrix<> yfore_pred_mat(nstore, n);
    Matrix<> yback_pred_mat(nstore, n);

    if (*change == 1){    
      // fixed effects change only
      Matrix <> beta(ns, k, betastart);
      Matrix<> beta_store(nstore, ns*k);
      Matrix<> Sigma_store(nstore, 1);
      MCMCPACK_PASSRNG2MODEL(MCMCintervention_fixed_impl, 
			     *m, 
			     *intervention, 
			     Y, X, 
			     beta, Sigma, P, s, b0, B0,   
			     *c0, *d0, A0,
			     *burnin, *mcmc, *thin, *verbose, 
			     *chib, *ar, 
			     beta_store, Sigma_store, 
			     P_store, ps_store, s_store, 
			     logmarglike, 
			     yhat_mat, 
			     yerror_mat, 
			     yfore_pred_mat, 
			     yback_pred_mat, 
			     acceptance);	
      for (int i = 0; i<(nstore*ns*k); ++i){
	betaout[i] = beta_store[i]; 
      }
      for (int i = 0; i<(nstore); ++i){
	Sigmaout[i] = Sigma_store[i]; 
      }  
    }
    else if (*change == 2){
      // random effects change only
      Matrix<> beta(k, 1, betastart);
      Matrix<> beta_store(nstore, k);
      Matrix<> Sigma_store(nstore, ns);
      MCMCPACK_PASSRNG2MODEL(MCMCintervention_random_impl, 
			     *m, 
			     *intervention, 
			     Y, X, 
			     beta, Sigma, P, s, b0, B0,   
			     *c0, *d0, A0,
			     *burnin, *mcmc, *thin, *verbose, 
			     *chib, *ar, 
			     beta_store, Sigma_store, 
			     P_store, ps_store, s_store, 
			     logmarglike, 
			     yhat_mat, 
			     yerror_mat, 
			     yfore_pred_mat, 
			     yback_pred_mat)  ;
      for (int i = 0; i<(nstore*k); ++i){
	betaout[i] = beta_store[i]; 
      }
      for (int i = 0; i<(nstore*ns); ++i){
	Sigmaout[i] = Sigma_store[i]; 
      }  
      acceptance = 1;
    }
    else {
      Matrix <> beta(ns, k, betastart);
      Matrix<> beta_store(nstore, ns*k);
      Matrix<> Sigma_store(nstore, ns);
      MCMCPACK_PASSRNG2MODEL(MCMCintervention_impl, 
			     *m, 
			     *intervention, 
			     Y, X, 
			     beta, Sigma, P, s, b0, B0,   
			     *c0, *d0, A0,
			     *burnin, *mcmc, *thin, *verbose, 
			     *chib, *ar, 
			     beta_store, Sigma_store, 
			     P_store, ps_store, s_store, 
			     logmarglike, 
			     yhat_mat, 
			     yerror_mat, 
			     yfore_pred_mat, 
			     yback_pred_mat, 
			     acceptance);
      // return output
      for (int i = 0; i<(nstore*ns*k); ++i){
	betaout[i] = beta_store[i]; 
      }
      for (int i = 0; i<(nstore*ns); ++i){
	Sigmaout[i] = Sigma_store[i]; 
      }
    }
    
    logmarglikeholder[0] = logmarglike;
    accept[0] = acceptance;

    for (int i = 0; i<(nstore*ns*ns); ++i){
      Pout[i] = P_store[i]; 
    }
    for (int i = 0; i<(n*ns); ++i){
      psout[i] = ps_store[i]; 
    }	   
    for (int i = 0; i<(nstore*n); ++i){
      sout[i] = s_store[i];
      yhatout[i] = yhat_mat[i];
      yerrorout[i] = yerror_mat[i];
      yforepredout[i] = yfore_pred_mat[i];
      ybackpredout[i] = yback_pred_mat[i];
    }
    
  }// end of MCMCpoissonChange
}//end extern "C"

#endif



  
