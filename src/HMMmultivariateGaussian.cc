//////////////////////////////////////////////////////////////////////////
// HMMmultivariateGaussian.cc is C++ code to estimate a Gaussian panel model with a structural break
// y_{it} = \x'_{it}\b + \varepsilon_{it}
// \varepsilon_{it} \sim \normdist{0}{\sigma^2}
// Parameters = {beta, sigma2, P_i}// static + changing
// 
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr

// Written 11/19/2008
// Modified 07/04/2011
//
//////////////////////////////////////////////////////////////////////////
#ifndef HMMMULTIVARIATEGAUSSIAN_CC
#define HMMMULTIVARIATEGAUSSIAN_CC

#include<vector>
#include<algorithm>

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
#include "lapack.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;

static double ln_invgamma(double theta, double a, double b) {
  double logf =  a * log(b) - lngammafn(a) + -(a+1) * log(theta) + 
                 -b/theta;
  return (logf);
  //pow(b, a) / gammafn(a) * pow(theta, -(a+1)) * exp(-b/theta);
}


// used to access 1d arrays holding R matrices like a 2d array 
#define M(ROW,COL,NROWS) (COL*NROWS+ROW)

template <typename RNGTYPE>
void MultivariateGaussian_impl(rng<RNGTYPE>& stream, 
			        int nsubj,  int ntime,  int nobs, 
			       const Matrix<int>& subjectid, const Matrix<int>& timeid, 
			       const Matrix<>& Y, const Matrix<>& X,  
			       const Matrix<>& YT, const Matrix<>& XT, 
			        int burnin,  int  mcmc, 
			        int thin,  int  verbose,  
			       double sigma2start,
			       const Matrix<>& b0, const Matrix<>& B0, 
			       const double c0, const double d0, 
			       const Matrix<>& time_groupinfo, 
			       const Matrix<>& subject_groupinfo, 
			       Matrix<>& beta_store, Matrix<>& sigma_store, 
			       double& logmarglike, double& loglike, 
			        int chib)
{ // redefine constants
  const int K = X.cols();; // ncol(X)
  const int NOBS = nobs;
  const Matrix<> B0inv = invpd(B0);

  const  int tot_iter = burnin + mcmc;
  const  int nstore = mcmc / thin; // number of draws to store
  
  // generate posk_arr and post_arr
  // Number of observations by group k
  int *nobsk = new int[nsubj];
  for (int k=0; k<nsubj; k++) {
    nobsk[k]=0;
    for (int n=0; n<NOBS; n++) {
      if (subjectid[n]==k+1) {
	nobsk[k]+=1;
      }
    }
  }

  // Position of each group in the data-set
  int **posk_arr = new int*[nsubj];
  for (int k=0; k<nsubj; k++) {
    posk_arr[k] = new int[nobsk[k]];
    int repk=0;
    for (int n=0; n<NOBS; n++) {
      if (subjectid[n]==k+1) {
	posk_arr[k][repk]=n;
	repk++;
      }
    }
  }

  // Small fixed matrices indexed on k for data access
  Matrix<double> *Yk_arr = new Matrix<double>[nsubj];
  Matrix<double> *Xk_arr = new Matrix<double>[nsubj];
  for(int k=0; k<nsubj; k++) {
    Xk_arr[k] = Matrix<double>(nobsk[k], K);
    Yk_arr[k] = Matrix<double>(nobsk[k], 1);
    for (int l=0; l<nobsk[k]; l++) {
      for (int p=0; p< K; p++) {
	Xk_arr[k](l, p) = X[p*NOBS + posk_arr[k][l]];
      }
      Yk_arr[k](l,0) = Y[posk_arr[k][l]];
    }
  } 
 
  Matrix<double> *cpXk_arr = new Matrix<double>[nsubj];
  Matrix<double> *tXYk_arr = new Matrix<double>[nsubj];
  for(int k=0; k<nsubj; k++) {
    cpXk_arr[k] = crossprod(Xk_arr[k]);
    tXYk_arr[k] = t(Xk_arr[k])*Yk_arr[k];
  }
  
  // starting values
  Matrix<> beta(K, 1);
  double sigma2 = sigma2start;
  
  // MCMC iterations start here 
  int sampcount = 0;

  Rprintf("\n ///////////////////////////////////////////////// \n");
  Rprintf("\n MCMC for Multivariate Gaussian loop starts! \n");
  Rprintf("\n ///////////////////////////////////////////////// \n");
  
  /////////////////////////////////////////////////
  // initialize Yhat for the first loop only
  /////////////////////////////////////////////////
  for ( int iter=0; iter < tot_iter; ++iter){
    //////////////////////
    // Step 1. Sample beta 
    //////////////////////
    Matrix<> XVX(K, K);
    Matrix<> XVY(K, 1);
    for( int s = 0; s<nsubj; ++s) {       
      XVX = XVX + cpXk_arr[s];
      XVY = XVY + tXYk_arr[s];       
    }
    Matrix<> beta_var = invpd(B0 + XVX/sigma2);
    Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma2);
    beta = stream.rmvnorm(beta_mean, beta_var);
    

    /////////////////////////////////////////////////
    // Step 2. Sample sigma2
    /////////////////////////////////////////////////
    double SSE = 0;   
    int counter = 0;
    for( int s=0;s<nsubj; ++s){
      int ntime_s = subject_groupinfo(s, 2);
      Matrix<> e = t(Yk_arr[s]-Xk_arr[s]*beta)*(Yk_arr[s] - Xk_arr[s]*beta);
      SSE = SSE + e[0];
      counter = counter + ntime_s;
    }
    double nu = (c0 + NOBS)/2;
    double delta = (d0 + SSE)/2;
    sigma2 = stream.rigamma(nu, delta);
   
    /////////////////////////////////////////////////
    // STORE
    /////////////////////////////////////////////////
    if (iter >= burnin && ((iter % thin) == 0)) { 
      for (int i=0; i< K; ++i){
	beta_store(sampcount, i) = beta(i);// stored by the order of (11, 12, 13, 21, 22, 23)
      }
      sigma_store(sampcount) = sigma2; 
      ++sampcount;
    }
  
      
    /////////////////////////////////////////////////
    // REPORT
    /////////////////////////////////////////////////
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n ----------------------------------------------------------------------- \n");
      Rprintf("Multivaraite Gaussian iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("\n beta is ");
      for(int i=0;i<K; ++i) {
	Rprintf("%10.5f", beta(i));
      }
      Rprintf("\n sigma2 = %10.5f\n", sigma2);
    }      
  }// END of MCMC loop
  
  //////////////////////////////////////
  if(chib == 1){
    //////////////////////////////////////
     Matrix<> beta_st(K, 1);
     beta_st(_, 0) = meanc(beta_store);       
     const double sigma2_st = mean(sigma_store);       
     
    Matrix<> density_beta(nstore, 1);  
    //////////////////////////////////////////////////////////////////
    // 1. pdf.beta | sigma2_g
    //////////////////////////////////////////////////////////////////
    for ( int iter = 0; iter<nstore; ++iter){
      Matrix<> XVX(K, K);
      Matrix<> XVY(K, 1);
      for( int s = 0; s<nsubj; ++s) {       
	XVX = XVX + cpXk_arr[s];
	XVY = XVY + tXYk_arr[s];       
      }
      Matrix<> beta_var = invpd(B0 + XVX/sigma_store(iter));
      Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma_store(iter));
      beta = stream.rmvnorm(beta_mean, beta_var);
      
      if (K == 1){
	density_beta(iter) = dnorm(beta_st(0), beta_mean[0], sqrt(beta_var[0]));   
      }
      else{
	density_beta(iter) = ::exp(lndmvn(beta_st, beta_mean, beta_var));
      }
    }
    double pdf_beta = log(mean(density_beta));   
      
    //////////////////////////////////////////////////////////////////
    // 2. pdf.sigma2
    ////////////////////////////////////////////////////////////////// 
    Matrix<> density_sigma2(nstore, 1);  
    for ( int iter = 0; iter<nstore; ++iter){
      double SSE = 0;   
      for( int s=0;s<nsubj; ++s){
	Matrix<> e = t(Yk_arr[s]-Xk_arr[s]*beta_st)*(Yk_arr[s] - Xk_arr[s]*beta_st);
	SSE = SSE + e[0];
      }
      double nu = (c0 + NOBS)/2;
      double delta = (d0 + SSE)/2;
      sigma2 = stream.rigamma(nu, delta);
      density_sigma2(iter) = ::exp(ln_invgamma(sigma2_st, nu, delta));         
    }
    double pdf_sigma2 = log(mean(density_sigma2));
  
    //////////////////////////////////////////////////////////////////
    // likelihood f(y|beta_st, D_st, sigma2_st, P_st)
    //////////////////////////////////////////////////////////////////  
    loglike = 0;
    for( int s = 0; s<nsubj; ++s) {       
      int ntime_s = subject_groupinfo(s, 2);     
      Matrix<> Sigma = sigma2_st*eye(ntime_s);
      Matrix<> Mu = Xk_arr[s]*beta_st;
      loglike += lndmvn(Yk_arr[s], Mu, Sigma); 
    }
    
    //////////////////////
    // log prior ordinate 
    //////////////////////
      double density_beta_prior = 0;
     if (K == 1){
      density_beta_prior = log(dnorm(beta_st(0), b0[0], sqrt(B0inv[0])));   
     }
     else{
       density_beta_prior = lndmvn(beta_st, b0, B0inv); 
     }   
     double density_sigma2_prior = ln_invgamma(sigma2_st, c0/2, d0/2);
     
     // compute marginal likelihood
     double logprior = density_beta_prior + density_sigma2_prior; 
     logmarglike = (loglike + logprior) - (pdf_beta + pdf_sigma2);
     if(verbose > 0){
       Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
       Rprintf("loglike = %10.5f\n", loglike);
       Rprintf("logprior = %10.5f\n", logprior);
       Rprintf("pdf_beta = %10.5f\n", pdf_beta);
       Rprintf("pdf_Sigma = %10.5f\n", pdf_sigma2);
       // Rprintf("pdf_P = %10.5f\n", pdf_P);
     }      
  }
  
  delete[] nobsk;
  for( int k=0; k<nsubj; k++) {
    delete[] posk_arr[k];
  }
  delete [] posk_arr;
  delete [] Xk_arr;
  delete [] Yk_arr;
  delete [] tXYk_arr;
  delete [] cpXk_arr;
  
}// end of MultivariateGaussian_impl


template <typename RNGTYPE>
void HMMmultivariateGaussian_impl(rng<RNGTYPE>& stream, 
				   int nsubj,  int ntime, 
				   int m,   int nobs, 
				  const Matrix<int>& subjectid, const Matrix<int>& timeid, 
				  const Matrix<>& Y, const Matrix<>& X,  
				  const Matrix<>& YT, const Matrix<>& XT, 
				   int burnin,  int  mcmc, 
				   int thin,  int  verbose,  
				  Matrix<>& betastart, double sigma2start,
				  const Matrix<>& b0, const Matrix<>& B0, 
				  const double c0, const double d0, 
				  const Matrix<>& P0,
				  const Matrix<>& time_groupinfo, 
				  const Matrix<>& subject_groupinfo, 
				  Matrix<>& beta_store, Matrix<>& sigma_store, 
				  Matrix<>& P_store,  Matrix<>& ps_store, Matrix<>& s_store, 
				  double& logmarglike, double& loglike, 
				   int chib)
{ // redefine constants
  const int K = X.cols();; // ncol(X)
  const int NOBS = nobs;
  const Matrix<> B0inv = invpd(B0);

  const  int tot_iter = burnin + mcmc;
  const  int nstore = mcmc / thin; // number of draws to store
  const int ns = m + 1; 
  
  // generate posk_arr and post_arr
  // Number of observations by group k
  int *nobsk = new int[nsubj];
  for (int k=0; k<nsubj; k++) {
    nobsk[k]=0;
    for (int n=0; n<NOBS; n++) {
      if (subjectid[n]==k+1) {
	nobsk[k]+=1;
      }
    }
  }

  // Position of each group in the data-set
  int **posk_arr = new int*[nsubj];
  for (int k=0; k<nsubj; k++) {
    posk_arr[k] = new int[nobsk[k]];
    int repk=0;
    for (int n=0; n<NOBS; n++) {
      if (subjectid[n]==k+1) {
	posk_arr[k][repk]=n;
	repk++;
      }
    }
  }

  // Small fixed matrices indexed on k for data access
  Matrix<double> *Yk_arr = new Matrix<double>[nsubj];
  Matrix<double> *Xk_arr = new Matrix<double>[nsubj];
  for(int k=0; k<nsubj; k++) {
    Xk_arr[k] = Matrix<double>(nobsk[k], K);
    Yk_arr[k] = Matrix<double>(nobsk[k], 1);
    for (int l=0; l<nobsk[k]; l++) {
      for (int p=0; p< K; p++) {
	Xk_arr[k](l, p) = X[p*NOBS + posk_arr[k][l]];
      }
      Yk_arr[k](l,0) = Y[posk_arr[k][l]];
    }
  } 
  // number of observations by time t
  int *nobst = new int[ntime];
  for (int k=0; k<ntime; k++) {
    nobst[k]=0;
    for (int n=0; n<NOBS; n++) {
      if (timeid[n]==k + 1) {
	nobst[k]+=1;
      }
    }
  }
  // Position of each group in the data-set
  int **post_arr = new int*[ntime];
  for (int k=0; k<ntime; k++) {
    post_arr[k] = new int[nobst[k]];
    int repk=0;
    for (int n=0; n<NOBS; n++) {
      if (timeid[n]==k+1) {
	post_arr[k][repk]=n;
	repk++;
      }
    }
  }
     
  // XTarr is data transformed for multivariate TS analsysi
  Matrix<>* Xt_arr = new Matrix<>[ntime];
  Matrix<>* Yt_arr = new Matrix<>[ntime];
  for(int k=0; k<ntime; k++) {
    if (nobst[k] > 0){
      Xt_arr[k] = Matrix<double>(nobst[k], K);
      Yt_arr[k] = Matrix<double>(nobst[k], 1);
      for (int l = 0; l<nobst[k]; l++) {
	for (int p = 0; p < K; p++) {
	  Xt_arr[k](l, p) = X[p*NOBS+post_arr[k][l]];
	}
	Yt_arr[k](l,0) = Y[post_arr[k][l]];
      }
    }
  } 

  Matrix<double> *cpXt_arr = new Matrix<double>[ntime];
  Matrix<double> *tXt_arr = new Matrix<double>[ntime];
  Matrix<double> *tXYt_arr = new Matrix<double>[ntime];
  for(int k=0; k<ntime; k++) {
    cpXt_arr[k] = crossprod(Xt_arr[k]);
    tXt_arr[k] = t(Xt_arr[k]);
    tXYt_arr[k] = t(Xt_arr[k])*Yt_arr[k];
  }
  
  // starting values
  Matrix<> beta(ns, K);
  Matrix<> sigma2(ns, 1);
  for (int j=0; j<ns; ++j){
    beta(j,_) = betastart;
    sigma2(j) = sigma2start;
  }
  Matrix<> P = P0;
  
  // MCMC iterations start here 
  int sampcount = 0;

  /////////////////////////////////////////////////
  // initialize Yhat for the first loop only
  /////////////////////////////////////////////////
  for ( int iter=0; iter < tot_iter; ++iter){
    //////////////////////
    // Step 1. Sample state
    //////////////////////  
    Matrix<> F(ntime, ns);
    Matrix<> pr1(ns, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);
    for ( int tt=0; tt<ntime ; ++tt) {
      if(nobst[tt]>0){
	int nsubj_s = time_groupinfo(tt, 2);
	for (int j=0; j<ns;++j){
	  Matrix<> Sigma = eye(nsubj_s)*sigma2[j];
	  Matrix<> Mu = Xt_arr[tt]*::t(beta(j,_));
	  py[j] = ::exp(lndmvn(Yt_arr[tt], Mu, Sigma));
	}	
	if (tt==0) pstyt1 = pr1;
	else {
	  pstyt1 =  ::t(F(tt-1,_)*P); 
	}	
	Matrix<> unnorm_pstyt = pstyt1%py;
	Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt); 
	for (int j=0; j<ns ; ++j) {
	  F(tt,j) = pstyt(j);
	}
      }	
    }// end of F matrix filtering
    Matrix<int> state(ntime, 1);        
    Matrix<> ps = Matrix<>(ntime, ns); 
    ps(ntime-1,_) = F(ntime-1,_);                      
    state(ntime-1) = ns;                                
    
    Matrix<> pstyn = Matrix<>(ns, 1);
    double pone = 0.0;
    int tt = ntime-2;
    
    while (tt >= 0){
      if(nobst[tt]>0){
	int st = state(tt+1);
	Matrix<> Pst_1 = ::t(P(_, st-1)); 
	Matrix<> unnorm_pstyn = F(tt,_)%Pst_1;
	pstyn = unnorm_pstyn/sum(unnorm_pstyn); 
	
	if (st==1){		
	  state(tt) = 1;                  	  
	}
	else{
	  pone = pstyn(st-2);
	  if(stream.runif() < pone) state(tt) = st-1;
	  else state(tt) = st;
	}
	ps(tt,_) = pstyn;
	--tt;
      }
    }// end of while loop
   
    //////////////////////
    // Step 2. Sample beta 
    //////////////////////
    int beta_count = 0;
    Matrix<int> nstate(ns, 1); 
    
    for (int j = 0; j<ns; ++j){
      for ( int i = 0; i<ntime; ++i){
	if (state[i] == j + 1) { 
	  nstate[j] = nstate[j] + 1;
	}// end of if
      }// end of int i<n
      beta_count = beta_count + nstate[j];        
 
      // BETA UPDATE
      Matrix<> XVX(K, K);
      Matrix<> XVY(K, 1);
      for(int tt = (beta_count - nstate[j]); tt<beta_count; ++tt) { 
	if (nobst[tt] > 0){
	  XVX = XVX + cpXt_arr[tt];
	  XVY = XVY + tXYt_arr[tt];       
	}
      }
      Matrix<> beta_var = invpd(B0 + XVX/sigma2[j]);
      Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma2[j]);
      beta(j,_) = stream.rmvnorm(beta_mean, beta_var);
    } 

    /////////////////////////////////////////////////
    // Step 3. Sample sigma2
    /////////////////////////////////////////////////
    beta_count = 0;
    Matrix<int> YN(ns, 1); 
    Matrix<> SSE(ns, 1);
    for (int j = 0; j<ns; ++j){
      beta_count = beta_count + nstate[j];        
      for( int s = 0; s<nsubj; ++ s) {    
	Matrix<> yj = Yk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> Xj = Xk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), K-1);  	
	Matrix<> e = ::t(yj - Xj*::t(beta(j,_)))*(yj - Xj*::t(beta(j,_)));
	YN[j] = YN[j] + yj.rows();
	SSE[j] = SSE[j] + e[0];
      }
      double nu = c0 + (double)YN[j]*0.5;
      double scale = d0 + SSE[j]*0.5;      
      sigma2[j] = stream.rigamma(nu, scale);                      
    }
    
    //////////////////////
    // Step 4. Sample P
    //////////////////////   
    double shape1 = 0;
    double shape2 = 0;    
    for ( int j =0; j<m; ++j){
      shape1 =  std::abs(P0(j,j) + nstate[j] - 1);
      shape2 =  P0(j,j+1) + 1; //       
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }
    P(m, m) = 1; //no jump at the last state
    
    /////////////////////////////////////////////////
    // STORE
    /////////////////////////////////////////////////
    if (iter >= burnin && ((iter % thin) == 0)) { 
      Matrix<> tbeta = ::t(beta); 
      for (int i=0; i<(ns*K); ++i){
	beta_store(sampcount, i) = tbeta[i];
      }
      for (int i=0; i<ns; ++i){
	sigma_store(sampcount, i) = sigma2[i]; 
      }
      for (int j=0; j<ns*ns; ++j){    
	P_store(sampcount,j)= P[j];
      }
      s_store(sampcount,_) = state;
      for ( int l=0; l<ntime ; ++l){           
	ps_store(l,_) = ps_store(l,_) + ps(l,_);          
      }
      ++sampcount;
    }
  
      
    /////////////////////////////////////////////////
    // REPORT
    /////////////////////////////////////////////////
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n ----------------------------------------------------------------------- \n");
      Rprintf("HMM Multivaraite Gaussian iteration %i of %i \n", (iter+1), tot_iter);
      for(int j=0;j<ns; ++j) {
	Rprintf("\n number of obs from state %i is %i", j, nstate[j]);     
      }
      for(int j=0;j<ns; ++j) {
	Rprintf("\n beta for state %i = ", j);
	for(int i=0;i<K; ++i) {
	  Rprintf("%10.5f", beta(j, i));
	}
      }
      for(int j=0;j<ns; ++j) {
	Rprintf("\n sigma2 for state %i is %10.5f", j, sigma2[j]);
      }      
    }
    
  }// END of MCMC loop
  
  //////////////////////////////////////
  if(chib == 1){
    //////////////////////////////////////
    
    // posterior ordinate    
    Matrix<double> betast = meanc(beta_store);       
    Matrix<double, Row> beta_st(ns, K);
    for (int j = 0; j< ns*K; ++j){   
      beta_st[j] = betast[j];
    }    
    
    Matrix<double> sigma2st = meanc(sigma_store);       
    Matrix<double, Row> sigma2_st(ns, 1);
    for (int j = 0; j< ns; ++j){   
      sigma2_st[j] = sigma2st[j];
    }  
    
    Matrix<double> P_vec_st = meanc(P_store);
    const Matrix<double> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){  
      P_st[j] = P_vec_st[j]; 
    }  
       
    Matrix<> density_beta(nstore, ns);  
    //////////////////////////////////////////////////////////////////
    // 1. pdf.beta | sigma2_g, P_g
    //////////////////////////////////////////////////////////////////
    for ( int iter = 0; iter<nstore; ++iter){
      int beta_count = 0;
      Matrix<int> nstate(ns, 1);    
      for (int j = 0; j<ns; ++j){
	for ( int i = 0; i<ntime; ++i){
	  if (s_store(iter, i) == j + 1) { 
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
	beta_count = beta_count + nstate[j];        
	
	Matrix<> XVX(K, K);
	Matrix<> XVY(K, 1);
	for(int tt = (beta_count - nstate[j]); tt<beta_count; ++tt) { 
	  XVX = XVX + cpXt_arr[tt];
	  XVY = XVY + tXYt_arr[tt];       
	}
	Matrix<> beta_var = invpd(B0 + XVX/sigma_store(iter, j));
	Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma_store(iter, j));
	beta(j,_) = stream.rmvnorm(beta_mean, beta_var);
	
	if (K == 1){
	  density_beta(iter, j) = ::exp(log(dnorm(beta_st(j), beta_mean[0], sqrt(beta_var[0]))));   
	}
	else{
	  density_beta(iter, j) = ::exp(lndmvn(::t(beta_st(j,_)), beta_mean, beta_var));
	}
      } 
    }
    double pdf_beta = log(prod(meanc(density_beta)));   
      
    //////////////////////////////////////////////////////////////////
    // 2. pdf.sigma2
     ////////////////////////////////////////////////////////////////// 
    Matrix<> density_sigma2(nstore, ns);  
    for ( int iter = 0; iter<nstore; ++iter){
      Matrix<> F(ntime, ns);
      Matrix<> pr1(ns, 1);
      pr1[0] = 1;
      Matrix<> py(ns, 1);
      Matrix<> pstyt1(ns, 1);
      
      for ( int tt=0; tt<ntime ; ++tt) {
	int nsubj_s = time_groupinfo(tt, 2);
	for (int j=0; j<ns;++j){
	  Matrix<> Sigma = eye(nsubj_s)*sigma2[j];
	  Matrix<> Mu = Xt_arr[tt]*::t(beta_st(j,_));
	  py[j]  =  ::exp(lndmvn(Yt_arr[tt], Mu, Sigma));
	}	  
	if (tt==0) pstyt1 = pr1;
	else {
	  pstyt1 =  ::t(F(tt-1,_)*P); 
	}
	
	Matrix<> unnorm_pstyt = pstyt1%py;
	Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt); 
	for (int j=0; j<ns ; ++j) {
	  F(tt,j) = pstyt(j);
	}
      }// end of F matrix filtering
      Matrix<int> state(ntime, 1);        
      Matrix<> ps = Matrix<>(ntime, ns); 
      ps(ntime-1,_) = F(ntime-1,_);                      
      state(ntime-1) = ns;                                   
      Matrix<> pstyn = Matrix<>(ns, 1);
      double pone = 0.0;
      int tt = ntime-2;    
      while (tt >= 0){
	int st = state(tt+1);
	Matrix<> Pst_1 = ::t(P(_, st-1)); 
	Matrix<> unnorm_pstyn = F(tt,_)%Pst_1;
	pstyn = unnorm_pstyn/sum(unnorm_pstyn); 	  
	if (st==1){		
	  state(tt) = 1;                  	  
	}
	else{
	  pone = pstyn(st-2);
	  if(stream.runif() < pone) state(tt) = st-1;
	  else state(tt) = st;
	}
	ps(tt,_) = pstyn;
	--tt;
	
      }// end of while loop
      
     /////////////////////////////////////////////////
      // 2.2. sigma2| beta_st, s, P
      /////////////////////////////////////////////////
      Matrix<int> nstate(ns, 1); 
      for (int j = 0; j<ns; ++j){
	for ( int i = 0; i<ntime; ++i){
	  if (state(i) == j + 1) { 
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
      }
      
      Matrix<> SSE(ns, 1); 
      int beta_count = 0;
      Matrix<int> YN(ns, 1); 
      for (int j = 0; j<ns; ++j){
	beta_count = beta_count + nstate[j];      
	for( int s = 0; s<nsubj; ++ s) {    
	  Matrix<> yj = Yk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  Matrix<> Xj = Xk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), K-1);  
	  YN[j] = YN[j] + yj.rows();
	  Matrix<> e = ::t(yj - Xj*::t(beta_st(j,_)))*(yj - Xj*::t(beta_st(j,_)));
	  SSE[j] = SSE[j] + e[0];
	}	
	double nu = (c0 + (double)YN[j])/2;
	double scale = (d0 + SSE[j])/2;      
	sigma2[j] = stream.rigamma(nu, scale);  
	
	density_sigma2(iter, j) = ::exp(ln_invgamma(sigma2_st[j], nu, scale));    
      }
      
      /////////////////////////////////////////////////
      // 2.3. P| beta_st, sigma2, s
      /////////////////////////////////////////////////
      double shape1 = 0;
      double shape2 = 0;    
      for ( int j =0; j<m; ++j){
	shape1 =  std::abs(P0(j,j) + nstate[j] - 1);
	shape2 =  P0(j,j+1) + 1; //       
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
      }
      P(m, m) = 1;       
    }
    double pdf_sigma2 = log(prod(meanc(density_sigma2)));
    
    //////////////////////
    // 3. pdf.P: 
    //////////////////////      
    Matrix<> density_P(nstore, ns);  
    //  Matrix<> density_local_P(ns, 1);  
    
    for ( int iter = 0; iter < nstore; ++iter){     
      Matrix<> F(ntime, ns);
      Matrix<> pr1(ns, 1);
      pr1[0] = 1;
      Matrix<> py(ns, 1);
      Matrix<> pstyt1(ns, 1);     
      for ( int tt=0; tt<ntime ; ++tt) {
	if(nobst[tt]>0){
	  int nsubj_s = time_groupinfo(tt, 2);
	  for (int j=0; j<ns;++j){
	    Matrix<> Sigma = eye(nsubj_s)*sigma2_st[j];
	    Matrix<> Mu = Xt_arr[tt]*::t(beta_st(j,_));
	    py[j]  =  ::exp(lndmvn(Yt_arr[tt], Mu, Sigma));
	  }	  
	  if (tt==0) pstyt1 = pr1;
	  else {
	    pstyt1 =  ::t(F(tt-1,_)*P); 
	  }	  
	  Matrix<> unnorm_pstyt = pstyt1%py;
	  Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt); 
	  for (int j=0; j<ns ; ++j) {
	    F(tt,j) = pstyt(j);
	  }
	}	
      }// end of F matrix filtering
      Matrix<int> state(ntime, 1);        
      Matrix<> ps = Matrix<>(ntime, ns);  
      ps(ntime-1,_) = F(ntime-1,_);                      
      state(ntime-1) = ns;                                   
      Matrix<> pstyn = Matrix<>(ns, 1);
      double pone = 0.0;
      int tt = ntime-2;
      
      while (tt >= 0){
	if(nobst[tt]>0){
	  int st = state(tt+1);
	  Matrix<> Pst_1 = ::t(P(_, st-1)); 
	  Matrix<> unnorm_pstyn = F(tt,_)%Pst_1;
	  pstyn = unnorm_pstyn/sum(unnorm_pstyn); 
	  if (st==1){		
	    state(tt) = 1;                  	  
	  }
	  else{
	    pone = pstyn(st-2);
	    if(stream.runif() < pone) state(tt) = st-1;
	    else state(tt) = st;
	  }
	  ps(tt,_) = pstyn;
	  --tt;
	}
      }// end of while loop
      
      Matrix<int> nstate(ns, 1); 
      for (int j = 0; j<ns; ++j){
	for ( int i = 0; i<ntime; ++i){
	  if (state(i) == j + 1) { 
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
      }
      double shape1 = 0;
      double shape2 = 0;    
      for ( int j =0; j<m; ++j){
	shape1 =  std::abs(P0(j,j) + nstate[j] - 1);
	shape2 =  P0(j,j+1) + 1; //       
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
	density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2); 
      }
      density_P(iter, ns-1) = 1; 
    }         
    
    double pdf_P = log(prod(meanc(density_P)));
    
    
    //////////////////////////////////////////////////////////////////
    // likelihood f(y|beta_st, D_st, sigma2_st, P_st)
    //////////////////////////////////////////////////////////////////  
    Matrix<> F(ntime, ns);
    Matrix<> pr1(ns, 1);
    Matrix<> like(ntime, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);
    for ( int tt=0; tt<ntime ; ++tt) {
      for (int j=0; j<ns;++j){
	Matrix<> Sigma = eye(nsubj)*sigma2_st[j];
	Matrix<> Mu = Xt_arr[tt]*::t(beta_st(j,_));
	py[j]  =  ::exp(lndmvn(Yt_arr[tt], Mu, Sigma));
      }
      if (tt==0) pstyt1 = pr1;
      else {
	pstyt1 =  ::t(F(tt-1,_)*P_st); 
      }   
      Matrix<> unnorm_pstyt = pstyt1%py;
      Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
      for (int j=0; j<ns ; ++j) {
	F(tt,j) = pstyt(j);
      }
      like[tt] = sum(unnorm_pstyt);
    }// end of F matrix filtering
    
    loglike = sum(log(like));
    
    //////////////////////
    // log prior ordinate 
    //////////////////////
    Matrix<double> density_beta_prior(ns, 1);
    Matrix<double> density_sigma2_prior(ns, 1);
    Matrix<double> density_P_prior(ns, 1);
    density_P_prior[ns-1] = 0; //
    
    for (int j=0; j<ns ; ++j){
      if (K == 1){
	density_beta_prior[j] = log(dnorm(beta_st(j), b0[0], sqrt(B0inv[0])));   
      }
      else{
	density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv);    
      }
      density_sigma2_prior[j] = ln_invgamma(sigma2_st[j], c0/2, d0/2);
    }
    for ( int j=0; j<m ; ++j){
      density_P_prior[j] = log(dbeta(P_st(j,j), P0(j,j), P0(j,j+1))); 
    }        
    
    // compute marginal likelihood
    double logprior = sum(density_beta_prior) + sum(density_sigma2_prior) + 
      sum(density_P_prior);
    logmarglike = (loglike + logprior) - (pdf_beta + pdf_P + pdf_sigma2);
    if(verbose > 0){
       Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
       Rprintf("loglike = %10.5f\n", loglike);
       Rprintf("logprior = %10.5f\n", logprior);
       Rprintf("pdf_beta = %10.5f\n", pdf_beta);
       Rprintf("pdf_Sigma = %10.5f\n", pdf_sigma2);
       Rprintf("pdf_P = %10.5f\n", pdf_P);
     }      
      
  }// end of marginal likelihood computation
  
  delete[] nobst;
  for(int k=0; k<ntime; k++) {
    delete[] post_arr[k];
  }
  delete[] nobsk;
  for(int k=0; k<nsubj; k++) {
    delete[] posk_arr[k];
  }
  delete [] posk_arr;
  delete [] Xk_arr;
  delete [] Yk_arr;
  delete [] Xt_arr;
  delete [] tXt_arr;
  delete [] Yt_arr;
  delete [] cpXt_arr;
  delete [] tXYt_arr;
  
  
}// end of HMMmultivariateGaussian_impl

extern "C" {
  void HMMmultivariateGaussian(double* betadata, const int* betarow, const int* betacol,
			       double* sigmadata, double *psout,
			       const int* nsubj, const int* ntime, 
			       const int* m, const int* nobs, const int* subjectid, const int* timeid, 
			       const double* Ydata, const int* Yrow, const int* Ycol,
			       const double* Xdata, const int* Xrow, const int* Xcol,
			       const double* YTdata, const double* XTdata,
			       const int* burnin, const int* mcmc, const int* thin, const int* verbose,
			       const int *uselecuyer, const int *seedarray, const int *lecuyerstream, 
			       const double* betastartdata, const double* sigma2start,
			       const double* b0data, const double* B0data, 
			       const double* c0, const double* d0, 
			       const double* P0data, const int* P0row, const int* P0col,
			       const double* subject_groupinfodata, const double* time_groupinfodata, 
			       double *logmarglikeholder, double *loglikeholder, 
			       const int *chib){
    // pull together Matrix objects
    Matrix<> Y(*Yrow, *Ycol, Ydata);
    Matrix<> X(*Xrow, *Xcol, Xdata);
    Matrix<> YT(*Yrow, *Ycol, YTdata);
    Matrix<> XT(*Xrow, *Xcol, XTdata);
    Matrix<> betastart(*Xcol, 1, betastartdata);
    Matrix<> b0(*Xcol, 1, b0data);
    Matrix<> B0(*Xcol, *Xcol, B0data);
    Matrix<> Pstart(*P0row, *P0col, P0data);
    Matrix<int> subjectid_mat(*nobs, 1, subjectid);
    Matrix<int> timeid_mat(*nobs, 1, timeid);
    Matrix<> subject_groupinfo(*nsubj, 3, subject_groupinfodata);
    Matrix<> time_groupinfo(*ntime, 3, time_groupinfodata);
    const int ns = *m + 1;   
    double logmarglike = 0.0;
    double loglike = 0.0;
    
    if (*m == 0) {
      Matrix<> beta_store(*betarow, *betacol);
      Matrix<> sigma_store(*betarow, 1);
      MCMCPACK_PASSRNG2MODEL(MultivariateGaussian_impl, 
			     *nsubj,  *ntime,  *nobs, 
			     subjectid_mat, timeid_mat, 
			     Y, X, YT, XT, 			   
			     *burnin, *mcmc, *thin, *verbose,			   
			     *sigma2start,  
			     b0, B0, *c0, *d0, 
			     time_groupinfo,  subject_groupinfo, 
			     beta_store, sigma_store,  
			     logmarglike, loglike, *chib);
      // store marginal likelihood
      logmarglikeholder[0] = logmarglike;
      loglikeholder[0] = loglike;    
      
      for (int i=0; i < (*betarow* *betacol); ++i){
	betadata[i] = beta_store(i);
      }
      for (int i=0; i < (*betarow*ns); ++i){
	sigmadata[i] = sigma_store(i);
      }
    }// end of if m == 0
    else {
      Matrix<> beta_store(*betarow, *betacol);
      Matrix<> sigma_store(*betarow, ns);
      Matrix<> P_store(*betarow,  ns*ns);
      Matrix<> s_store(*betarow,  *ntime);
      Matrix<> ps_store(*ntime, ns); 
      MCMCPACK_PASSRNG2MODEL(HMMmultivariateGaussian_impl, 
			     *nsubj,  *ntime,  *m, *nobs, 
			     subjectid_mat, timeid_mat, 
			     Y, X, YT, XT, 			   
			     *burnin, *mcmc, *thin, *verbose,			   
			     betastart, *sigma2start,  
			     b0, B0, *c0, *d0, 
			     Pstart,  time_groupinfo,  subject_groupinfo, 
			     beta_store, sigma_store, P_store, ps_store, s_store, 
			     logmarglike, loglike, *chib);
      // store marginal likelihood
      logmarglikeholder[0] = logmarglike;
      loglikeholder[0] = loglike;    
      
      for (int i=0; i < (*betarow* *betacol); ++i){
	betadata[i] = beta_store(i);
      }
      for (int i=0; i < (*betarow*ns); ++i){
	sigmadata[i] = sigma_store(i);
      }
      for (int i = 0; i<(*ntime *ns); ++i){
	psout[i] = ps_store[i]; 
      }
    
    }// end of m>0
  }// end of HMMmultivariateGaussian_CC
}// end of extern C

#endif /*HMMmultivariateGaussian_CC  */
