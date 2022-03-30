//////////////////////////////////////////////////////////////////////////
// HMMpanelRE.cc is C++ code to estimate a Gaussian panel model with a structural break
// y_{it} = \x'_{it}\b + \w'_{it}\bi_i + \z'_{it}\d_{s_{it}} +\varepsilon_{it}
// \varepsilon_{it} \sim \normdist{0}{\sigma^2}
// \bi_i \sim \normdist{0}{\D}
// Parameters = {beta, bi, D, sigma2, || {s_it}, delta, P_i}// static + changing
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr

// Modified (consistent with MCMChregress by Ghislain Vieilledent) on 07/04/2011
// Included in MCMCpack on 09/2011
//////////////////////////////////////////////////////////////////////////
#ifndef CHMMPANELRE_CC
#define CHMMPANELRE_CC

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

// used to access 1d arrays holding R matrices like a 2d array
#define M(ROW,COL,NROWS) (COL*NROWS+ROW)

static double lndwish (const Matrix<>& W,
	        int v,
	       const Matrix<>& S) {
  const int K = S.rows();
  double gammapart = 1;
  for (int i=0; i<K; ++i){
    double gammain = (v - i)/2;
    gammapart = gammapart * gammafn(gammain);// since i starts from 0!
  }
  double logdenom = log(gammapart) + (v * K/2)*log(2.0) + (K * (K - 1)/4)*log(M_PI);
  double detS = 0;
  double detW = 0;
  if(K==1){
    detS = S[0];
    detW = W[0];
  }
  else{
    detS = det(S);
    detW = det(W);
  }
  Matrix<> hold = invpd(S) * W;
  Matrix<> diaghold(K, 1);
  diaghold(_,0) = diag(hold);
  double tracehold = sum(diaghold);
  double lognum = -1*(v *.5)*log(detS) + (v - K - 1)/2*log(detW) - 1/2 * tracehold;
  return(lognum - logdenom);
}

static double lndinvgamma(double theta, double a, double b) {
  double logf =  a * log(b) - lngammafn(a) + -(a+1) * log(theta) +
                 -b/theta;
  return (logf);
  //pow(b, a) / gammafn(a) * pow(theta, -(a+1)) * exp(-b/theta);
}

template <typename RNGTYPE>
void GaussianPanelRE_impl(rng<RNGTYPE>& stream,
			  const int nsubj,  int ntime,
			   int nobs,
			  const Matrix<int>& subjectid, const Matrix<int>& timeid,
			  const Matrix<>& Y, const Matrix<>& X, const Matrix<>& W,
			   int burnin,  int  mcmc,
			   int thin,  int  verbose,
			  const double sigma2start, Matrix<>& D,
			  const Matrix<>& b0, const Matrix<>& B0,
			  const double c0, const double d0,
			   int r0, const Matrix<>& R0,
			  const Matrix<>& time_groupinfo,
			  const Matrix<>& subject_groupinfo,
			  Matrix<>& beta_store, Matrix<>& sigma_store,
			  Matrix<>& D_store,
			  double& logmarglike, double& loglike,
			   int chib)
{ // redefine constants
  const int K = X.cols();; // ncol(X)
  const int Q = W.cols(); // ncol(W)
  const int NOBS = nobs;
  double sigma2 = sigma2start;

  const Matrix<> R0inv = invpd(R0);
  const Matrix<> B0inv = invpd(B0);

  const int tot_iter = burnin + mcmc;
  const int nstore = mcmc / thin; // number of draws to store

  Matrix<double> Dinv = invpd(D);

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
  Matrix<double> *Wk_arr = new Matrix<double>[nsubj];
  for(int k=0; k<nsubj; k++) {
    Xk_arr[k] = Matrix<double>(nobsk[k], K);
    Wk_arr[k] = Matrix<double>(nobsk[k], Q);
    Yk_arr[k] = Matrix<double>(nobsk[k], 1);
    for (int l=0; l<nobsk[k]; l++) {
      for (int p=0; p< K; p++) {
	Xk_arr[k](l, p) = X[p*NOBS + posk_arr[k][l]];
	// Rprintf("\n Xk_arr ( %i,  %i) is %10.5f", l, p, Xk_arr[k](l, p));
      }
      for (int q=0; q < Q; q++) {
	Wk_arr[k](l, q) = W[q*NOBS + posk_arr[k][l]];
      }
      Yk_arr[k](l,0) = Y[posk_arr[k][l]];
    }
  }

  Matrix<double> *tXk_arr = new Matrix<double>[nsubj];
  Matrix<double> *tWk_arr = new Matrix<double>[nsubj];
  Matrix<double> *tXWk_arr = new Matrix<double>[nsubj];
  Matrix<double> *tWXk_arr = new Matrix<double>[nsubj];
  Matrix<double> *tWYk_arr = new Matrix<double>[nsubj];
  Matrix<double> *tXYk_arr = new Matrix<double>[nsubj];
  Matrix<double> *cpXk_arr = new Matrix<double>[nsubj];
  Matrix<double> *cpWk_arr = new Matrix<double>[nsubj];
  for(int k=0; k<nsubj; k++) {
    tXk_arr[k] = t(Xk_arr[k]);
    tXYk_arr[k] = t(Xk_arr[k])*Yk_arr[k];
    tXWk_arr[k] = t(Xk_arr[k])*Wk_arr[k];
    tWk_arr[k] = t(Wk_arr[k]);
    tWXk_arr[k] = t(Wk_arr[k])*Xk_arr[k];
    tWYk_arr[k] = t(Wk_arr[k])*Yk_arr[k];
    cpWk_arr[k] = crossprod(Wk_arr[k]);
    cpXk_arr[k] = crossprod(Xk_arr[k]);
  }
  // MCMC iterations start here
  int sampcount = 0;

  Rprintf("\n ///////////////////////////////////////////////// \n");
  Rprintf("\n MCMC for HMM Gaussian Panel Randome Effects loop starts! \n");
  Rprintf("\n ///////////////////////////////////////////////// \n");

  Matrix<> beta(K, 1);

  /////////////////////////////////////////////////
  // initialize Yhat for the first loop only
  /////////////////////////////////////////////////
  for (int iter=0; iter < tot_iter; ++iter){
    /////////////////////////////////////////////////
     // Step 1. Sample beta (fixed-effects coef)
     /////////////////////////////////////////////////
     Matrix<double> XVX(K, K);
     Matrix<double> XVY(K, 1);
     for(int s = 0; s<nsubj; ++s) {
       XVX += (1/sigma2)*cpXk_arr[s]-pow(1/sigma2,2)*tXWk_arr[s]*invpd(Dinv + cpWk_arr[s]*(1/sigma2))*tWXk_arr[s];
       XVY += (1/sigma2)*tXYk_arr[s]-pow(1/sigma2,2)*tXWk_arr[s]*invpd(Dinv + cpWk_arr[s]*(1/sigma2))*tWYk_arr[s];
     }
     Matrix<> beta_var = invpd(B0 + XVX/sigma2);
     Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma2);
     beta = stream.rmvnorm(beta_mean, beta_var);

    /////////////////////////////////////////////////
    // Step 2. Sample bi (random-effects coef)
    /////////////////////////////////////////////////
     Matrix<double> bi(Q, nsubj);
     for(int s = 0; s<nsubj; ++ s) {
       Matrix<double> b_var = invpd(Dinv + cpWk_arr[s]/sigma2);
       Matrix<double> b_mean = b_var*tWk_arr[s]*(Yk_arr[s] - Xk_arr[s]*beta)/sigma2;
       bi(_,s) = stream.rmvnorm(b_mean, b_var);
     }

     /////////////////////////////////////////////////
    // Step 3. Sample sigma2
    /////////////////////////////////////////////////
    double SSE = 0;
    for(int s=0;s<nsubj; ++s){
      Matrix<> e = t(Yk_arr[s]-Xk_arr[s]*beta - Wk_arr[s]*bi(_,s))*
	(Yk_arr[s] - Xk_arr[s]*beta - Wk_arr[s]*bi(_,s));
      SSE = SSE + e[0];
    }
    double nu = (c0 + NOBS)/2;
    double delta = (d0 + SSE)/2;
    sigma2 = stream.rigamma(nu, delta);

    /////////////////////////////////////////////////
    // Step 4. Sample D
    /////////////////////////////////////////////////
    Matrix<double> SSB = bi*t(bi);
    Matrix<double> D_scale = invpd(R0inv + SSB);
    int D_dof = r0 + nsubj;
    Dinv = stream.rwish(D_dof, D_scale);
    D = inv(Dinv);

     /////////////////////////////////////////////////
    // STORE
    /////////////////////////////////////////////////
    if (iter >= burnin && ((iter % thin) == 0)) {
      for(int j=0;j<K; ++j) {
	beta_store(sampcount,j) = beta(j);
      }
      for(int j=0;j<(Q*Q); ++j) {
	D_store(sampcount, j) = D(j);
      }
      sigma_store(sampcount) = sigma2;
      ++sampcount;
    }
    /////////////////////////////////////////////////
    // REPORT
    /////////////////////////////////////////////////
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n ----------------------------------------------------------------------- \n");
      Rprintf("\n Gaussian Panel iteration %i of %i ", (iter+1), tot_iter);
      // store global estimates
      Rprintf("\n beta = ");
      for(int j=0;j<K; ++j) {
	Rprintf("%10.5f", beta(j));
      }
      Rprintf("\n D[,] = ");
      for(int j=0;j<(Q*Q); ++j) {
	Rprintf("%10.5f", D(j));
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
     Matrix<> Dst = meanc(D_store);
     Matrix<> D_st(Q, Q);
     for (int k = 0; k<(Q*Q); ++k){
       D_st[k] = Dst[k];
     }
     const Matrix <> Dinv_st = invpd(D_st);
     Matrix<> beta_density(nstore, 1);

     //////////////////////////////////////////////////////////////////
     // 1. pdf.beta | D_g, sigma2_g
     //////////////////////////////////////////////////////////////////
     for (int iter=0; iter < nstore; ++iter){
       for(int j=0;j<(Q*Q); ++j) {
	 D(j) = D_store(iter,j);
       }
       Dinv = invpd(D);
       Matrix<double> XVX(K, K);
       Matrix<double> XVY(K, 1);
       for(int s = 0; s<nsubj; ++s) {
	 XVX += (1/sigma_store(iter))*cpXk_arr[s]-pow(1/sigma_store(iter),2)*tXWk_arr[s]*invpd(Dinv + cpWk_arr[s]*(1/sigma_store(iter)))*tWXk_arr[s];
	 XVY += (1/sigma_store(iter))*tXYk_arr[s]-pow(1/sigma_store(iter),2)*tXWk_arr[s]*invpd(Dinv + cpWk_arr[s]*(1/sigma_store(iter)))*tWYk_arr[s];
       }
       Matrix<> beta_var = invpd(B0 + XVX/sigma_store(iter));
       Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma_store(iter));
       if (K == 1){
	 beta_density(iter) = dnorm(beta_st(0), beta_mean[0], sqrt(beta_var[0]));
       }
       else{
	 beta_density(iter) = ::exp(lndmvn(beta_st, beta_mean, beta_var));
       }
     }
     double pdf_beta = log(mean(beta_density));

     //////////////////////////////////////////////////////////////////
     // 2. pdf.D
      //////////////////////////////////////////////////////////////////
     Matrix <> D_density(nstore, 1);
     for (int iter=0; iter < nstore; ++iter){

       /////////////////////////////////////////////////
       // Step 2.1 Sample bi (random-effects coef)
       /////////////////////////////////////////////////
       Matrix<double> bi(Q, nsubj);
       for(int s = 0; s<nsubj; ++ s) {
	 Matrix<double> b_var = invpd(Dinv + cpWk_arr[s]/sigma2);
	 Matrix<double> b_mean = b_var*t(Wk_arr[s])*(Yk_arr[s] - Xk_arr[s]*beta_st)/sigma2;
	 bi(_,s) = stream.rmvnorm(b_mean, b_var);
       }

       /////////////////////////////////////////////////
       // Step 2.2 Sample sigma2
       /////////////////////////////////////////////////
       double SSE = 0;
       for(int s=0;s<nsubj; ++s){
	 Matrix<> e = t(Yk_arr[s]-Xk_arr[s]*beta_st - Wk_arr[s]*bi(_,s))*
	   (Yk_arr[s] - Xk_arr[s]*beta_st - Wk_arr[s]*bi(_,s));
	 SSE = SSE + e[0];
       }
       double nu = (c0 + NOBS)/2;
       double delta = (d0 + SSE)/2;
       sigma2 = stream.rigamma(nu, delta);

       /////////////////////////////////////////////////
       // Step 2.3 Sample D
       /////////////////////////////////////////////////
       Matrix<double> SSB = bi*t(bi);
       Matrix<double> D_scale = invpd(R0inv + SSB);
       int D_dof = r0 + nsubj;
       Dinv = stream.rwish(D_dof, D_scale);
       D_density(iter) = ::exp(lndwish(invpd(D_st), D_dof, D_scale));
       D = inv(Dinv);
     }
     double pdf_D = log(mean(D_density));

     //////////////////////////////////////////////////////////////////
     // 3. pdf.sigma2
     //////////////////////////////////////////////////////////////////
     Matrix <> sigma_density(nstore, 1);
     for (int iter=0; iter < nstore; ++iter){

       /////////////////////////////////////////////////
       // Step 3.1 Sample bi (random-effects coef)
       /////////////////////////////////////////////////
       Matrix<double> bi(Q, nsubj);
       for(int s = 0; s<nsubj; ++ s) {
	 Matrix<double> b_var = invpd(Dinv_st + cpWk_arr[s]/sigma2);
	 Matrix<double> b_mean = b_var*tWk_arr[s]*(Yk_arr[s] - Xk_arr[s]*beta_st)/sigma2;
	 bi(_,s) = stream.rmvnorm(b_mean, b_var);
       }

       /////////////////////////////////////////////////
       // Step 3.2 Sample sigma2
       /////////////////////////////////////////////////
       double SSE = 0;
       for(int s=0;s<nsubj; ++s){
	 Matrix<> e = t(Yk_arr[s]-Xk_arr[s]*beta_st - Wk_arr[s]*bi(_,s))*
	   (Yk_arr[s] - Xk_arr[s]*beta_st - Wk_arr[s]*bi(_,s));
	SSE = SSE + e[0];
       }
       double nu = (c0 + NOBS)/2;
       double delta = (d0 + SSE)/2;
       sigma2 = stream.rigamma(nu, delta);
       sigma_density(iter) = ::exp(lndinvgamma(sigma2_st, nu, delta));
     }
     double pdf_sigma2 = log(mean(sigma_density));
     //////////////////////////////////////////////////////////////////
     // likelihood f(y|beta_st, D_st, sigma2_st)
     //////////////////////////////////////////////////////////////////
     loglike = 0;
     for(int s = 0; s<nsubj; ++s) {
       int ntime_s = subject_groupinfo(s, 2);
       Matrix<> Sigma = sigma2_st*eye(ntime_s) + Wk_arr[s] *D_st* t(Wk_arr[s]);
       Matrix<> Mu = Xk_arr[s]*beta_st;
       loglike += lndmvn(Yk_arr[s], Mu, Sigma);
     }

     //////////////////////
     // log prior ordinate
     //////////////////////
     double density_beta_prior = 0;
     if (K == 1){
       density_beta_prior =log(dnorm(beta_st(0), b0[0], sqrt(B0inv[0])));
     }
     else{
       density_beta_prior = lndmvn(beta_st, b0, B0inv);
     }
     double density_D_prior = lndwish(Dinv_st, r0, R0);
     double density_sigma2_prior = lndinvgamma(sigma2_st, c0/2, d0/2);

     // compute marginal likelihood
     double logprior = density_beta_prior + density_sigma2_prior + density_D_prior;
     logmarglike = (loglike + logprior) - (pdf_beta + pdf_sigma2 + pdf_D);

     if (verbose >0){
       Rprintf("\n ----------------------------------------------------------------------- \n");
       Rprintf("\n logmarglike %10.5f", logmarglike);
       Rprintf("\n loglike %10.5f", loglike, "\n");
       Rprintf("\n log_prior %10.5f", logprior, "\n");
       Rprintf("\n pdf_beta is %10.5f", pdf_beta, "\n");
       Rprintf("\n pdf_D is %10.5f", pdf_D, "\n");
       Rprintf("\n pdf_sigma2 is %10.5f", pdf_sigma2, "\n");
      }
   }// end of marginal likelihood computation

   delete [] Xk_arr;
   delete [] Yk_arr;
   delete [] Wk_arr;
   delete [] cpWk_arr;
   delete [] cpXk_arr;
   delete [] tWk_arr;
   delete [] tXk_arr;
   delete [] tXYk_arr;
   delete [] tXWk_arr;
   delete [] tWXk_arr;
   delete [] tWYk_arr;

}


template <typename RNGTYPE>
void HMMpanelRE_impl(rng<RNGTYPE>& stream,
		     int nsubj, int ntime,
		     int m,  int nobs,
		     const Matrix<int>& subjectid, const Matrix<int>& timeid,
		     const Matrix<>& Y, const Matrix<>& X, const Matrix<>& W,
		     const Matrix<>& YT, const Matrix<>& XT, const Matrix<>& WT,
		      int burnin,  int  mcmc,
		      int thin,  int  verbose,
		     Matrix<>& betastart, double sigma2start, Matrix<>& Dstart,
		     const Matrix<>& b0, const Matrix<>& B0,
		     const double c0, const double d0,
		      int r0, const Matrix<>& R0, const Matrix<>& P0,
		     const Matrix<>& time_groupinfo,
		     const Matrix<>& subject_groupinfo,
		     Matrix<>& beta_store, Matrix<>& sigma_store,
		     Matrix<>& D_store,
		     Matrix<>& P_store,  Matrix<>& ps_store, Matrix<>& s_store,
		     double& logmarglike, double& loglike,
		      int chib)
{ // redefine constants
  const int K = X.cols();
  const int Q = W.cols();
  const int NOBS = nobs;

  const Matrix<> R0inv = invpd(R0);
  const Matrix<> B0inv = invpd(B0);

  const int tot_iter = burnin + mcmc;
  const int nstore = mcmc / thin; // number of draws to store
  const int ns = m + 1;

  Matrix<>* D = new Matrix<>[ns];
  Matrix<>* Dinv = new Matrix<>[ns];
  for (int j=0; j<ns; ++j){
    D[j] = Dstart;
    Dinv[j] = invpd(Dstart);
  }

  Matrix<> D_record(ns, Q*Q);// for record D in D_store
  for (int j = 0; j<ns; ++j){
    for (int i=0; i<(Q*Q); ++i){
      D_record(j,i) = D[j](i);
    }
  }


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
  Matrix<double> *Wk_arr = new Matrix<double>[nsubj];
  for(int k=0; k<nsubj; k++) {
    Xk_arr[k] = Matrix<double>(nobsk[k], K);
    Wk_arr[k] = Matrix<double>(nobsk[k], Q);
    Yk_arr[k] = Matrix<double>(nobsk[k], 1);
    for (int l=0; l<nobsk[k]; l++) {
      for (int p=0; p< K; p++) {
	Xk_arr[k](l, p) = X[p*NOBS + posk_arr[k][l]];
      }
      for (int q=0; q < Q; q++) {
	Wk_arr[k](l, q) = W[q*NOBS + posk_arr[k][l]];
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
  Matrix<>* Wt_arr = new Matrix<>[ntime];
  Matrix<>* Yt_arr = new Matrix<>[ntime];
  for(int k=0; k<ntime; k++) {
    if (nobst[k] > 0){
      Xt_arr[k] = Matrix<double>(nobst[k], K);
      Wt_arr[k] = Matrix<double>(nobst[k], Q);
      Yt_arr[k] = Matrix<double>(nobst[k], 1);
      for (int l = 0; l<nobst[k]; l++) {
	for (int p = 0; p < K; p++) {
	  Xt_arr[k](l, p) = X[p*NOBS+post_arr[k][l]];
	}
	for (int q = 0; q < Q; q++) {
	  Wt_arr[k](l, q) = W[q*NOBS+post_arr[k][l]];
	}
	Yt_arr[k](l,0) = Y[post_arr[k][l]];
      }
    }
  }

  Matrix<double> *cpWt_arr = new Matrix<double>[ntime];
  Matrix<double> *cpXt_arr = new Matrix<double>[ntime];
  Matrix<double> *tXWt_arr = new Matrix<double>[ntime];
  Matrix<double> *tWXt_arr = new Matrix<double>[ntime];
  Matrix<double> *tWYt_arr = new Matrix<double>[ntime];
  Matrix<double> *tXYt_arr = new Matrix<double>[ntime];
  Matrix<double> *tXt_arr = new Matrix<double>[ntime];
  Matrix<double> *tWt_arr = new Matrix<double>[ntime];
  for(int k=0; k<ntime; k++) {
    cpXt_arr[k] = crossprod(Xt_arr[k]);
    cpWt_arr[k] = crossprod(Wt_arr[k]);
    tXt_arr[k] = t(Xt_arr[k]);
    tWt_arr[k] = t(Wt_arr[k]);
    tXWt_arr[k] = t(Xt_arr[k])*Wt_arr[k];
    tWXt_arr[k] = t(Wt_arr[k])*Xt_arr[k];
    tWYt_arr[k] = t(Wt_arr[k])*Yt_arr[k];
    tXYt_arr[k] = t(Xt_arr[k])*Yt_arr[k];
  }

  // starting values
  Matrix<double> *bk_run = new Matrix<double>[nsubj]; // Random effects
  for (int k=0;k<nsubj;k++) {
    bk_run[k] = Matrix<double>(Q,1);
  }
  Matrix<>* bi = new Matrix<>[ns];
  Matrix<> beta(ns, K);
  Matrix<> sigma2(ns, 1);
  for (int j=0; j<ns; ++j){
    bi[j] = Matrix<>(Q, nsubj);
    beta(j,_) = betastart;
    sigma2(j) = sigma2start;
  }
  Matrix<> P = P0;

  // MCMC iterations start here
  int sampcount = 0;
  // initialize Yhat for the first loop only
  for ( int iter=0; iter < tot_iter; ++iter){
    //////////////////////
    // Step 1. Sample state
    //////////////////////
    Matrix<> F(ntime, ns);
    Matrix<> pr1(ns, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);
    for (int tt=0; tt<ntime ; ++tt) {
      int nsubj_s = time_groupinfo(tt, 2);
      // Matrix<> Wbi(nsubj_s, 1);
      for (int j=0; j<ns;++j){
	Matrix<> Sigma = eye(nsubj_s)*sigma2[j];
	Matrix<> WDW = Wt_arr[tt]*D[j]*tWt_arr[tt];
	Matrix<> Mu = Xt_arr[tt]*::t(beta(j,_));
	py[j] = ::exp(lndmvn(Yt_arr[tt], Mu, WDW + Sigma));
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

    //////////////////////
    // Step 2. Sample beta
    //////////////////////
    Matrix<int> nstate(ns, 1); // refresh the container for numbers of each state
    int beta_count = 0;
    for (int j = 0; j<ns; ++j){
      for (int i = 0; i<ntime; ++i){
	if (state[i] == j + 1) {
	  nstate[j] = nstate[j] + 1;
	}
      }
      beta_count = beta_count + nstate[j];
      Matrix<> XVX(K, K);
      Matrix<> XVY(K, 1);
      for(int tt = (beta_count - nstate[j]); tt<beta_count; ++tt) {
	XVX += (1/sigma2[j])*cpXt_arr[tt]-pow(1/sigma2[j],2)*tXWt_arr[tt]*invpd(Dinv[j] + cpWt_arr[tt]*(1/sigma2[j]))*tWXt_arr[tt];
	XVY += (1/sigma2[j])*tXYt_arr[tt]-pow(1/sigma2[j],2)*tXWt_arr[tt]*invpd(Dinv[j] + cpWt_arr[tt]*(1/sigma2[j]))*tWYt_arr[tt];
      }
      Matrix<> beta_var = invpd(B0 + XVX/sigma2[j]);
      Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma2[j]);
      beta(j,_) = stream.rmvnorm(beta_mean, beta_var);
    }

    /////////////////////////////////////////////////
    // Step 3. Sample bi (random-effects coef)
    /////////////////////////////////////////////////
    beta_count = 0;
    Matrix<int> YN(ns, 1);
    Matrix<> SSE(ns, 1);
    for (int j = 0; j<ns; ++j){
      beta_count = beta_count + nstate[j];
      for(int s = 0; s<nsubj; ++ s) {
	Matrix<> yj = Yk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> Xj = Xk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), K-1);
	Matrix<> Wj = Wk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), Q-1);
	Matrix<> b_var = invpd(Dinv[j]+ ::t(Wj)*Wj/sigma2[j]);
	Matrix<> b_mean = b_var*::t(Wj)*(yj - Xj*::t(beta(j,_)))/sigma2[j];
	bi[j](_,s) = stream.rmvnorm(b_mean, b_var);

	// FOR SIGMA
	YN[j] = YN[j] + yj.rows();
	Matrix<> e = ::t(yj - Xj*::t(beta(j,_)) - Wj*bi[j](_,s))*(yj - Xj*::t(beta(j,_)) - Wj*bi[j](_,s));
	SSE[j] = SSE[j] + e[0];
      }
    }

    /////////////////////////////////////////////////
    // Step 4. Sample sigma2
    /////////////////////////////////////////////////
    for (int j = 0; j<ns; ++j){
      double nu = (c0 + (double)YN[j])*0.5;
      double scale = (d0 + SSE[j])*0.5;
      sigma2[j] = stream.rigamma(nu, scale);
    }

    /////////////////////////////////////////////////
    // Step 5. Sample D
    /////////////////////////////////////////////////
    for (int j = 0; j<ns; ++j){
      Matrix<> SSB = bi[j]*t(bi[j]);
      Matrix<> D_scale = invpd(R0inv + SSB);
      int D_dof = r0 + nsubj;
      Dinv[j] = stream.rwish(D_dof, D_scale);
      D[j] = invpd(Dinv[j]);

      for (int i=0; i<(Q*Q); ++i){
	D_record(j,i) = D[j](i);
      }
    }

    //////////////////////
    // Step 6. Sample P
    //////////////////////
    double shape1 = 0;
    double shape2 = 0;
    for (int j =0; j<m; ++j){
      shape1 =  std::abs(P0(j,j) + nstate[j] - 1);
      shape2 =  P0(j,j+1) + 1; //
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }
    P(m, m) = 1;

    /////////////////////////////////////////////////
    // STORE
    /////////////////////////////////////////////////
    if (iter >= burnin && ((iter % thin) == 0)) {
      Matrix<> tbeta = ::t(beta); //transpose beta for R output
      for (int i=0; i<(ns*K); ++i){
	beta_store(sampcount, i) = tbeta[i];// stored by the order of (11, 12, 13, 21, 22, 23)
      }
      for (int i=0; i<ns; ++i){
	sigma_store(sampcount, i) = sigma2[i];
      }
      Matrix<> tD_record = ::t(D_record);
      for (int i=0; i<(ns*Q*Q); ++i){
	D_store(sampcount,i) = tD_record(i);// stored by the order of (D11, D12, D13, 21, 22, 23)
      }
      for (int j=0; j<ns*ns; ++j){
	P_store(sampcount,j)= P[j];
      }
      s_store(sampcount,_) = state;
      for (int l=0; l<ntime ; ++l){
	ps_store(l,_) = ps_store(l,_) + ps(l,_);
      }
      ++sampcount;
    }


    /////////////////////////////////////////////////
    // REPORT
    /////////////////////////////////////////////////
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n ----------------------------------------------------------------------- \n");
      Rprintf("HMMpanelRE iteration %i of %i \n", (iter+1), tot_iter);
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

    Matrix<double> Dst = meanc(D_store);
    Matrix<>* Dinv_st = new Matrix<>[ns];
    Matrix<>* D_st = new Matrix<>[ns];
    int Dcount = 0;
    for (int j=0; j<ns; ++j){
      Matrix<> Dtemp(Q, Q);
      for (int k = 0; k<(Q*Q); ++k){
	Dtemp[k] = Dst[Dcount + k];
      }
      Dcount = Dcount + Q*Q;
      D_st[j] = Dtemp;
      Dinv_st[j] = invpd(Dtemp);
    }
    Matrix<> Dst_record(ns, Q*Q);// for record D in D_store
    for (int j = 0; j<ns; ++j){
      for (int i=0; i<(Q*Q); ++i){
	Dst_record(j,i) = D_st[j](i);
      }
    }

    Matrix<double> P_vec_st = meanc(P_store);
    const Matrix<double> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){
      P_st[j] = P_vec_st[j];
    }

    Matrix<> density_beta(nstore, ns);
    Matrix<> density_local_beta(ns, 1);
    //////////////////////////////////////////////////////////////////
    // 1. pdf.beta | D_g, sigma2_g, P_g, bi_g
    //////////////////////////////////////////////////////////////////
    Matrix<> Dmcmc(Q, Q);

    for (int iter = 0; iter<nstore; ++iter){
      int Dcount = 0;
      int beta_count = 0;
      Matrix<int> nstate(ns, 1);
      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<ntime; ++i){
	  if (s_store(iter, i) == j + 1) {
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
	beta_count = beta_count + nstate[j];

        for (int h = 0; h< (Q*Q); ++h){
	  Dmcmc[h] = D_store(iter, Dcount + h);
	}
	Dcount = Dcount + Q*Q;
	Matrix<> Dinv_j = invpd(Dmcmc);
	Matrix<> XVX(K, K);
	Matrix<> XVY(K, 1);
	for(int tt = (beta_count - nstate[j]); tt<beta_count; ++tt) {
	  XVX += (1/sigma_store(iter, j))*cpXt_arr[tt]-pow(1/sigma_store(iter, j),2)*tXWt_arr[tt]*invpd(Dinv_j + cpWt_arr[tt]*(1/sigma_store(iter, j)))*tWXt_arr[tt];
	  XVY += (1/sigma_store(iter, j))*tXYt_arr[tt]-pow(1/sigma_store(iter, j),2)*tXWt_arr[tt]*invpd(Dinv_j + cpWt_arr[tt]*(1/sigma_store(iter, j)))*tWYt_arr[tt];
	}
	Matrix<> beta_var = invpd(B0 + XVX/sigma_store(iter, j));
	Matrix<> beta_mean = beta_var*(B0*b0 + XVY/sigma_store(iter, j));
	if (K == 1){
	  density_beta(iter, j) = dnorm(beta_st(j), beta_mean[0], sqrt(beta_var[0]));
	}
	else{
	  density_beta(iter, j) = ::exp(lndmvn(::t(beta_st(j,_)), beta_mean, beta_var));
	}
      }
    }
    double pdf_beta = log(prod(meanc(density_beta)));

    //////////////////////////////////////////////////////////////////
    // 2. pdf.D
     //////////////////////////////////////////////////////////////////
    Matrix<> density_D(nstore, ns);
    for (int iter = 0; iter<nstore; ++iter){
      Matrix<> SSE(ns, 1);
      int beta_count = 0;
      Matrix<int> YN(ns, 1);
      Matrix<int> nstate(ns, 1);
      /////////////////////////////////////////////////
      // 2.1. s | beta_st, y, bi, sigma, D, P
      /////////////////////////////////////////////////
      Matrix<> F(ntime, ns);
      Matrix<> pr1(ns, 1);
      pr1[0] = 1;
      Matrix<> py(ns, 1);
      Matrix<> pstyt1(ns, 1);
      for (int tt=0; tt<ntime ; ++tt) {
	int nsubj_s = time_groupinfo(tt, 2);
	Matrix<> Wbi(nsubj_s, 1);
	for (int j=0; j<ns;++j){
	  Matrix<> Sigma = eye(nsubj_s)*sigma2[j];
	  Matrix<> WDW = Wt_arr[tt]*D[j]*tWt_arr[tt];
	  Matrix<> Mu = Xt_arr[tt]*::t(beta_st(j,_));
	  py[j] = ::exp(lndmvn(Yt_arr[tt], Mu, WDW+Sigma));
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
      ps(ntime-1, _) = F(ntime-1, _);
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
      // 2.2. P | beta_st, y, bi, sigma, s, D
      /////////////////////////////////////////////////
      for (int j = 0; j<ns; ++j){
      	for (int i = 0; i<ntime; ++i){
      	  if (state(i) == j + 1) { // since j starts from 0
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
      }
      double shape1 = 0;
      double shape2 = 0;
      for (int j =0; j<m; ++j){
	shape1 =  std::abs(P0(j,j) + nstate[j] - 1);
	shape2 =  P0(j,j+1) + 1; //
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
      }
      P(m, m) = 1;
      /////////////////////////////////////////////////
      // 2.3. bi | beta_st, y, D, sigma, s, P
      /////////////////////////////////////////////////
      beta_count = 0;
      for (int j = 0; j<ns; ++j){
	beta_count = beta_count + nstate[j];
	for(int s = 0; s<nsubj; ++ s) {
	  Matrix<> yj = Yk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  Matrix<> Xj = Xk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), K-1);
	  Matrix<> Wj = Wk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), Q-1);
	  Matrix<> b_var = invpd(Dinv[j]+ ::t(Wj)*Wj/sigma2[j]);
	  Matrix<> b_mean = b_var*::t(Wj)*(yj - Xj*::t(beta_st(j,_)))/sigma2[j];
	  bi[j](_,s) = stream.rmvnorm(b_mean, b_var);

	  // FOR SIGMA
	  YN[j] = YN[j] + yj.rows();
	  Matrix<> e = ::t(yj - Xj*::t(beta_st(j,_)) - Wj*bi[j](_,s))*(yj - Xj*::t(beta_st(j,_)) - Wj*bi[j](_,s));
	  SSE[j] = SSE[j] + e[0];
	}
      }
      /////////////////////////////////////////////////
      // 2.4. D | beta_st, y, bi, sigma, s, P
      /////////////////////////////////////////////////
      for (int j = 0; j<ns; ++j){
	Matrix<> SSB = bi[j]*t(bi[j]);
	Matrix<> D_scale = invpd(R0inv + SSB);
	int D_dof = r0 + nsubj;
	Dinv[j] = stream.rwish(D_dof, D_scale);
	density_D(iter, j) = ::exp(lndwish(Dinv_st[j], D_dof, D_scale));
	D[j] = inv(Dinv[j]);
      }
      /////////////////////////////////////////////////
      // 2.5. sigma | beta_st, y, bi, D, s, P
      /////////////////////////////////////////////////
      for (int j = 0; j<ns; ++j){
	double nu = (c0 + (double)YN[j])/2;
	double scale = (d0 + SSE[j])/2;
	sigma2[j] = stream.rigamma(nu, scale);
      }
    }// end of reduced run for step 2

    double pdf_D = log(prod(meanc(density_D)));


    //////////////////////////////////////////////////////////////////
    // 3. pdf.sigma2 : sigma2| beta_st, D_st, y, bi, s, P
    //////////////////////////////////////////////////////////////////
    Matrix<> density_sigma2(nstore, ns);
    for (int iter = 0; iter<nstore; ++iter){
      Matrix<> SSE(ns, 1);
      int beta_count = 0;
      Matrix<int> YN(ns, 1);
      Matrix<int> nstate(ns, 1);

      /////////////////////////////////////////////////
      // 3.1. s| beta_st, D_st, y, bi, sigma2, P
      /////////////////////////////////////////////////
      Matrix<> F(ntime, ns);
      Matrix<> pr1(ns, 1);
      pr1[0] = 1;
      Matrix<> py(ns, 1);
      Matrix<> pstyt1(ns, 1);

      for (int tt=0; tt<ntime ; ++tt) {
	int nsubj_s = time_groupinfo(tt, 2);
	Matrix<> Wbi(nsubj_s, 1);
	for (int j=0; j<ns;++j){
	  Matrix<> Sigma = eye(nsubj_s)*sigma2[j];
	  Matrix<> WDW = Wt_arr[tt]*D_st[j]*tWt_arr[tt];
	  Matrix<> Mu = Xt_arr[tt]*::t(beta_st(j,_));
	  py[j]  =  ::exp(lndmvn(Yt_arr[tt], Mu, WDW+Sigma));
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
      // 3.2. bi| beta_st, D_st, y, sigma2, s, P
      /////////////////////////////////////////////////
      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<ntime; ++i){
	  if (state(i) == j + 1) { // since j starts from 0
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
      }
      for (int j = 0; j<ns; ++j){
	beta_count = beta_count + nstate[j];
	for(int s = 0; s<nsubj; ++ s) {
	  Matrix<> yj = Yk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	  Matrix<> Xj = Xk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), K-1);
	  Matrix<> Wj = Wk_arr[s]((beta_count - nstate[j]), 0, (beta_count - 1), Q-1);
	  Matrix<> b_var = inv(Dinv_st[j]+ ::t(Wj)*Wj/sigma2[j]);
	  Matrix<> b_mean = b_var*::t(Wj)*(yj - Xj*::t(beta_st(j,_)))/sigma2[j];
	  bi[j](_,s) = stream.rmvnorm(b_mean, b_var);

	  // FOR SIGMA
	  YN[j] = YN[j] + yj.rows();
	  Matrix<> e = ::t(yj - Xj*::t(beta_st(j,_)) - Wj*bi[j](_,s))*(yj - Xj*::t(beta_st(j,_)) - Wj*bi[j](_,s));
	  SSE[j] = SSE[j] + e[0];
	}

	/////////////////////////////////////////////////
	// 3.3. sigma2| beta_st, D_st, y, bi, s, P
	/////////////////////////////////////////////////
	double nu = (c0 + (double)YN[j])/2;
	double scale = (d0 + SSE[j])/2;
	sigma2[j] = stream.rigamma(nu, scale);
	density_sigma2(iter, j) = ::exp(lndinvgamma(sigma2_st[j], nu, scale));
      }

      /////////////////////////////////////////////////
      // 3.4. P| beta_st, D_st, y, bi, sigma2, s
      /////////////////////////////////////////////////
      double shape1 = 0;
      double shape2 = 0;
      for (int j =0; j<m; ++j){
	shape1 =  std::abs(P0(j,j) + nstate[j] - 1);
	shape2 =  P0(j,j+1) + 1; //
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
      }
      P(m, m) = 1;
    }// end of reduced run
    double pdf_sigma2 = log(prod(meanc(density_sigma2)));

    //////////////////////
    // 4. pdf.P:
    //////////////////////
    Matrix<> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      Matrix<int> nstate(ns, 1);
      /////////////////////////////////////////////////
      // 4.1. s| y, P, beta_st, sigma_st, D_st
      /////////////////////////////////////////////////
      Matrix<> F(ntime, ns);
      Matrix<> pr1(ns, 1);
      pr1[0] = 1;
      Matrix<> py(ns, 1);
      Matrix<> pstyt1(ns, 1);

      for (int tt=0; tt<ntime ; ++tt) {
	int nsubj_s = time_groupinfo(tt, 2);
	Matrix<> Wbi(nsubj_s, 1);
	for (int j=0; j<ns;++j){
	  Matrix<> Sigma = eye(nsubj_s)*sigma2_st[j];
	  Matrix<> WDW = Wt_arr[tt]*D_st[j]*t(Wt_arr[tt]);
	  Matrix<> Mu = Xt_arr[tt]*::t(beta_st(j,_));
	  py[j]  =  ::exp(lndmvn(Yt_arr[tt], Mu, WDW+Sigma));
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

      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<ntime; ++i){
	  if (state(i) == j + 1) {
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
      }
      /////////////////////////////////////////////////
      // 4.2. P| y, s, beta_st, sigma_st, D_st
      /////////////////////////////////////////////////
      double shape1 = 0;
      double shape2 = 0;
      for (int j =0; j<m; ++j){
	shape1 =  std::abs(P0(j,j) + nstate[j] - 1);
	shape2 =  P0(j,j+1) + 1; //
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
	density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2);
      }
      P(m, m) = 1;
      density_P(iter, ns-1) = 1; //without this, there will be an error
    }// end of pdf.P MCMC run

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
    for (int tt=0; tt<ntime ; ++tt) {
      Matrix<> Wbi(nsubj, 1);
      for (int j=0; j<ns;++j){
	Matrix<> Sigma = eye(nsubj)*sigma2_st[j];
	Matrix<> WDW = Wt_arr[tt]*D_st[j]*t(Wt_arr[tt]);
	Matrix<> Mu = Xt_arr[tt]*::t(beta_st(j,_));
	py[j]  =  ::exp(lndmvn(Yt_arr[tt], Mu, WDW+Sigma));
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
    }

    loglike = sum(log(like));

    //////////////////////
    // log prior ordinate
    //////////////////////
    Matrix<double> density_beta_prior(ns, 1);
    Matrix<double> density_sigma2_prior(ns, 1);
    Matrix<double> density_D_prior(ns, 1);
    Matrix<double> density_P_prior(ns, 1);
    density_P_prior[ns-1] = 0; //

    for (int j=0; j<ns ; ++j){
      if (K == 1){
	density_beta_prior[j] = log(dnorm(beta_st(j), b0[0], sqrt(B0inv[0])));
      }
      else{
	density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv);
      }
      density_D_prior[j] = lndwish(Dinv_st[j], r0, R0);
      density_sigma2_prior[j] = lndinvgamma(sigma2_st[j], c0/2, d0/2);
    }
    for (int j=0; j<m ; ++j){
      density_P_prior[j] = log(dbeta(P_st(j,j), P0(j,j), P0(j,j+1)));
    }

    // compute marginal likelihood
    double logprior = sum(density_beta_prior) + sum(density_sigma2_prior) + sum(density_D_prior) +
      sum(density_P_prior);
    logmarglike = (loglike + logprior) - (pdf_beta + pdf_P + pdf_sigma2 + pdf_D);

    if (verbose >0){
      Rprintf("\n ----------------------------------------------------------------------- \n");
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("logprior = %10.5f\n", logprior);
      Rprintf("pdf_beta = %10.5f\n", pdf_beta);
      Rprintf("pdf_Sigma = %10.5f\n", pdf_sigma2);
      Rprintf("pdf_D = %10.5f\n", pdf_D);
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
  delete [] Wk_arr;
  delete [] Xt_arr;
  delete [] tXt_arr;
  delete [] Yt_arr;
  delete [] cpWt_arr;
  delete [] cpXt_arr;
  delete [] Wt_arr;
  delete [] tWt_arr;
  delete [] tXYt_arr;
  delete [] tXWt_arr;
  delete [] tWXt_arr;
  delete [] tWYt_arr;
  delete [] bi;
}// end of HMMpanelRE_impl

extern "C" {
  void cHMMpanelRE(double* betadata, const int* betarow, const int* betacol,
		  double* sigmadata, double* Ddata, double *psout, double *sout,
		  const int* nsubj, const int* ntime,
		  const int* m, const int* nobs, const int* subjectid, const int* timeid,
		  const double* Ydata, const int* Yrow, const int* Ycol,
		  const double* Xdata, const int* Xrow, const int* Xcol,
		  const double* Wdata, const int* Wrow, const int* Wcol,
		  const double* YTdata, const double* XTdata, const double* WTdata,
		  const int* burnin, const int* mcmc, const int* thin, const int* verbose,
		  const int *uselecuyer, const int *seedarray, const int *lecuyerstream,
		  const double* betastartdata, const double* sigma2start,
		  const double *Pstart,
		  const double* b0data, const double* B0data,
		  const double* c0, const double* d0, const int* r0, const double* R0data,
		  const double* subject_groupinfodata, const double* time_groupinfodata,
		  double *logmarglikeholder, double *loglikeholder,
		  const int *chib){
    // pull together Matrix objects
    Matrix<> Y(*Yrow, *Ycol, Ydata);
    Matrix<> X(*Xrow, *Xcol, Xdata);
    Matrix<> W(*Wrow, *Wcol, Wdata);
    Matrix<> YT(*Yrow, *Ycol, YTdata);
    Matrix<> XT(*Xrow, *Xcol, XTdata);
    Matrix<> WT(*Wrow, *Wcol, WTdata);
    Matrix<> betastart(*Xcol, 1, betastartdata);

    Matrix<> b0(*Xcol, 1, b0data);
    Matrix<> B0(*Xcol, *Xcol, B0data);
    Matrix<> R0(*Wcol, *Wcol, R0data);

    Matrix<> Dstart = invpd(R0);
    Matrix<int> subjectid_mat(*nobs, 1, subjectid);
    Matrix<int> timeid_mat(*nobs, 1, timeid);
    Matrix<> subject_groupinfo(*nsubj, 3, subject_groupinfodata);
    Matrix<> time_groupinfo(*ntime, 3, time_groupinfodata);
    const int mns = *m + 1;
    Matrix<> beta_store(*betarow, *betacol);
    Matrix<> sigma_store(*betarow, mns);
    Matrix<> D_store(*betarow, *Wcol* *Wcol * mns);
    double logmarglike = 0.0;
    double loglike = 0.0;

    if (*m == 0){
      MCMCPACK_PASSRNG2MODEL(GaussianPanelRE_impl,
			     *nsubj,  *ntime,  *nobs,
			     subjectid_mat, timeid_mat,
			     Y, X, W,
			     *burnin, *mcmc, *thin, *verbose,
			     *sigma2start, Dstart,
			     b0, B0, *c0, *d0, *r0, R0,
			     time_groupinfo,  subject_groupinfo,
			     beta_store, sigma_store, D_store,
			     logmarglike, loglike, *chib);


      // store marginal likelihood
      logmarglikeholder[0] = logmarglike;
      loglikeholder[0] = loglike;

      for (int i=0; i < (*betarow* *betacol); ++i){
	betadata[i] = beta_store(i);
      }
      for (int i=0; i < (*betarow); ++i){
	sigmadata[i] = sigma_store(i);
      }
      for (int i=0; i < (*betarow*mns* *Wcol* *Wcol); ++i){
	Ddata[i] = D_store(i);
      }
    }
    else {
      Matrix <> P(mns, mns, Pstart);
      Matrix<> P_store(*betarow,  mns*mns);
      Matrix<> s_store(*betarow,  *ntime*mns);
      Matrix<> ps_store(*ntime, mns);
      MCMCPACK_PASSRNG2MODEL(HMMpanelRE_impl,
			     *nsubj,  *ntime,  *m, *nobs,
			     subjectid_mat, timeid_mat,
			     Y, X, W, YT, XT, WT,
			     *burnin, *mcmc, *thin, *verbose,
			     betastart, *sigma2start, Dstart,
			     b0, B0, *c0, *d0, *r0, R0,
			     P,  time_groupinfo,  subject_groupinfo,
			     beta_store, sigma_store, D_store, P_store, ps_store, s_store,
			     logmarglike, loglike, *chib);
      // store marginal likelihood
      logmarglikeholder[0] = logmarglike;
      loglikeholder[0] = loglike;

      for (int i=0; i < (*betarow* *betacol); ++i){
	betadata[i] = beta_store(i);
      }
      for (int i=0; i < (*betarow*mns); ++i){
	sigmadata[i] = sigma_store(i);
      }
      for (int i=0; i < (*betarow*mns* *Wcol* *Wcol); ++i){
	Ddata[i] = D_store(i);
      }
      for (int i = 0; i<(*ntime *mns); ++i){
	psout[i] = ps_store[i];
      }
      for (int i = 0; i<(*betarow* *ntime *mns); ++i){
	sout[i] = s_store[i];
      }

    }
  }// end of cHMMpanelRE

}// end of extern C

#endif /* CHMMPANELRE_CC  */
