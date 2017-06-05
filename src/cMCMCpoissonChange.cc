//////////////////////////////////////////////////////////////////////////
// cMCMCpoissonRegChange.cc is C++ code to estimate a Poisson changepoint model
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr
//
// 07/06/2007
// 12/20/2007 included in MCMCpack
// 03/09/2012 fixed a bug in a random uniform draw (thanks to Matt Blackwell)
//////////////////////////////////////////////////////////////////////////

#ifndef CMCMCPOISSONCHANGE_CC
#define CMCMCPOISSONCHANGE_CC

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

#include <R.h>
#include <R_ext/Utils.h>

using namespace std;
using namespace scythe;

//tau and component sampler
template <typename RNGTYPE>
Matrix<> tau_comp_sampler(rng<RNGTYPE>& stream,
				const int m,
				const int totcomp,
				const Matrix<>& Y,
				const Matrix<>& X,
				const Matrix<>& wr,
				const Matrix<>& mr,
				const Matrix<>& sr,
				const Matrix<>& beta,
				const Matrix<>& s){

  // itialize
  const int n = Y.rows();
  const int k = X.cols();
  Matrix<int> component(totcomp, 1);
  Matrix<> tau(totcomp, 1);
  Matrix<> post_taut_mat(totcomp, 5);
  Matrix<> post_tau_mat(totcomp, 5);
  Matrix<> norm_post_tau(totcomp, 5);
  Matrix<> cumsum_post_tau(totcomp, 5);
  Matrix<> xt(1, k);
  Matrix<> wr_mat = eye(5);
  for (int i=0; i<5; ++i) {
    wr_mat(i,i) = wr[i];
  }

  int tau_start = 0;
  int tau_end = 0;

  for(int t=0; t<n; ++t){
    int yt = (int)Y[t];
    xt = X(t,_);
    tau_start = tau_end;
    tau_end = tau_start  + yt + 1;

    int st = (int)s[t];
    Matrix<> mu_t = exp(xt* ::t(beta(st-1,_)));
    double mut = mu_t[0];

    if (yt == 0){
      double tau_double = 1 + stream.rexp(mut);
      tau[tau_end - 1] = tau_double;
      for (int h=0; h<5; ++h){
	double first = 1/(sr[h]*tau_double);
	double second = (log(tau_double) + log(mut) - mr[h])/sr[h];
	post_taut_mat(tau_end-1, h) = first*exp(-0.5*second*second);
      }
    }

    else {
      Matrix<> ut = stream.runif(yt, 1);
      // redraw the uniform if there are any repeats
      // thanks to Matt Blackwell
      while (unique(ut).size() != ut.size()) {
	ut = stream.runif(yt, 1);
      }
      Matrix<> sort_ut = sort(ut);
      Matrix<> tau_tj(yt, 1);
      for(int i=1; i<yt; ++i){
	tau_tj[i] = sort_ut[i] - sort_ut[i-1];
      }
      tau_tj[0] = sort_ut[0];

      double sum_tau_tj = sum(tau_tj);
      double tau_last = 1 - sum_tau_tj + stream.rexp(mut);

      Matrix<> tau_mat(yt+1, 1);
      tau_mat(0, 0, yt-1, 0) = tau_tj(_,0);
      tau_mat[yt] = tau_last;

      for (int i = 0; i<(yt+1); ++i){
	tau[i + tau_start] = tau_mat[i];
	for (int h=0; h<5; ++h){
	  double first = 1/(sr[h]*tau_mat[i]);
	  double second = (log(tau_mat[i]) + log(mut) - mr[h])/sr[h];
	  post_taut_mat(i+tau_start, h) = first*exp(-0.5*second*second);
	}
      }
    }

  }
  post_tau_mat = post_taut_mat*wr_mat;

  for(int i = 0; i<totcomp ; ++i){
    norm_post_tau(i,_) = post_tau_mat(i,_)/sum(post_tau_mat(i,_));
    cumsum_post_tau(i,0) = norm_post_tau(i,0);
    for (int j=1; j<5; ++j){
      cumsum_post_tau(i,j) = cumsum_post_tau(i,j-1) + norm_post_tau(i,j);
    }
    double U = stream.runif();
    if (U < cumsum_post_tau(i,0)){
      component[i] = 1;
    }
    if (cumsum_post_tau(i,0)<=U&&U<cumsum_post_tau(i,1)){
      component[i] = 2;
    }
    if (cumsum_post_tau(i,1)<=U&&U<cumsum_post_tau(i,2)){
      component[i] = 3;
    }
    if (cumsum_post_tau(i,2)<=U&&U<cumsum_post_tau(i,3)){
      component[i] = 4;
    }
    if (cumsum_post_tau(i,3)<=U&&U<=1){
      component[i] = 5;
    }
  }

  Matrix<> TAUout(totcomp, 2);
  TAUout(_, 0) = tau(_, 0);
  TAUout(_, 1) = component(_, 0);

  return TAUout;

}

template <typename RNGTYPE>
Matrix<> poisson_state_sampler(rng<RNGTYPE>& stream,
			       const int& m,
			       const Matrix<>& Y,
			       const Matrix<>& lambda,
			       const Matrix<>& P){
  const int ns = m + 1;
  const int n = Y.rows();
  Matrix<> F(n, ns);
  Matrix<> pr1(ns, 1);
  pr1[0] = 1;
  Matrix<> py(ns, 1);
  Matrix<> pstyt1(ns, 1);

  for (int t=0; t<n ; ++t){
    int yt = (int) Y[t];
    for (int j=0; j<ns ; ++j){
      py[j]  =  dpois(yt, lambda[j]);
    }
    if (t==0) pstyt1 = pr1;
    else {
      pstyt1 =  ::t(F(t-1,_)*P);
    }
    Matrix<> unnorm_pstyt = pstyt1%py;
    Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
    for (int j=0; j<ns ; ++j){
      F(t,j) = pstyt(j);
    }
  }
  Matrix<int> s(n, 1);
  Matrix<> ps(n, ns);
  ps(n-1,_) = F(n-1,_);
  s(n-1) = ns;

  Matrix<> pstyn(ns, 1);
  double pone = 0.0;
  int t = n-2;
  while (t >= 0){
    int st = s(t+1);
    Matrix<> Pst_1 = ::t(P(_,st-1));
    Matrix<> unnorm_pstyn = F(t,_)%Pst_1;
    pstyn = unnorm_pstyn/sum(unnorm_pstyn);
    if (st==1)   s(t) = 1;
    else{
      pone = pstyn(st-2);
      if(stream.runif () < pone) s(t) = st-1;
      else s(t) = st;
    }
    ps(t,_) = pstyn;
    --t;
  }// end of while loop

  Matrix<> Sout(n, ns+1);
  Sout(_, 0) = s(_,0);
  for (int j = 0; j<ns; ++j){
    Sout(_,j+1) = ps(_, j);
  }

  return Sout;

}

template <typename RNGTYPE>
Matrix<> poisson_reg_state_sampler(rng<RNGTYPE>& stream,
				   const int m,
				   const Matrix<>& Y,
				   const Matrix<>& X,
				   const Matrix<>& beta,
				   const Matrix<>& P){

  const int ns = m + 1;
  const int n = Y.rows();

  Matrix<> F(n, ns);
  Matrix<> pr1(ns, 1);
  pr1[0] = 1;
  Matrix<> py(ns, 1);
  Matrix<> pstyt1(ns, 1);

  for (int t=0; t<n ; ++t){
    int yt = (int) Y[t];
    Matrix<> lambda = exp(X(t,_)*::t(beta));
    for (int j = 0; j< ns; ++j){
      py[j]  =  dpois(yt, lambda[j]);
    }
    if (t==0) pstyt1 = pr1;
    else {
      pstyt1 =  ::t(F(t-1,_)*P);
    }
    Matrix<> unnorm_pstyt = pstyt1%py;

    const Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
    for (int j=0; j<ns ; ++j) F(t,j) = pstyt(j);

  }

  Matrix<int> s(n, 1);
  Matrix<> ps = Matrix<>(n, ns);
  ps(n-1,_) = F(n-1,_);
  s(n-1) = ns;

  Matrix<> pstyn = Matrix<>(ns, 1);
  double pone = 0.0;
  int t = n-2;
  while (t >= 0){
    int st = s(t+1);
    Matrix<> Pst_1 = ::t(P(_,st-1));
    Matrix<> unnorm_pstyn = F(t,_)%Pst_1;
    pstyn = unnorm_pstyn/sum(unnorm_pstyn);
    if (st==1)   s(t) = 1;
    else{
      pone = pstyn(st-2);
      if(stream.runif() < pone) s(t) = st-1;
      else s(t) = st;
    }
    ps(t,_) = pstyn;
    --t;
  }

  Matrix<> Sout(n, ns+1);
  Sout(_, 0) = s(_,0);
  for (int j = 0; j<ns; ++j){
    Sout(_,j+1) = ps(_, j);
  }
  return Sout;
}


////////////////////////////////////////////
// cMCMCpoissonChangepoint implementation.
////////////////////////////////////////////
template <typename RNGTYPE>
void MCMCpoissonChangepoint_impl(rng<RNGTYPE>& stream,
				 double *betaout,
				 double *Pout,
				 double *psout,
				 double *sout,
				 const double *Ydata,
				 const int *Yrow,
				 const int *Ycol,
				 const int *m,
				 const double *c0,
				 const double *d0,
				 const int *burnin,
				 const int *mcmc,
				 const int *thin,
				 const int *verbose,
				 const double *betastart,
				 const double *Pstart,
				 const double *a,
				 const double *b,
				 const double *A0data,
				 double *logmarglikeholder,
				 double *loglikeholder,
				 const int *chib)

{
  const Matrix <> Y(*Yrow, *Ycol, Ydata);
  const int tot_iter = *burnin + *mcmc;
  const int nstore = *mcmc / *thin;
  const int n = Y.rows();
  const int ns = *m + 1;
  const Matrix <> A0(ns, ns, A0data);

  Matrix <> lambda(ns, 1, betastart);
  Matrix <> P(ns, ns, Pstart);

  Matrix<> lambda_store(nstore, ns);
  Matrix<> P_store(nstore, ns*ns);
  Matrix<> ps_store(n, ns);
  Matrix<> s_store(nstore, n);

  //MCMC loop
   int count = 0;
  for (int iter = 0; iter < tot_iter; ++iter){

    //////////////////////
    // 1. Sample s
    //////////////////////
    Matrix<> Sout = poisson_state_sampler(stream, *m, Y, lambda, P);
    Matrix<> s = Sout(_, 0);
    Matrix<> ps(n, ns);
    for (int j = 0; j<ns; ++j){
      ps(_,j) = Sout(_,j+1);
    }

    //////////////////////
    // 2. Sample lambda
    //////////////////////
    Matrix<> addY(ns, 1);
    Matrix<> addN(ns, 1);

    for (int j = 0; j<ns; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == (j+1)) {
	  addN[j] = addN[j] + 1;
	  addY[j] = addY[j] + Y[i];
	}
      }
      double c1 = addY[j] + *c0;
      double d1 = addN[j] + 1/ *d0;
      lambda[j] = stream.rgamma(c1, d1);
    }

    //////////////////////
    // 3. Sample P
    //////////////////////
    double shape1 = 0;
    double shape2 = 0;
    P(ns-1, ns-1) = 1;

    for (int j =0; j< (ns-1); ++j){
      shape1 =  A0(j,j) + addN[j] - 1;
      shape2 =  A0(j,j+1) + 1;
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }

    if (iter >= *burnin && ((iter % *thin)==0)){
      for (int i=0; i<ns; ++i)
	lambda_store(count,i) = lambda[i];
      for (int j=0; j<ns*ns; ++j)
	P_store(count,j)= P[j];
      s_store(count,_) = s;
      for (int l=0; l<n ; ++l)
	ps_store(l,_) = ps_store(l,_) + ps(l,_);

      ++count;

    }

    if(*verbose > 0 && iter % *verbose == 0){
      Rprintf("\n\n MCMCpoissonChange iteration %i of %i", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
	Rprintf("\n The number of observations in state %i is %10.5f", j + 1, addN[j]);
      }
      for (int j=0; j<ns; ++j)
        Rprintf("\n lambda in state %i is %10.5f", j + 1, lambda[j]);
    }

    R_CheckUserInterrupt();

  }

  if(*chib ==1){

    Matrix<> lambda_st = meanc(lambda_store);
    Matrix<> P_vec_st = meanc(P_store);
    const Matrix<> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){
      P_st[j] = P_vec_st[j];
    }

    //////////////////////
    // lambda
    //////////////////////
    Matrix<> density_lambda(nstore, ns);
    for (int iter = 0; iter<nstore; ++iter){

      Matrix<> addY(ns, 1);
      Matrix<> addN(ns, 1);

      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<n; ++i){
            if (s_store(iter, i) == (j+1)) {
	      addN[j] = addN[j] + 1;
	      addY[j] = addY[j] + Y[i];
	    }
	}
        double c1 = addY[j] + *c0;
        double d1 = addN[j] + *d0;
        density_lambda(iter, j) = dgamma(lambda_st[j], c1, 1/d1);
      }

    }
    double pdf_lambda = log(prod(meanc(density_lambda)));

    //////////////////////
    // P
    //////////////////////
    Matrix<> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      Matrix<> Sout = poisson_state_sampler(stream, *m, Y, lambda_st, P);
      Matrix <> s = Sout(_, 0);
      Matrix <> ps(n, ns);
      for (int j = 0; j<ns; ++j){
	ps(_,j) = Sout(_,j+1);
      }

      double shape1 = 0;
      double shape2 = 0;
      P(ns-1, ns-1) = 1;
      Matrix <> P_addN(ns, 1);
      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<n; ++i){
	  if (s[i] == (j+1)) {
	    P_addN[j] = P_addN[j] + 1;
	  }
	}
      }

      for (int j =0; j< (ns-1); ++j){
	shape1 =  A0(j,j) + P_addN[j] - 1;
	shape2 =  A0(j,j+1) + 1;
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
	density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2);
      }
      density_P(iter, ns-1) = 1;

    }
    double pdf_P = log(prod(meanc(density_P)));

    //////////////////////
    // likelihood
    //////////////////////
    Matrix<> F(n, ns);
    Matrix<> like(n, 1);
    Matrix<> pr1(ns, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);

    for (int t=0; t<n ; ++t){
      int yt = (int) Y[t];
      for (int j=0; j<ns ; ++j){
	py[j]  =  dpois(yt, lambda_st[j]);
      }
      if (t==0) pstyt1 = pr1;
      else {
	pstyt1 =  ::t(F(t-1,_)*P_st);
      }
      Matrix<> unnorm_pstyt = pstyt1%py;
      Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
      for (int j=0; j<ns ; ++j){
        F(t,j) = pstyt(j);
      }
      like[t] = sum(unnorm_pstyt);
    }

    double loglike = sum(log(like));

    //////////////////////
    // prior
    //////////////////////
    Matrix<> density_lambda_prior(ns, 1);
    Matrix<> density_P_prior(ns, 1);
    density_P[ns-1] = 1; //

    for (int j=0; j<ns ; ++j){
      density_lambda_prior[j] = log(dgamma(lambda_st[j], *c0, *d0));
    }

    for (int j =0; j< (ns-1); ++j){
      density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1)));
    }

    // compute marginal likelihood
    const double logprior = sum(density_lambda_prior) + sum(density_P_prior);
    const double logmarglike = (loglike + logprior) - (pdf_lambda + pdf_P);

    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;
  }

  R_CheckUserInterrupt();

  for (int i = 0; i<(nstore*ns); ++i){
    betaout[i] = lambda_store[i];
  }
  for (int i = 0; i<(nstore*ns*ns); ++i){
    Pout[i] = P_store[i];
  }
  for (int i = 0; i<(n*ns); ++i){
    psout[i] = ps_store[i];
  }
  for (int i = 0; i<(nstore*n); ++i){
    sout[i] = s_store[i];
  }


}

////////////////////////////////////////////
// cMCMCpoissonRegChange implementation.
////////////////////////////////////////////
template <typename RNGTYPE>
void MCMCpoissonRegChangepoint_impl(rng<RNGTYPE>& stream,
				    double *betaout,
				    double *Pout,
				    double *psout,
				    double *sout,
				    const double *Ydata,
				    const int *Yrow,
				    const int *Ycol,
				    const double *Xdata,
				    const int *Xrow,
				    const int *Xcol,
				    const int *m,
				    const int *burnin,
				    const int *mcmc,
				    const int *thin,
				    const int *verbose,
				    const double *betastart,
				    const double *Pstart,
				    const double *taustart,
				    const double *componentstart,
				    const double *a,
				    const double *b,
				    const double *b0data,
				    const double *B0data,
				    const double *A0data,
				    double *logmarglikeholder,
				    double *loglikeholder,
				    const double* wrin,
				    const double* mrin,
				    const double* srin,
				    const int *chib)
{
  const Matrix <> Y(*Yrow, *Ycol, Ydata);
  const Matrix <> X(*Xrow, *Xcol, Xdata);

  const int tot_iter = *burnin + *mcmc;
  const int nstore = *mcmc / *thin;
  const int n = Y.rows();
  const int k = X.cols();
  const int ns = *m + 1;
  const int totcomp = n + (int) sum(Y);
  const Matrix <> b0(k, 1, b0data);
  const Matrix <> B0(k, k, B0data);
  const Matrix <> B0inv = invpd(B0);
  Matrix <> wr(5, 1, wrin);
  Matrix <> mr(5, 1, mrin);
  Matrix <> sr(5, 1, srin);
  const Matrix <> A0(ns, ns, A0data);

  Matrix <> beta(ns, k, betastart);
  Matrix <> tau(totcomp, 1, taustart);
  Matrix <> component(totcomp, 1, componentstart);
  Matrix <> P(ns, ns, Pstart);

  Matrix<> beta_store(nstore, ns*k);
  Matrix<> P_store(nstore, ns*ns);
  Matrix<> ps_store(n, ns);
  Matrix<> s_store(nstore, n);
  Matrix<> component_store(nstore, totcomp);
  Matrix<> tau_store(nstore, totcomp);

  Matrix<> y_tilde(n ,1);
  Matrix<> Sigma_inv_sum(n, 1);

  //MCMC loop
  int count = 0;
  for (int iter = 0; iter < tot_iter; ++iter){

    int y_count = 0;
    for (int t = 0; t<n ; ++t){
      int yt = (int) Y[t];
      y_count = y_count + yt + 1;

      Matrix<> Yt_over_Sigma(yt + 1, 1);
      Matrix<> Sigma_inv(yt + 1, 1);

      for(int j = (y_count - yt - 1); j< y_count; ++j){
	int jone = (int) component[j] - 1 ; //zero base in C!
	Sigma_inv[j-(y_count-yt-1)] = 1/(sr[jone]*sr[jone]);
	Yt_over_Sigma[j-(y_count-yt-1)] = (log(tau[j])- mr[jone])*Sigma_inv[j-(y_count-yt-1)];
      }
      y_tilde[t] = sum(Yt_over_Sigma);
      Sigma_inv_sum[t] = sum(Sigma_inv);
    }
    //////////////////////
    // 1. Sample s
    //////////////////////
    Matrix <> Sout = poisson_reg_state_sampler(stream, *m, Y, X, beta, P);
    Matrix <> s = Sout(_, 0);
    Matrix <> ps(n, ns);
    for (int j=0; j<ns; ++j){
      ps(_,j) = Sout(_,j+1);
    }

    //////////////////////
    // 2. Sample beta
    //////////////////////
    int beta_count = 0;
    Matrix<int> nstate(ns, 1);

    for (int j = 0; j <ns ; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == (j+1)) {
	  nstate[j] = nstate[j] + 1;
	}
      }
      beta_count = beta_count + nstate[j];
      Matrix<> yj = y_tilde((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      Matrix<> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
      Matrix<> wi = Sigma_inv_sum((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      Matrix<> Xwj(nstate[j], k);
      for (int h = 0; h<nstate[j]; ++h){
	Xwj(h, _) = Xj(h,_)*wi[h];
      }
      Matrix<> Bn = invpd(B0 + ::t(Xj)*Xwj);
      Matrix<> bn = Bn*gaxpy(B0, b0,  -1*::t(Xj)*yj);
      beta(j,_) = stream.rmvnorm(bn, Bn);
   }

    //////////////////////
    // 3. Sample P
    //////////////////////
    double shape1 = 0;
    double shape2 = 0;
    P(ns-1, ns-1) = 1;

    for (int j =0; j< (ns-1); ++j){
      shape1 =  A0(j,j) + nstate[j] - 1;
      shape2 =  A0(j,j+1) + 1;
      P(j,j) =  stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }

    //////////////////////
    // 4. Sample tau
    //////////////////////
    Matrix <> TAUout = tau_comp_sampler(stream, *m, totcomp, Y, X, wr, mr, sr, beta, s);
    tau = TAUout(_, 0);
    component = TAUout(_, 1);

    if (iter >= *burnin && ((iter % *thin)==0)){
      Matrix<> tbeta = ::t(beta);
      for (int i=0; i<(ns*k); ++i){
	beta_store(count,i) = tbeta[i];
      }
      for (int j=0; j<ns*ns; ++j){
	P_store(count,j)= P[j];
      }
      for (int l=0; l<n ; ++l){
	ps_store(l,_) = ps_store(l,_) + ps(l,_);
      }
      s_store(count,_) = s(_, 0);
      tau_store(count,_) = tau(_, 0);
      component_store(count,_) = component(_, 0);

      ++count;

    }

    if(*verbose > 0 && iter % *verbose == 0){
      Rprintf("\n\n MCMCpoissonChange iteration %i of %i \n", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
	Rprintf("\n The number of observations in state %i is %10.5f", j+1, static_cast<double>(nstate[j]));
      }
      for (int i = 0; i<ns; ++i){
	for (int j = 0; j<k; ++j){
	  Rprintf("\n beta(%i) in state %i is %10.5f", j+1, i+1, beta(i, j));
	}
      }
    }

  }// end MCMC loop


  //////////////////////////////////////
  if(*chib ==1){
    //////////////////////////////////////

    Matrix<> betast = meanc(beta_store);
    Matrix<double, Row> beta_st(ns, k);
    for (int j = 0; j< ns*k; ++j){
      beta_st[j] = betast[j];
    }

    Matrix<> P_vec_st = meanc(P_store);
    const Matrix<> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){
      P_st[j] = P_vec_st[j];
    }

    //////////////////////
    // beta
    //////////////////////
    Matrix<> density_beta(nstore, ns);

    for (int iter = 0; iter<nstore; ++iter){
      int y_count = 0;
      for (int t = 0; t<n ; ++t){
	int yt = (int) Y[t];
	y_count = y_count + yt + 1;
	Matrix<> Yt_over_Sigma(yt + 1, 1);
	Matrix<> Sigma_inv(yt + 1, 1);
	for(int j = (y_count - yt - 1); j< y_count; ++j){
	  int jone = (int)component_store(iter, j) - 1 ;
	  Sigma_inv[j-(y_count-yt-1)] = 1/(sr[jone]*sr[jone]);
	  Yt_over_Sigma[j-(y_count-yt-1)] = (log(tau_store(iter, j))- mr[jone])*Sigma_inv[j-(y_count-yt-1)];
	}
	y_tilde[t] = sum(Yt_over_Sigma);
	Sigma_inv_sum[t] = sum(Sigma_inv);
      }

     int beta_count = 0;
      Matrix<int> nstate(ns, 1);

      for (int j = 0; j <ns ; ++j){
	for (int i = 0; i<n; ++i){
	  if (s_store(iter, i) == (j+1)) {
	    nstate[j] = nstate[j] + 1;
	  }
	}
	beta_count = beta_count + nstate[j];
	Matrix<> yj = y_tilde((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	Matrix<> wi = Sigma_inv_sum((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> Xwj(nstate[j], k);
	for (int h = 0; h<nstate[j]; ++h){
	  Xwj(h, _) = Xj(h,_)*wi[h];
	}
	  Matrix<> Bn = invpd(B0 + ::t(Xj)*Xwj);
	  Matrix<> bn = Bn*gaxpy(B0, b0,  -1*::t(Xj)*yj);
	  density_beta(iter, j) = exp(lndmvn(::t(beta_st(j,_)), bn, Bn));
      }
    }

    double pdf_beta = log(prod(meanc(density_beta)));

    //////////////////////
    // P
    //////////////////////
    Matrix<> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      Matrix <> Sout = poisson_reg_state_sampler(stream, *m, Y, X, beta_st, P);
      Matrix <> s = Sout(_, 0);
      Matrix <> ps(n, ns);
      for (int j = 0; j<ns; ++j){
	ps(_,j) = Sout(_,j+1);
      }
      double shape1 = 0;
      double shape2 = 0;
      P(ns-1, ns-1) = 1;

      Matrix <> P_addN(ns, 1);
      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<n; ++i){
	  if (s[i] == (j+1)) {
	    P_addN[j] = P_addN[j] + 1;
	  }
	}
      }
      for (int j =0; j< (ns-1); ++j){
	shape1 =  A0(j,j) + P_addN[j] - 1;
	shape2 =  A0(j,j+1) + 1;
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
	density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2);
      }
      density_P(iter, ns-1) = 1;

    }
    double pdf_P = log(prod(meanc(density_P)));

    //////////////////////
    // likelihood
    //////////////////////
    Matrix<> F = Matrix<>(n, ns);
    Matrix<> like(n, 1);
    Matrix<> pr1 = Matrix<>(ns, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);

    for (int t=0; t<n ; ++t){
      int yt = (int) Y[t];
      Matrix<> lambda = exp(X(t,_)*::t(beta_st));
      for (int j = 0; j< ns; ++j){
	py[j]  =  dpois(yt, lambda[j]);
      }
      if (t==0) pstyt1 = pr1;
      else {
	pstyt1 =  ::t(F(t-1,_)*P_st);
      }
      Matrix<> unnorm_pstyt = pstyt1%py;
      Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
      for (int j=0; j<ns ; ++j){
	F(t,j) = pstyt(j);
      }
      like[t] = sum(unnorm_pstyt);
    }

    const double loglike = sum(log(like));

    //////////////////////
    // prior
    //////////////////////
    Matrix<> density_beta_prior(ns, 1);
    Matrix<> density_P_prior(ns, 1);
    density_P[ns-1] = 1; //

    for (int j=0; j<ns ; ++j){
      density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv);
    }

    for (int j =0; j< (ns-1); ++j){
      density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1)));
    }

    // compute marginal likelihood
    const double logprior = sum(density_beta_prior) + sum(density_P_prior);
    const double logmarglike = (loglike + logprior) - (pdf_beta + pdf_P);
    if (*verbose > 0){
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("log_prior = %10.5f\n", logprior);
      Rprintf("log_beta = %10.5f\n", pdf_beta);
      Rprintf("log_P = %10.5f\n", pdf_P);
    }
    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;

  }

  R_CheckUserInterrupt();

  for (int i = 0; i<(nstore*ns*k); ++i){
    betaout[i] = beta_store[i];
  }
  for (int i = 0; i<(nstore*ns*ns); ++i){
    Pout[i] = P_store[i];
  }
  for (int i = 0; i<(n*ns); ++i){
    psout[i] = ps_store[i];
  }
  for (int i = 0; i<(nstore*n); ++i){
    sout[i] = s_store[i];
  }

}

extern "C" {
  void cMCMCpoissonChange( double *betaout,
			  double *Pout,
			  double *psout,
			  double *sout,
			  const double *Ydata,
			  const int *Yrow,
			  const int *Ycol,
			  const double *Xdata,
			  const int *Xrow,
			  const int *Xcol,
			  // const double *logoffset
			  const int *m,
			  const int *burnin,
			  const int *mcmc,
			  const int *thin,
			  const int *verbose,
			  const double *betastart,
			  const double *Pstart,
			  const double *taustart,
			  const double *componentstart,
			  const double *a,
			  const double *b,
			  const double *c0,
			  const double *d0,
			  const int* uselecuyer,
			  const int* seedarray,
			  const int* lecuyerstream,
			  const double *b0data,
			  const double *B0data,
			  const double *A0data,
			  double *logmarglikeholder,
			  double *loglikeholder,
			  const double *wrin,
			  const double *mrin,
			  const double *srin,
			  const int *chib){
    if(*Xcol == 1){
      MCMCPACK_PASSRNG2MODEL(MCMCpoissonChangepoint_impl,
			     betaout, Pout, psout, sout,
			     Ydata, Yrow, Ycol,
			     m, c0, d0,
			     burnin, mcmc, thin, verbose,
			     betastart, Pstart,
			     a, b, A0data,
			     logmarglikeholder, loglikeholder, chib)
	}
    else{
      MCMCPACK_PASSRNG2MODEL(MCMCpoissonRegChangepoint_impl,
			     betaout, Pout, psout, sout,
			     Ydata, Yrow, Ycol,
			     Xdata, Xrow, Xcol,
			     m, burnin, mcmc, thin, verbose,
			     betastart, Pstart,
			     taustart, componentstart,
			     a, b, b0data, B0data, A0data,
			     logmarglikeholder, loglikeholder,
			     wrin, mrin, srin, chib);
    }
  } // end MCMC

} // end extern "C"


#endif


