////////////////////////////////////////////////////////////////////
// cMCMCbinaryChange.cc is C++ code to estimate a binary changepoint model
// with a beta prior
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr
//
// 03/03/2009 Written
////////////////////////////////////////////////////////////////////

#ifndef CMCMCBINARYCHANGE_CC
#define CMCMCBINARYCHANGE_CC


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


template <typename RNGTYPE>
// bianry state sampler
Matrix<> binary_state_sampler(rng<RNGTYPE>& stream,
				    const int m,
				    const Matrix<>& Y,
				    const Matrix<>& theta,
				    const Matrix<>& P){

  const int ns = m + 1;
  const int n = Y.rows();
  Matrix<> F = Matrix<>(n, ns);
  Matrix<> pr1 = Matrix<>(ns, 1);
  pr1[0] = 1;
  Matrix<> py(ns, 1);
  Matrix<> pstyt1(ns, 1);

  for (int t=0; t<n ; ++t){
    int yt = (int) Y[t];
    for (int j=0; j<ns ; ++j){
      py[j]  =  dbinom(yt, 1, theta[j]);
    }
    if (t==0) pstyt1 = pr1;
    else {
      pstyt1 = ::t(F(t-1,_)*P);
    }
    Matrix<> unnorm_pstyt = pstyt1%py;
    const Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
    for (int j=0; j<ns ; ++j){
      F(t,j) = pstyt(j);
    }
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
// cMCMCbinaryChangepoint implementation.
////////////////////////////////////////////
template <typename RNGTYPE>
void MCMCbinaryChange_impl(rng<RNGTYPE>& stream, const Matrix<>& Y,
			   Matrix<>& phi, Matrix<>& P, const Matrix<>& A0,
			   const double m, const double c0, const double d0,
			   int burnin,  int mcmc,  int thin,
			   int verbose, bool chib,
			   Matrix<>& phi_store, Matrix<>& P_store,
			   Matrix<>& ps_store, Matrix<>& s_store,
			   double& logmarglike)
{
  const int tot_iter = burnin + mcmc;
  const int nstore = mcmc / thin;
  const int n = Y.rows();
  const int ns = m + 1;

  //MCMC loop
   int count = 0;
  for (int iter = 0; iter < tot_iter; ++iter){

    //////////////////////
    // 1. Sample s
    //////////////////////
    Matrix<> Sout = binary_state_sampler(stream, m, Y, phi, P);
    Matrix<> s = Sout(_, 0);
    Matrix<> ps(n, ns);
    for (int j = 0; j<ns; ++j){
      ps(_,j) = Sout(_,j+1);
    }

    //////////////////////
    // 2. Sample phi
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
      double c1 = addY[j] + c0;
      double d1 = addN[j] - addY[j] + d0;
      phi[j] = stream.rbeta(c1, d1);
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

    if (iter >= burnin && ((iter % thin)==0)){
      for (int i=0; i<ns; ++i)
	phi_store(count,i) = phi[i];
      for (int j=0; j<ns*ns; ++j)
	P_store(count,j)= P[j];
      s_store(count,_) = s;
      for (int l=0; l<n ; ++l)
	ps_store(l,_) = ps_store(l,_) + ps(l,_);

      ++count;

    }
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n\n MCMCbinaryChange iteration %i of %i", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
	Rprintf("\n The number of observations in state %i is %10.5f", j + 1, addN[j]);
      }
      for (int j=0; j<ns; ++j)
        Rprintf("\n phi in state %i is %10.5f", j + 1, phi[j]);
    }

    R_CheckUserInterrupt();

  }// end MCMC loop

  if(chib ==1){

    Matrix<> phi_st = meanc(phi_store);
    Matrix<> P_vec_st = meanc(P_store);
    const Matrix<> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){
      P_st[j] = P_vec_st[j];
    }

    //////////////////////
    // phi
    //////////////////////
    Matrix<> density_phi(nstore, ns);
    // Matrix<> density_phi_j(ns, 1);
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
        double c1 = addY[j] + c0;
        double d1 = addN[j] - addY[j] + d0;
        density_phi(iter, j) = dbeta(phi_st[j], c1, d1);
      }
    }
    double pdf_phi = log(prod(meanc(density_phi)));

    //////////////////////
    // P
    //////////////////////
    Matrix<> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      Matrix<> Sout = binary_state_sampler(stream, m, Y, phi_st, P);
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
	py[j]  =  dbinom(yt, 1, phi_st[j]);
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
    Matrix<> density_phi_prior(ns, 1);
    Matrix<> density_P_prior(ns, 1);

    for (int j=0; j<ns ; ++j){
      density_phi_prior[j] = scythe::lndbeta1(phi_st[j], c0, d0);
    }

    for (int j =0; j< (ns-1); ++j){
      density_P_prior[j] = scythe::lndbeta1(P_st(j,j), A0(j,j), A0(j,j+1));
    }

    // compute marginal likelihood
    double logprior = sum(density_phi_prior) + sum(density_P_prior);
    logmarglike = (loglike + logprior) - (pdf_phi + pdf_P);
    if(verbose > 0){
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("log_prior = %10.5f\n", logprior);
      Rprintf("log_phi = %10.5f\n", pdf_phi);
      Rprintf("log_P = %10.5f\n", pdf_P);
    }

  }
}
////////////////////////////////////////////
// Start cMCMCbinaryChangepoint function
///////////////////////////////////////////
extern "C"{
  void cMCMCbinaryChange(double *phiout, double *Pout, double *psout, double *sout,
			const double *Ydata, const int *Yrow, const int *Ycol, const int *m,
			const int *burnin, const int *mcmc, const int *thin, const int *verbose,
			const int *uselecuyer, const int *seedarray, const int *lecuyerstream,
			const double *phistart, const double *Pstart,
			const double *a, const double *b, const double *c0, const double *d0,
			const double *A0data, double *logmarglikeholder,
			const int *chib){

    // pull together Matrix objects
    const Matrix <> Y(*Yrow, *Ycol, Ydata);
    const int nstore = *mcmc / *thin;
    const int n = Y.rows();
    const int ns = *m + 1;

    // generate starting values
    Matrix <> phi(ns, 1, phistart);
    const Matrix <> A0(ns, ns, A0data);
    Matrix <> P(ns, ns, Pstart);

    double logmarglike;

    // storage matrices
    Matrix<> phi_store(nstore, ns);
    Matrix<> P_store(nstore, ns*ns);
    Matrix<> ps_store(n, ns);
    Matrix<> s_store(nstore, n);

    MCMCPACK_PASSRNG2MODEL(MCMCbinaryChange_impl, Y,
			   phi, P, A0, *m, *c0, *d0,
			   *burnin, *mcmc, *thin, *verbose, *chib,
			   phi_store, P_store, ps_store, s_store,
			   logmarglike)

      logmarglikeholder[0] = logmarglike;

    // return output
    for (int i = 0; i<(nstore*ns); ++i){
      phiout[i] = phi_store[i];
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
}

#endif

