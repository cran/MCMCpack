////////////////////////////////////////////////////////////////////
// cMCMCprobitChange.cc is C++ code to estimate a probit changepoint model
// with a beta prior
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr
//
// Written 07/06/2007
// Modifieed 11/02/2009
// Included 1/21/2011
//////////////////////////////////////////////////////////////////////////

#ifndef CMCMCPROBITCHANGE_CC
#define CMCMCPROBITCHANGE_CC

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

// probit state sampler
template <typename RNGTYPE>

Matrix<> probit_state_sampler(rng<RNGTYPE>& stream,
			      const int m,
			      const Matrix<double>& Y,
			      const Matrix<double>& X,
			      const Matrix<double>& beta,
			      const Matrix<double>& P){

  const int ns = m + 1;
  const int n = Y.rows();
  Matrix<double> F = Matrix<double>(n, ns);
  Matrix<double> pr1 = Matrix<double>(ns, 1);
  pr1[0] = 1;
  Matrix<double> py(ns, 1);
  Matrix<double> pstyt1(ns, 1);

  for (int t=0; t<n ; ++t){
    int yt = (int) Y[t];
    Matrix<double> mu = X(t,_)*::t(beta);
    for (int j=0; j<ns ; ++j){
      py[j]  =  dbinom(yt, 1, pnorm(mu[j], 0, 1));
    }
    if (t==0) pstyt1 = pr1;
    else {
      pstyt1 = ::t(F(t-1,_)*P);
    }
    Matrix<double> unnorm_pstyt = pstyt1%py;
    const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
    for (int j=0; j<ns ; ++j){
      F(t,j) = pstyt(j);
    }
  }

  Matrix<int> s(n, 1);
  Matrix<double> ps = Matrix<double>(n, ns);
  ps(n-1,_) = F(n-1,_);
  s(n-1) = ns;

  Matrix<double> pstyn = Matrix<double>(ns, 1);
  double pone = 0.0;
  int t = n-2;
  while (t >= 0){
    int st = s(t+1);
    Matrix<double> Pst_1 = ::t(P(_,st-1));
    Matrix<double> unnorm_pstyn = F(t,_)%Pst_1;
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

  Matrix<double> Sout(n, ns+1);
  Sout(_, 0) = s(_,0);
  for (int j = 0; j<ns; ++j){
    Sout(_,j+1) = ps(_, j);
  }

  return Sout;

}

////////////////////////////////////////////
// cMCMCprobitChangepoint implementation.
////////////////////////////////////////////
template <typename RNGTYPE>
void MCMCprobitChange_impl(rng<RNGTYPE>& stream,
			   const int m,
			   const Matrix<>& Y, const Matrix<>& X,
			   Matrix<>& beta, Matrix<>& P,
			   Matrix<>& b0, Matrix<>& B0,
			   const Matrix<>& A0,
			    int burnin,  int mcmc,  int thin,
			    int verbose, bool chib,
			   Matrix<>& beta_store, Matrix<>& Z_store,
			   Matrix<>& P_store, Matrix<>& ps_store, Matrix<int>& s_store,
			   double& logmarglike, double& loglike)
{
  const int tot_iter = burnin + mcmc;
  const int nstore = mcmc / thin;
  const int n = Y.rows();
  const int ns = m + 1;
  const int k = X.cols();
  const Matrix<> B0inv = invpd(B0);
  Matrix<> Z(n, 1);

   int count = 0;
  for (int iter = 0; iter < tot_iter; ++iter){

    // 1. Sample s
    Matrix<> Sout = probit_state_sampler(stream, m, Y, X, beta, P);
    Matrix<int> s = Sout(_, 0);
    Matrix <double> ps(n, ns);
    for (int j = 0; j<ns; ++j){
      ps(_,j) = Sout(_,j+1);
    }

    // 2. Sample Z
    for ( int i = 0; i < n; ++i) {
      const Matrix<> mu = X(i,_)*t(beta(s[i]-1,_));
      double muj = mu[0];
      if(muj>200){
	muj = 200;
      }
      if(muj<-200){
	muj = -200;
      }
      if (Y[i] == 1.0)
	Z[i] = stream.rtbnorm_combo(muj, 1.0, 0);
      if (Y[i] == 0.0)
	Z[i] = stream.rtanorm_combo(muj, 1.0, 0);
    }

    // 3. Sample beta
    int beta_count = 0;
    Matrix<int> nstate(ns, 1);
    for (int j = 0; j <ns ; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == (j+1)) {
	  nstate[j] = nstate[j] + 1;
	}
      }
      beta_count = beta_count + nstate[j];
      const Matrix<double> Zj = Z((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      const Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
      const Matrix<double> XpX = t(Xj)*Xj;
      const Matrix<double> XpZ = t(Xj)*Zj;
      beta(j,_) = NormNormregress_beta_draw(XpX, XpZ, b0, B0, 1.0, stream);
    }

     // 4. Sample P
    double shape1 = 0;
    double shape2 = 0;
    P(ns-1, ns-1) = 1;

    for (int j =0; j< (ns-1); ++j){
      shape1 =  A0(j,j) + nstate[j] - 1;
      shape2 =  A0(j,j+1) + 1; //
      P(j,j) =  stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }

    // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){
      Matrix<double> tbeta = ::t(beta);
      for (int i=0; i<(ns*k); ++i)
	beta_store(count,i) = tbeta[i];
      for (int j=0; j<ns*ns; ++j)
	P_store(count,j)= P[j];
      for (int l=0; l<n ; ++l)
	ps_store(l,_) = ps_store(l,_) + ps(l,_);
      s_store(count,_) = s;
      Z_store(count,_) = Z;

      ++count;

    }   // end of if(iter...)


    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n\nMCMCprobitChange iteration %i of %i \n", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
	Rprintf("\n The number of observations in state %i is %10.5d", j+1, nstate[j]);
      }
      for (int j = 0; j<ns; ++j){
	Rprintf("\n beta  %i is ", j);
	for (int i = 0; i<k; ++i){
	  Rprintf("%10.5f", beta(j, i));
	}
      }
    }

    R_CheckUserInterrupt();

  }// end MCMC loop

  if(chib==1){
    Matrix<double> betast = meanc(beta_store);
    Matrix<double, Row> beta_st(ns, k);
    for (int j = 0; j<ns*k; ++j){
      beta_st[j] = betast[j];
    }
    Matrix<double> P_vec_st = meanc(P_store);
    const Matrix<double> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){
      P_st[j] = P_vec_st[j];
    }

    // 1. beta
    Matrix<double> density_beta(nstore, ns);
    for (int iter = 0; iter<nstore; ++iter){
      Matrix<int> nstate(ns, 1);
      int beta_count = 0;
      const Matrix<double> Z(n, 1);
      Z(_,0) = Z_store(iter,_);
       for (int j = 0; j<ns ; ++j){
	for (int i = 0; i<n; ++i){
	  if (s_store(iter, i) == (j+1)) {
	    nstate[j] = nstate[j] + 1;
	  }
	}
	beta_count = beta_count + nstate[j];
	const Matrix<double> Zj = Z((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	const Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	const Matrix<double> XpX = (::t(Xj)*Xj);
	const Matrix<double> XpZ = (::t(Xj)*Zj);
	const Matrix<double> Bn = invpd(B0 + XpX);
	const Matrix<double> bn = Bn*gaxpy(B0, b0, XpZ);
	density_beta(iter, j) = exp(lndmvn(::t(beta_st(j,_)), bn, Bn));
      }
    }

    double pdf_beta = log(prod(meanc(density_beta)));

    // 2. P
    Matrix<double> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      Matrix <double> Sout = probit_state_sampler(stream, m, Y, X, beta_st, P);
      Matrix <double> s = Sout(_, 0);
      Matrix <double> ps(n, ns);
      for (int j = 0; j<ns; ++j){
	ps(_,j) = Sout(_,j+1);
      }

      double shape1 = 0;
      double shape2 = 0;
      P(ns-1, ns-1) = 1;
      Matrix <double> P_addN(ns, 1);
      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<n; ++i){
	  if (s[i] == (j+1)) {
	    P_addN[j] = P_addN[j] + 1;
	  }
	}
      }
      for (int j =0; j< (ns-1); ++j){
	shape1 =  A0(j,j) + P_addN[j] - 1;
	shape2 =  A0(j,j+1) + 1; //
	P(j,j) = stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
	density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2);
      }
      density_P(iter, ns-1) = 1;

    }
    double pdf_P = log(prod(meanc(density_P)));

    // likelihood
    Matrix<double> F = Matrix<double>(n, ns);
    Matrix<double> like(n, 1);
    Matrix<double> pr1 = Matrix<double>(ns, 1);
    pr1[0] = 1;
    Matrix<double> py(ns, 1);
    Matrix<double> pstyt1(ns, 1);

    for (int t=0; t<n ; ++t){
      int yt = (int) Y[t];
      Matrix<double> mu = X(t,_)*::t(beta_st);
      for (int j=0; j<ns ; ++j){
	py[j]  =  dbinom(yt, 1, pnorm(mu[j], 0, 1));
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
    }

    loglike = sum(log(like));

    //////////////////////
    // log prior ordinate
    //////////////////////
    Matrix<double> density_beta_prior(ns, 1);
    Matrix<double> density_P_prior(ns, 1);
    density_P[ns-1] = 1; //

    for (int j=0; j<ns ; ++j){
      density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv);
    }

    for (int j =0; j< (ns-1); ++j){
      density_P_prior[j] = scythe::lndbeta1(P_st(j,j), A0(j,j), A0(j,j+1));
    }

    // compute marginal likelihood
    double logprior = sum(density_beta_prior) + sum(density_P_prior);
    logmarglike = (loglike + logprior) - (pdf_beta + pdf_P);

    if (verbose >0 ){
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("log_prior = %10.5f\n", logprior);
      Rprintf("log_beta = %10.5f\n", pdf_beta);
      Rprintf("log_P = %10.5f\n", pdf_P);
    }
  } // end of marginal likelihood
}//end extern "C"

extern "C"{
  void cMCMCprobitChange(double *betaout, double *Pout, double *psout, double *sout,
			const double *Ydata, const int *Yrow, const int *Ycol,
			const double *Xdata, const int *Xrow, const int *Xcol,
			const int *m,
			const int *burnin, const int *mcmc, const int *thin, const int *verbose,
			const int *uselecuyer, const int *seedarray, const int *lecuyerstream,
			const double *betastart,  const double *Pstart,
			const double *a, const double *b,
			const double *b0data, const double *B0data,
			const double *A0data,
			double *logmarglikeholder, double *loglikeholder,
			const int *chib){

    // pull together Matrix objects
    const Matrix <> Y(*Yrow, *Ycol, Ydata);
    const Matrix <> X(*Xrow, *Xcol, Xdata);
    const  int nstore = *mcmc / *thin;
    const int n = Y.rows();
    const int k = X.cols();
    const int ns = *m + 1;

    // generate starting values
    Matrix <> beta(ns, k, betastart);
    Matrix <> P(ns, ns, Pstart);
    Matrix <> b0(k, 1, b0data);
    Matrix <> B0(k, k, B0data);
    const Matrix <> A0(ns, ns, A0data);
    double logmarglike;
    double loglike;

    // storage matrices
    Matrix<> beta_store(nstore, ns*k);
    Matrix<> Z_store(nstore, n);
    Matrix<> P_store(nstore, ns*ns);
    Matrix<> ps_store(n, ns);
    Matrix<int> s_store(nstore, n);

    MCMCPACK_PASSRNG2MODEL(MCMCprobitChange_impl, *m,
			     Y, X, beta,  P, b0, B0, A0,
			     *burnin, *mcmc, *thin, *verbose, *chib,
			     beta_store, Z_store,
			     P_store, ps_store, s_store,
			     logmarglike, loglike);
    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;

    // return output
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
}
#endif


