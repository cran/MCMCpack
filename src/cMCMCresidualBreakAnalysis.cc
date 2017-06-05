////////////////////////////////////////////////////////////////////
// cMCMCresidualBreakAnalysis.cc
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr
//
// Written 03/03/2009
//
////////////////////////////////////////////////////////////////////

#ifndef CMCMCRESIDUALBREAKANALYSIS_CC
#define CMCMCRESIDUALBREAKANALYSIS_CC

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
static double dinvgamma(double theta, double a, double b) {
  double logf =  a * log(b) - lngammafn(a) + -(a+1) * log(theta) +
                 -b/theta;
  return exp(logf);
}
Matrix<double> loglike_fn2(const int m,
			  const Matrix<double>& Y,
			  const Matrix<double>& X,
			  const Matrix<double>& beta,
			  const Matrix<double>& Sigma,
			  const Matrix<double>& P){

  const int ns = m + 1;
  const int n = Y.rows();
  Matrix<double> F = Matrix<double>(n, ns);
  Matrix<double> like(n, 1);
  Matrix<double> pr1 = Matrix<double>(ns, 1);
  pr1[0] = 1;
  Matrix<double> py(ns, 1);
  Matrix<double> pstyt1(ns, 1);

  for (int t=0; t<n ; ++t){
    //int yt = (int) Y[t];
    Matrix<double> mu = X(t,_)*::t(beta); //k by 1 vector
    for (int j = 0; j< ns; ++j){
      py[j]  =  dnorm(Y[t], mu[j], sqrt(Sigma[j]));
    }
    if (t==0) pstyt1 = pr1;
    else {
      pstyt1 =  ::t(F(t-1,_)*P); // make it an ns by 1 matrix
    }
    Matrix<double> unnorm_pstyt = pstyt1%py;
    Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
    for (int j=0; j<ns ; ++j){
      F(t,j) = pstyt(j);
    }
    like[t] = sum(unnorm_pstyt);
  }// end of likelihood computation

  return like;
}


////////////////////////////////////////////
// Start cMCMCresidualBreakAnalysispoint function
///////////////////////////////////////////
template <typename RNGTYPE>
void MCMCresidualBreakAnalysis_impl(rng<RNGTYPE>& stream,
				    const double m,
				    const Matrix<>& Y, Matrix<>& beta,
				    Matrix<>& Sigma, Matrix<>& P, Matrix<int>& s,
				    const double b0, const double B0,
				    const double c0, const double d0,
				    const Matrix<>& A0,
				     int burnin,  int mcmc,  int thin,
				     int verbose, bool marginalrun, bool sos,
				    Matrix<>& beta_store, Matrix<>& Sigma_store,
				    Matrix<>& P_store, Matrix<>& ps_store, Matrix<int>& s_store,
				    Matrix<>& y_loglike_store,
				    double& logmarglike, double&loglike)
{
  // define constants and form cross-product matrices
  const int tot_iter = burnin + mcmc;
  const int nstore = mcmc / thin;
  const int n = Y.rows();
  const int ns = m + 1;
  const double B0inv = 1/B0;
  const Matrix<double> X(n, 1, true, 1);
  Matrix<> sigma(ns, 1);

  //MCMC loop
   int count = 0;
  for ( int iter = 0; iter < tot_iter; ++iter){

    //////////////////////
    // 1. Sample beta and Sigma
    //////////////////////
    int beta_count = 0;
    Matrix<int> nstate(ns, 1);

    for (int j = 0; j<ns; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == j + 1) {
	  nstate[j] = nstate[j] + 1;
	}// end of if
      }// end of int i<n
      beta_count = beta_count + nstate[j];

      // BETA UPDATE
      Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      double Bn = 1/(B0 + (double)nstate[j]/Sigma[j]);
      double bn = Bn*(B0*b0 + sum(yj)/Sigma[j]);
      beta(j) = stream.rnorm(bn, sqrt(Bn));

      // SIGMA UPDATE
      double shape = (c0 + (double)nstate[j])/2;
      const Matrix<> ej(nstate[j], 1);
      for (int i = 0; i<nstate[j]; ++i){
	ej(i) = yj(i) - beta(j);
      }
      const Matrix<> SSE = crossprod (ej);
      double scale =(d0 + SSE[0])/2;

      Sigma[j] = 1/stream.rgamma(shape, scale);
      sigma[j] = sqrt(Sigma[j]);
    }

    //////////////////////
    // 2. Sample P
    //////////////////////
    double shape1 = 0;
    double shape2 = 0;
    P(ns-1, ns-1) = 1;

    for (int j =0; j<(ns-1); ++j){
      shape1 =  A0(j,j) + (double)nstate[j] - 1;
      shape2 =  A0(j,j+1) + 1;
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
    Matrix<double> ps = Matrix<double>(n, ns);
    for (int tt=0; tt<n ; ++tt){
      for (int j = 0; j< ns; ++j){
	py[j]  =  dnorm(Y[tt], beta[j], sigma[j]);
      }
      if (tt==0) pstyt1 = pr1;
      else {
	pstyt1 =  ::t(F(tt-1,_)*P);
      }
      Matrix<double> unnorm_pstyt = pstyt1%py;
      const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
      for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);

    }
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
	else s(tt) = st;// stay
      }

      // state singleobsstate case
      // we randomly sample all states
      // [TODO] It would be more efficient to sample neighboring states
      // JHP "Thu Mar  9 09:49:26 2017"
      if (sos == 1){
	int singleobsstate = 0.0;
	for (int i = 1; i<n; ++i){
	  if (s[i] != s[i-1]) {
	    singleobsstate = singleobsstate + 1;
	  }// end of if
	}// end of int i<n

	if(singleobsstate > 0){
	  // Rprintf("\n Single Observation State (SOS) occurs and hence states are resampled randomly. \n");
	  for (int i = 1; i<(n - 1); ++i){
	    s[i] = stream.rbinom(m, 0.5) + 1;
	  }
	  s = sort(s);
	  // s <- sort(sample(1:ns, T, replace=TRUE, prob=rep(1/ns, ns)))
	}
      }

      ps(tt,_) = pstyn;
      --tt;
    }// end of while loop

    // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){
      for (int i=0; i<ns; ++i)
	beta_store(count,i) = beta[i];
      for (int i=0; i<ns; ++i)
	Sigma_store(count,i) = Sigma[i];
      for (int j=0; j<ns*ns; ++j)
	P_store(count,j)= P[j];
      for (int l=0; l<n ; ++l)
	ps_store(l,_) = ps_store(l,_) + ps(l,_);
      s_store(count,_) = s;
      y_loglike_store(count,_) = loglike_fn2(m, Y, X, beta, Sigma, P);

      ++count;

    }


    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n testpanelSubjectBreak iteration %i of %i \n", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
	Rprintf("\n The number of observations in state %i is %10.5f", j+1, static_cast<double>(nstate[j]));
      }
      Rprintf("\n beta \n");
      for (int i = 0; i<ns; ++i){
	Rprintf("%10.5f\n", beta(i));
      }
      Rprintf("\n sigma^2 \n");
      for (int i = 0; i<ns; ++i){
	Rprintf("%10.5f\n", Sigma(i));
      }
    }

  }// end MCMC loop

  if(marginalrun ==1){

    Matrix<double> betast = meanc(beta_store);
    Matrix<double, Row> beta_st(ns, 1);
    for (int j = 0; j<ns; ++j){
      beta_st[j] = betast[j];
    }

    Matrix<double> Sigma_st = meanc(Sigma_store);
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
	const double precision = 1.0/Sigma_store(iter, j);
	const double Bn = 1/(B0 + (double)nstate[j]*precision);
	double bn = Bn*(B0*b0 + sum(yj)*precision);
	density_beta(iter, j) = dnorm(beta_st(j), bn, sqrt(Bn));
      }
    }// end of pdf.beta

    double pdf_beta = log(prod(meanc(density_beta)));

    //////////////////////
    // 2. pdf.Sigma
    //////////////////////
    Matrix<double> density_Sigma(nstore, ns);
    for (int iter = 0; iter<nstore; ++iter){
      Matrix<double> F(n, ns);
      Matrix<double> pr1(ns, 1);
      pr1[0] = 1;
      Matrix<double> py(ns, 1);
      Matrix<double> pstyt1(ns, 1);
      Matrix<double> ps = Matrix<double>(n, ns);
      for (int tt=0; tt<n ; ++tt){
	for (int j = 0; j< ns; ++j){
	  py[j]  =  dnorm(Y[tt], beta_st[j], sigma[j]);
	}
	if (tt==0) pstyt1 = pr1;
	else {
	  pstyt1 =  ::t(F(tt-1,_)*P);
	}
	Matrix<double> unnorm_pstyt = pstyt1%py;
	const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
	for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);

      }
      ps(n-1,_) = F(n-1,_);


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

      int beta_count = 0;
      Matrix<int> nstate(ns, 1);

      for (int j = 0; j <ns ; ++j){
	for (int i = 0; i<n; ++i){
	  if (s[i] == (j+1)) {
	    nstate[j] = nstate[j] + 1;
	  }// end of if
	}// end of int i<n
	beta_count = beta_count + nstate[j];

	Matrix<double> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	const Matrix<> ej(nstate[j], 1);
	for (int i = 0; i<nstate[j]; ++i){
	  ej(i) = yj(i) - beta_st(j);
	}
	const Matrix<> SSE = crossprod (ej);
	double scale =(d0 + SSE[0])/2;
	double shape = (c0 + (double)nstate[j])/2;

	Sigma[j] = stream.rigamma(shape, scale);
	sigma[j] = sqrt(Sigma[j]);
	density_Sigma(iter, j) = dinvgamma(Sigma_st[j], shape, scale);
      }

      double shape1 = 0;
      double shape2 = 0;
      P(ns-1, ns-1) = 1;

      for (int j =0; j< (ns-1); ++j){
	shape1 =  A0(j,j) + (double)nstate[j] - 1;
	shape2 =  A0(j,j+1) + 1; //
	P(j,j) =  stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
      }

    }// end of pdf.Sigma

    double pdf_Sigma = log(prod(meanc(density_Sigma)));

    // 3. pdf.P|beta_st, Sigma_st, S
    Matrix<double> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){

      Matrix<double> F(n, ns);
      Matrix<double> pr1(ns, 1);
      pr1[0] = 1;
      Matrix<double> py(ns, 1);
      Matrix<double> pstyt1(ns, 1);
      Matrix<double> ps = Matrix<double>(n, ns);
      for (int tt=0; tt<n ; ++tt){
	for (int j = 0; j< ns; ++j){
	  py[j]  =  dnorm(Y[tt], beta_st[j], sqrt(Sigma_st[j]));
	}
	if (tt==0) pstyt1 = pr1;
	else {
	  pstyt1 =  ::t(F(tt-1,_)*P);
	}
	Matrix<double> unnorm_pstyt = pstyt1%py;
	const Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt);
	for (int j=0; j<ns ; ++j) F(tt,j) = pstyt(j);

      }
      ps(n-1,_) = F(n-1,_);

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
	  else s(tt) = st;// stay
	}
	ps(tt,_) = pstyn;
	--tt;
      }// end of while loop

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

    }
    double pdf_P = log(prod(meanc(density_P)));

    //////////////////////
    // likelihood
    //////////////////////
    loglike = sum(loglike_fn2(m, Y, X, beta_st, Sigma_st, P_st));

    //////////////////////
    // log prior ordinate
    //////////////////////
    Matrix<double> density_beta_prior(ns, 1);
    Matrix<double> density_Sigma_prior(ns, 1);
    Matrix<double> density_P_prior(ns, 1);
    density_P[ns-1] = 1; //
    for (int j=0; j<ns ; ++j){
	density_beta_prior[j] = log(dnorm(beta_st(j), b0, sqrt(B0inv)));
    }
    for (int j=0; j<ns ; ++j){
      density_Sigma_prior[j] = log(dinvgamma(Sigma_st[j], c0/2, d0/2));
    }
    for (int j =0; j< (ns-1); ++j){
      density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1)));
    }

    // compute marginal likelihood
    double logprior = sum(density_beta_prior) + sum(density_Sigma_prior) + sum(density_P_prior);

    logmarglike = (loglike + logprior) - (pdf_beta + pdf_Sigma + pdf_P);
    if(verbose > 0){
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("logprior = %10.5f\n", logprior);
      Rprintf("pdf_beta = %10.5f\n", pdf_beta);
      Rprintf("pdf_Sigma = %10.5f\n", pdf_Sigma);
      Rprintf("pdf_P = %10.5f\n", pdf_P);
    }
  }// end of marginal likelihood computation

}

extern "C"{
  void cMCMCresidualBreakAnalysis(double *betaout, double *Sigmaout, // double *Pout,
				 double *psout,  double *sout, double *yloglike,
				 const double *Ydata, const int *Yrow, const int *Ycol,
				 const int *m,
				 const int *burnin, const int *mcmc, const int *thin, const int *verbose,
				 const int *uselecuyer, const int *seedarray, const int *lecuyerstream,
				 const double *betastart, const double *Sigmastart, const double *Pstart,
				 const int *statestart,
				 const double *a, const double *b,
				 const double *b0data, const double *B0data,
				 const double *c0, const double *d0,
				 const double *A0data,
				 double *logmarglikeholder, double *loglikeholder,
				 const int *marginalrun, const int *sos){

    // pull together Matrix objects
    const Matrix <double> Y(*Yrow, *Ycol, Ydata);
    const  int nstore = *mcmc / *thin;
    const int n = Y.rows();
    const int ns = *m + 1;

    // generate starting values
    Matrix <> beta(ns, 1, betastart);
    Matrix <> Sigma(ns, 1, Sigmastart);
    Matrix <> P(ns, ns, Pstart);
    Matrix <int> s(n, 1, statestart);
    const Matrix <> A0(ns, ns, A0data);
    double logmarglike;
    double loglike;

    // storage matrices
    Matrix<> beta_store(nstore, ns);
    Matrix<> Sigma_store(nstore, ns);
    Matrix<> P_store(nstore, ns*ns);
    Matrix<> ps_store(n, ns);
    Matrix<int> s_store(nstore, n);
    Matrix<> y_loglike_store(nstore, n);

    MCMCPACK_PASSRNG2MODEL(MCMCresidualBreakAnalysis_impl, *m,
			   Y, beta, Sigma, P, s,
			   *b0data, *B0data, *c0, *d0, A0,
			   *burnin, *mcmc, *thin, *verbose,
			   *marginalrun, *sos,
			   beta_store, Sigma_store,
			   P_store, ps_store, s_store, y_loglike_store,
			   logmarglike, loglike);
    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;

    // return output
    for (int i = 0; i<(nstore*ns); ++i){
      betaout[i] = beta_store[i];
    }
    for (int i = 0; i<(nstore*ns); ++i){
      Sigmaout[i] = Sigma_store[i];
    }
    // for (int i = 0; i<(nstore*ns*ns); ++i){
    //   Pout[i] = P_store[i];
    // }
    for (int i = 0; i<(n*ns); ++i){
      psout[i] = ps_store[i];
    }
    for (int i = 0; i<(nstore*n); ++i){
      sout[i] = s_store[i];
      yloglike[i] = y_loglike_store[i];
    }

  }// end of cMCMCresidualBreakAnalysis
}//end extern "C"

#endif




