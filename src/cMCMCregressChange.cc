////////////////////////////////////////////////////////////////////
// cMCMCregressChange.cc is a C++ code to estimate
// linear regression changepoint model
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr
//
// 03/03/2009 Written
////////////////////////////////////////////////////////////////////

#ifndef CMCMCREGRESSCHANGE_CC
#define CMCMCREGRESSCHANGE_CC

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

double lndinvgamma_jhp (const double x, const double shape, const double scale){
  double log_density = shape *::log(scale) - lngammafn(shape) - (shape + 1) * ::log(x) - (scale/x);
  return (log_density);
}
Matrix<double> loglike_fn(const int m,
			  const Matrix<double>& Y,
			  const Matrix<double>& X,
			  const Matrix<double>& beta,
			  const Matrix<double>& Sigma,
			  const Matrix<double>& P){

  const int ns = m + 1;
  const int n = Y.rows();
  Matrix<double> F = Matrix<double>(n, ns);
  Matrix<double> loglike(n, 1);
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
    loglike[t] = log(sum(unnorm_pstyt));
    //like[t] = sum(unnorm_pstyt);
  }// end of likelihood computation

  return loglike;
}


template <typename RNGTYPE>
Matrix<double> gaussian_state_sampler(rng<RNGTYPE>& stream,
				      const int m, const int sos,
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
  Matrix<int> s(n, 1);
  Matrix<double> ps = Matrix<double>(n, ns);

  for (int tt=0; tt<n ; ++tt){
    Matrix<double> mu = X(tt,_)*::t(beta);
    for (int j = 0; j< ns; ++j){
      py[j]  =  dnorm(Y[tt], mu[j], sqrt(Sigma[j]));
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
      else s(tt) = st;
    }
    ps(tt,_) = pstyn;
    --tt;
  }// end of while loop

  // Single Observation State (SOS)
  // we randomly sample all states
  // [TODO] It would be more efficient to sample neighboring states
  // JHP "Thu Mar  9 09:49:26 2017"
  if(sos == 1){
     int singleobsstate = 0.0;
    for (int i = 1; i<(n - 1); ++i){
      if ((s[i-1] != s[i]) && (s[i] != s[i+1])) {
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
  Matrix<double> Sout(n, ns+1);
  Sout(_, 0) = s(_,0);
  for (int j = 0; j<ns; ++j){
    Sout(_,j+1) = ps(_, j);
  }

  return Sout;

} // end of state sampler

////////////////////////////////////////////
// cMCMCregressChange implementation.
////////////////////////////////////////////
template <typename RNGTYPE>
void MCMCregressChange_impl(rng<RNGTYPE>& stream,
			    const double m,
			    const Matrix<>& Y, const Matrix<>& X,
			    Matrix<>& beta, Matrix<>& Sigma, Matrix<>& P, Matrix<int>& s,
			    Matrix<>& b0, Matrix<>& B0,
			    const double c0, const double d0,
			    const Matrix<>& A0,
			     int burnin,  int mcmc,  int thin,
			     int verbose, bool marginalrun, bool sos,
			    Matrix<>& beta_store, Matrix<>& Sigma_store,
			    Matrix<>& P_store, Matrix<>& ps_store, Matrix<int>& s_store,
			    Matrix<>& y_loglike_store,
			    double& logmarglike,
			    double& loglike)
{
  // define constants and form cross-product matrices
  const int tot_iter = burnin + mcmc; //total iterations
  const int nstore = mcmc / thin; // number of draws to store
  const int n = Y.rows();
  const int ns = m + 1;                 // number of states
  const int k = X.cols();
  const Matrix<> B0inv = invpd(B0);
  Matrix <> sigma(ns, 1);

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
      Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
      Matrix<double> Bn = invpd(B0 + t(Xj)*Xj/Sigma[j]);
      Matrix<double> bn = Bn*(B0*b0 + t(Xj)*yj/Sigma[j]);
      beta(j,_) = stream.rmvnorm(bn, Bn);

      // SIGMA UPDATE
      double shape = (c0 + (double)nstate[j])/2;
      Matrix<> ysimul_j = Xj*::t(beta(j,_)); //
      Matrix<double> ej = yj - ysimul_j;
      Matrix<double> sum_ej = t(ej)*ej;
      double scale =(d0 + sum_ej[0])/2;

      Sigma[j] = 1/stream.rgamma(shape, scale);
      sigma[j] = sqrt(Sigma[j]);
    }// end of sampling beta and Sigma

    //////////////////////
    // 2. Sample P
    //////////////////////
    double shape1 = 0;
    double shape2 = 0;
    P(ns-1, ns-1) = 1;

    for (int j =0; j<(ns-1); ++j){
      shape1 =  A0(j,j) + (double)nstate[j] - 1;
      shape2 =  A0(j,j+1) + 1; // SS(j,j+1);
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }

    //////////////////////
    // 3. Sample s
    //////////////////////
    Matrix<double> Sout = gaussian_state_sampler(stream, m, sos, Y, X, beta, Sigma, P);
    Matrix<double> s = Sout(_, 0);
    Matrix<> ps(n, ns);
    for (int j = 0; j<ns; ++j){
      ps(_,j) = Sout(_,j+1);
    }

    // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){
      Matrix<double> tbeta = ::t(beta); //transpose beta for R output
      for (int i=0; i<(ns*k); ++i)
	beta_store(count,i) = tbeta[i];// stored by the order of (11, 12, 13, 21, 22, 23)
      for (int i=0; i<ns; ++i)
	Sigma_store(count,i) = Sigma[i];
      for (int j=0; j<ns*ns; ++j)
	P_store(count,j)= P[j];
      for (int l=0; l<n ; ++l){
	ps_store(l,_) = ps_store(l,_) + ps(l,_);           // add two matrices
      }
      s_store(count,_) = s;
      y_loglike_store(count,_) = loglike_fn(m, Y, X, beta, Sigma, P);

      ++count; // count when (iter % *thin)==0

    }   // end of if(iter...)


    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\nMCMCregressChange iteration %i of %i \n", (iter+1), tot_iter);
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

  if(marginalrun ==1){

    Matrix<double> betast = meanc(beta_store); //meanc(beta_store)=(11, 12, 13, 21, 22, 23)
    Matrix<double, Row> beta_st(ns, k);
    for (int j = 0; j<ns*k; ++j){
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
	const double precision = 1.0/Sigma_store(iter, j);
	const Matrix<double> XpX = (::t(Xj)*Xj);
	const Matrix<double> XpY = (::t(Xj)*yj);
	const Matrix<double> Bn = invpd(B0 + XpX*precision);
	const Matrix<double> bn = Bn*gaxpy(B0, b0, XpY*precision);
	if (k == 1){
	  density_beta(iter, j) = dnorm(beta_st(j), bn(0), sqrt(Bn(0)));
	}
	else{
	  density_beta(iter, j) = ::exp(lndmvn(::t(beta_st(j,_)), bn, Bn));
	}
      }

    }

    double pdf_beta = log(prod(meanc(density_beta)));

    //////////////////////
    // 2. pdf.Sigma
    //////////////////////
    Matrix<double> density_Sigma(nstore, ns);
    for (int iter = 0; iter<nstore; ++iter){
      Matrix <double> Sout = gaussian_state_sampler(stream, m, sos, Y, X, beta_st, Sigma, P);
      Matrix <double> s = Sout(_, 0);
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
	Matrix<double> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	Matrix<double> ej = yj - Xj*::t(beta_st(j,_));
	Matrix<double> sum_ej = ::t(ej)*ej;
	double shape = (c0 + (double)nstate[j])/2;
	double scale =(d0 + sum_ej[0])/2;

	Sigma[j] = stream.rigamma(shape, scale);
	density_Sigma(iter, j) = ::exp(lndinvgamma_jhp(Sigma_st[j], shape, scale));

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
      Matrix <double> Sout = gaussian_state_sampler(stream, m, sos, Y, X, beta_st, Sigma_st, P);
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
      density_P(iter, ns-1) = 1; //without this, there will be an error

    }// end of pdf.P MCMC run
    double pdf_P = log(prod(meanc(density_P)));

    //////////////////////
    // likelihood
    //////////////////////
    loglike = sum(loglike_fn(m, Y, X, beta_st, Sigma_st, P_st));

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
      density_Sigma_prior[j] = lndinvgamma_jhp(Sigma_st[j], c0/2, d0/2);
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

////////////////////////////////////////////
// Start cMCMCregressChangepoint function
///////////////////////////////////////////
extern "C"{
  void cMCMCregressChange(double *betaout, double *Sigmaout, // double *Pout,
			 double *psout, double *sout,  double *yloglike,
			 const double *Ydata, const int *Yrow, const int *Ycol,
			 const double *Xdata, const int *Xrow, const int *Xcol,
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
    const Matrix <double> X(*Xrow, *Xcol, Xdata);
    // const  int tot_iter = *burnin + *mcmc; //total iterations
    const  int nstore = *mcmc / *thin; // number of draws to store
    const int n = Y.rows();
    const int k = X.cols();
    const int ns = *m + 1;                 // number of states

    // generate starting values
    Matrix <> beta(ns, k, betastart);
    Matrix <> Sigma(ns, 1, Sigmastart);
    Matrix <> P(ns, ns, Pstart);
    Matrix <int> s(n, 1, statestart);
    Matrix <> b0(k, 1, b0data);
    Matrix <> B0(k, k, B0data);
    const Matrix <> A0(ns, ns, A0data);
    double logmarglike;
    double loglike;

    // storage matrices
    Matrix<> beta_store(nstore, ns*k);
    Matrix<> Sigma_store(nstore, ns);
    Matrix<> P_store(nstore, ns*ns);
    Matrix<> ps_store(n, ns);
    Matrix<int> s_store(nstore, n);
    Matrix<> y_loglike_store(nstore, n);

    MCMCPACK_PASSRNG2MODEL(MCMCregressChange_impl, *m,
			   Y, X,
			   beta, Sigma, P, s, b0, B0,
			   *c0, *d0, A0,
			   *burnin, *mcmc, *thin, *verbose,
			   *marginalrun, *sos,
			   beta_store, Sigma_store,
			   P_store, ps_store, s_store,
			   y_loglike_store,
			   logmarglike, loglike);
    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;


    // return output
    for (int i = 0; i<(nstore*ns*k); ++i){
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

  }// end of MCMCpoissonChange
}//end extern "C"

#endif




