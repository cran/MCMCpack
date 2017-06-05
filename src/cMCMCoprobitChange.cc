////////////////////////////////////////////////////////////////////
// cMCMCoprobitChange.cc is C++ code to estimate a oprobit changepoint model
// with linear approximation
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr
//
// 07/06/2007 Written
// 11/02/2009 Modified
////////////////////////////////////////////////////////////////////

#ifndef CMCMCOPROBITCHANGE_CC
#define CMCMCOPROBITCHANGE_CC

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

// density function for truncated normal
static double dtnormLX(const double x,
		     const double mean,
		     const double sd,
		     const double lower,
		     const double upper){
  double out = 0.0;
  if (x>lower && x<upper){
    double numer = dnorm(x, mean, sd);
    double denom = pnorm(upper, mean, sd) - pnorm(lower, mean, sd);
    out = numer/denom;
  }
  else{
    Rprintf("\n x input for dtnormLX() %10.5f is out of bounds %10.5f %10.5f ", x, lower, upper, "\n");
    out = 1;
  }
  return out;
}

// likelihood of oprobit
static double oprobit_pdfLX(const int ncat,
			    const int Y,
			    double Xbeta,
			    const Matrix<>& gamma){

  Matrix<> cat_prob(1, ncat-1);
  Matrix<> prob(1, ncat);
  for (int j=0; j< ncat-1; ++j){
    cat_prob(0, j) = pnorm(gamma[j+1] - Xbeta, 0.0, 1.0);
  }
  prob(0, ncat-1) = 1 - cat_prob(0, ncat-2);
  prob(0, 0) = cat_prob(0, 0);
  for (int j=1; j<(ncat-1); ++j){
    prob(0, j) = cat_prob(0,j) - cat_prob(0, j-1);
  }
  double like = prob(0,Y-1);

  return like;
}

static double oprobit_log_postLX( int j,
				 const int ncat,
				 const Matrix<>& gamma_p,
				 const Matrix<>& gamma,
				 const Matrix<>& Y,
				 const Matrix<>& X,
				 const Matrix<>& beta,
				 const Matrix<>& tune,
				 const int gammafixed){
  const int N = Y.rows();
  double loglikerat = 0.0;
  double loggendenrat = 0.0;
  Matrix<> Xbeta = X*t(beta(j,_));

  if (gammafixed==1){
    for ( int i=0; i<N; ++i){
      int yi = Y[i];
      if (yi == ncat){
	loglikerat = loglikerat
	  + log(1.0  - pnorm(gamma_p(yi-1) - Xbeta[i], 0.0, 1.0) )
	  - log(1.0 - pnorm(gamma(yi-1) - Xbeta[i], 0.0, 1.0) );
      }
      else if (yi == 1){
	loglikerat = loglikerat + log(pnorm(gamma_p(yi) - Xbeta[i], 0.0, 1.0)  )
	  - log(pnorm(gamma(yi)  - Xbeta[i], 0.0, 1.0) );
      }
      else{
	loglikerat = loglikerat
	  + log(pnorm(gamma_p(yi) - Xbeta[i], 0.0, 1.0) -
		pnorm(gamma_p(yi-1) - Xbeta[i], 0.0, 1.0) )
	  - log(pnorm(gamma(yi) - Xbeta[i], 0.0, 1.0) -
		pnorm(gamma(yi - 1)  - Xbeta[i], 0.0, 1.0) );
      }
    }
    for ( int k =2; k <ncat; ++k){
      loggendenrat = loggendenrat +
	log(pnorm(gamma(k+1), gamma(k), tune[j]) -
	    pnorm(gamma_p(k-1), gamma(k), tune[j]))  -
	log(pnorm(gamma_p(k+1), gamma_p(k), tune[j]) +
	    pnorm(gamma(k-1), gamma_p(k), tune[j]));
    }
  }
  else{
    for ( int i=0; i<N; ++i){
      int yi = Y[i];
      if (yi == ncat){
	loglikerat = loglikerat
	  + log(1.0  - pnorm(gamma_p(j, yi-1) - Xbeta[i], 0.0, 1.0) )
	  - log(1.0 - pnorm(gamma(j, yi-1) - Xbeta[i], 0.0, 1.0) );
      }
      else if (yi == 1){
	loglikerat = loglikerat + log(pnorm(gamma_p(j, yi) - Xbeta[i], 0.0, 1.0)  )
	  - log(pnorm(gamma(j, yi)  - Xbeta[i], 0.0, 1.0) );
      }
      else{
	loglikerat = loglikerat
	  + log(pnorm(gamma_p(j, yi) - Xbeta[i], 0.0, 1.0) -
		pnorm(gamma_p(j, yi-1) - Xbeta[i], 0.0, 1.0) )
	  - log(pnorm(gamma(j, yi) - Xbeta[i], 0.0, 1.0) -
		pnorm(gamma(j, yi - 1)  - Xbeta[i], 0.0, 1.0) );
      }
    }
    for ( int k =2; k <ncat; ++k){
      loggendenrat = loggendenrat +
	log(pnorm(gamma(j ,k+1), gamma(j, k), tune[j]) -
	    pnorm(gamma_p(j, k-1), gamma(j, k), tune[j]))  -
	log(pnorm(gamma_p(j, k+1), gamma_p(j, k), tune[j]) +
	    pnorm(gamma(j, k-1), gamma_p(j, k), tune[j]));
    }
  }
  return loglikerat + loggendenrat;
}

template <typename RNGTYPE>
Matrix<> gaussian_ordinal_state_sampler_fixedsigma(rng<RNGTYPE>& stream,
						   const int m,
						   const Matrix<>& Y,
						   const Matrix<>& X,
						   const Matrix<>& beta,
						   const Matrix<>& Sigma,
						   const Matrix<>& P){

  const int ns = m + 1;
  const int n = Y.rows();
  Matrix<> F(n, ns);
  Matrix<> pr1(ns, 1);
  pr1[0] = 1;
  Matrix<> py(ns, 1);
  Matrix<> pstyt1(ns, 1);

  for (int t=0; t<n ; ++t){
    Matrix<> mu = X(t,_)*::t(beta);
    for (int j = 0; j< ns; ++j){
      py[j]  =  dnorm(Y[t], mu[j], sqrt(Sigma[0]));
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
// cMCMCoprobitChangeLinearAppxpoint implementation.
////////////////////////////////////////////
template <typename RNGTYPE>
void MCMCoprobitChange_impl(rng<RNGTYPE>& stream,
			    const int m, const int ncat,
			    const Matrix<>& Y, const Matrix<>& X,
			    Matrix<>& beta, Matrix<>& beta_linear, Matrix<>& gamma, Matrix<>& P,
			    Matrix<>& Sigma,
			    Matrix<>& b0, Matrix<>& B0, const Matrix<>& A0,
			     int burnin,  int mcmc,  int thin,
			     int verbose, const Matrix<>& tune, // const Matrix<int>& tdf,
			    bool chib,  bool gammafixed,
			    Matrix<>& beta_store,  Matrix<>& beta_linear_store,
			    Matrix<>& gamma_store, Matrix<>& Z_store,
			    Matrix<>& P_store, Matrix<>& ps_store, Matrix<int>& s_store,
			    double& logmarglike, double& loglike)
{
  const  int tot_iter = burnin + mcmc;
  const  int nstore = mcmc / thin;
  const  int N = Y.rows();
  const  int ns = m + 1;
  const  int k = X.cols();
  const  int gk = ncat + 1;
  const Matrix<> B0inv = invpd(B0);
  Matrix<> Z(N, 1);
  Matrix<> accepts(ns, 1);
  Matrix<> gamma_p = gamma;

  //MCMC loop
   int count = 0;
  for (int iter = 0; iter < tot_iter; ++iter){
    // 1. Sample s
    Matrix<> Sout = gaussian_ordinal_state_sampler_fixedsigma(stream, m, Y, X, beta_linear, Sigma, P);
    Matrix<int> s = Sout(_, 0);
    Matrix <double> ps(N, ns);
    for (int j = 0; j<ns; ++j){
      ps(_,j) = Sout(_,j+1);
    }
    // 2. Sample Z
    for ( int i = 0; i<N; ++i) {
      Matrix<> mu = X(i,_)*t(beta(s[i]-1,_));
       int yi = Y[i];
      Z[i] = stream.rtnorm_combo(mu[0], 1.0, gamma(s[i]-1, yi-1), gamma(s[i]-1, yi));
    }

    // 3. Sample beta
    int beta_count = 0;
    Matrix<int> beta_count_storage(ns, 1);
    Matrix<int> nstate(ns, 1);
    for (int j = 0; j<ns; ++j){
      for (int i = 0; i<N; ++i){
	if (s[i] == j + 1) {
	  nstate[j] = nstate[j] + 1;
	}
      }
      beta_count = beta_count + nstate[j];
      Matrix<> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      Matrix<> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
      Matrix<> Zj = Z((beta_count - nstate[j]), 0, (beta_count - 1), 0);
      Matrix<> XpX = t(Xj)*Xj;
      Matrix<> XpZ = t(Xj)*Zj;
      Matrix<> XpY = t(Xj)*yj;
      Matrix<> Bn = invpd(B0 + XpX/Sigma[0]);
      Matrix<> bn = Bn*(B0*b0 + XpY/Sigma[0]);
      beta_linear(j,_) = stream.rmvnorm(bn, Bn);
      Matrix<> Bn2 = invpd(B0 + XpX);
      Matrix<> bn2 = Bn2*(B0*b0 + XpZ);
      beta(j,_) = stream.rmvnorm(bn2, Bn2);
      beta_count_storage[j] = beta_count;
    }

    // 4. Sample gamma
    for (int j = 0; j < ns ; ++j){
      for (int i = 2; i< ncat; ++i){
	if (i==(ncat-1)){
	  gamma_p(j, i) = stream.rtbnorm_combo(gamma(j, i), ::pow(tune[j], 2.0),
					       gamma_p(j, i-1));
	}
	else {
	  gamma_p(j, i) = stream.rtnorm_combo(gamma(j, i), ::pow(tune[j], 2.0),
					      gamma_p(j, i-1), gamma(j, i+1));
	}
      }
      Matrix<> Yj = Y((beta_count_storage[j] - nstate[j]), 0, (beta_count_storage[j] - 1), 0);
      Matrix<> Xj = X((beta_count_storage[j] - nstate[j]), 0, (beta_count_storage[j] - 1), k-1);
      double alpha = oprobit_log_postLX(j, ncat, gamma_p, gamma, Yj, Xj, beta, tune, gammafixed);

      if (stream.runif() <= exp(alpha)){
	gamma(j,_) = gamma_p(j,_);
	accepts[j] = accepts[j] + 1;
      }
    }

    // 5. Sample P
    double shape1 = 0;
    double shape2 = 0;
    P(ns-1, ns-1) = 1;

    for (int j =0; j< (ns-1); ++j){
      shape1 =  A0(j,j) +  (double)nstate[j] - 1;
      shape2 =  A0(j,j+1) + 1; //
      P(j,j) =  stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }

     // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){
      Matrix<> tbeta = ::t(beta);
      Matrix<> tbetaLX = ::t(beta_linear);
      for (int i=0; i<(ns*k); ++i){
	beta_store(count,i) = tbeta[i];
	beta_linear_store(count,i) = tbetaLX[i];
      }
      Matrix<> tgamma = ::t(gamma);
      for (int i=0; i<(ns*gk); ++i)
	gamma_store(count,i) = tgamma[i];
      for (int j=0; j<ns*ns; ++j)
	P_store(count,j)= P[j];
      for (int l=0; l<N ; ++l)
	ps_store(l,_) = ps_store(l,_) + ps(l,_);
      s_store(count,_) = s;
      Z_store(count,_) = Z;

      ++count;

    }   // end of if(iter...)


    // print output to stdout
    if(iter > 1 && verbose > 0 && iter % verbose == 0){
      Rprintf("\n\nMCMCoprobitChange iteration %i of %i \n", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
      	double acceptancerate = accepts[j]/iter;
	Rprintf("\n\n Acceptance rate for state %i is %10.5f \n", j+1, acceptancerate);
      }
      for (int j = 0;j<ns; ++j){
	Rprintf("\n The number of observations in state %i is %10.5d", j+1, nstate[j]);
      }
      for (int j = 0; j<ns; ++j){
	Rprintf("\n beta %i = ", j);
	for (int i = 0; i<k; ++i){
	  Rprintf("%10.5f", beta(j, i));
	}
      }
      // for (int j = 0; j<ns; ++j){
      //	Rprintf("\n beta_linear %i = ", j);
      //	for (int i = 0; i<k; ++i){
      //	  Rprintf("%10.5f", beta_linear(j, i));
      //	}
      // }
      for (int j = 0; j<ns; ++j){
	Rprintf("\n gamma %i = ", j);
	for (int i = 0; i<(ncat-2); ++i){
	  Rprintf("%10.5f", gamma(j, i+2));
	}
      }
    }

    R_CheckUserInterrupt();

  }// end MCMC loop

  if(chib==1){
    Matrix<> betast = meanc(beta_store);
    Matrix<> betastLX = meanc(beta_linear_store);
    Matrix<double, Row> beta_st(ns, k);
    Matrix<double, Row> beta_linear_st(ns, k);
    for (int j = 0; j<ns*k; ++j){
      beta_st[j] = betast[j];
      beta_linear_st[j] = betastLX[j];
    }
    Matrix<> gammast = meanc(gamma_store);
    Matrix<double, Row> gamma_st(ns, gk);
    for (int j = 0; j<ns*gk; ++j){
      gamma_st[j] = gammast[j];
    }

    // for (int j = 0; j<ns; ++j){
    //  Rprintf("\n gamma_st %i = ", j);
    //  for (int i = 0; i<gk; ++i){
    //	Rprintf("%10.5f", gamma_st(j, i));
    //  }
    // }

    Matrix<> P_vec_st = meanc(P_store);
    const Matrix<> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){
      P_st[j] = P_vec_st[j];
    }
    // storage
    Matrix<> pdf_numer_store(nstore, 1);
    Matrix<> pdf_alpha_store(nstore, 1);
    Matrix<> pdf_P_store(nstore, ns);

     // 1. gamma
    Matrix<> densityq(nstore, ns);
    Matrix<> alpha(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      int beta_count = 0;
       Matrix<int> nstate(ns, 1);

      Matrix<double, Row> gamma_g(ns, gk);
      for (int h = 0; h<(ns*gk); ++h){
	gamma_g[h] = gamma_store(iter, h);
      }

      Matrix<double, Row> beta_g(ns, k);
      for (int h = 0; h<(ns*k); ++h){
	beta_g[h] = beta_store(iter, h);
      }
      Matrix<> pdf_numer(ns, 1);
      for (int j = 0; j <ns ; ++j){
	for (int i = 0; i<N; ++i){
	  if (s_store(iter, i) == (j+1)) {
	    nstate[j] = nstate[j] + 1;
	  }
	}
	beta_count = beta_count + nstate[j];

 	Matrix<int> Yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	pdf_numer(j) = oprobit_log_postLX(j, ncat, gamma_st, gamma_g, Yj, Xj, beta_g,
						      tune, gammafixed);
	for (int h = 2; h<ncat; ++h){
	  if (h == (ncat-1)){
	    densityq(iter, j) = densityq(iter, j) + log(dtnormLX(gamma_st(j, h), gamma_g(j, h), tune[j],
								 gamma_st(j, h-1), 300));
	  }
	  else {
	    densityq(iter, j) = densityq(iter, j) + log(dtnormLX(gamma_st(j, h), gamma_g(j, h), tune[j],
								 gamma_st(j, h-1), gamma_g(j, h+1)));
	  }
	}
      }

     if (sum(pdf_numer) > 0){
	pdf_numer_store(iter) = 0;
      }
     else{
       pdf_numer_store(iter) = sum(pdf_numer);
     }

    }
    double numerator = sum(meanc(pdf_numer_store))  + sum(meanc(densityq));

    for (int iter = 0; iter < nstore; ++iter){
      Matrix<> Sout = gaussian_ordinal_state_sampler_fixedsigma(stream, m, Y, X, beta_linear, Sigma, P);
      Matrix<int> s = Sout(_, 0);
      for ( int i = 0; i<N; ++i) {
	const Matrix<> mu = X(i,_)*t(beta(s[i]-1,_));
	Z[i] = stream.rtnorm_combo(mu[0], 1.0, gamma_st(s[i]-1, Y[i]-1), gamma_st(s[i]-1, Y[i]));
      }
      int beta_count = 0;
      Matrix<int> beta_count_storage(ns, 1);
      Matrix<int> nstate(ns, 1);
      for (int j = 0; j <ns ; ++j){
	for (int i = 0; i<N; ++i){
	  if (s[i] == (j+1)) {
	    nstate[j] = nstate[j] + 1;
	  }
	}
	beta_count = beta_count + nstate[j];
	Matrix<> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	Matrix<> Zj = Z((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> XpX = t(Xj)*Xj;
	Matrix<> XpZ = t(Xj)*Zj;
	Matrix<> XpY = t(Xj)*yj;
	Matrix<> Bn = invpd(B0 + XpX/Sigma[0]);
	Matrix<> bn = Bn*(B0*b0 + XpY/Sigma[0]);
	beta_linear(j,_) = stream.rmvnorm(bn, Bn);
	Matrix<> Bn2 = invpd(B0 + XpX);
	Matrix<> bn2 = Bn2*(B0*b0 + XpZ);
	beta(j,_) = stream.rmvnorm(bn2, Bn2);
	beta_count_storage[j] = beta_count;
     }
      // Sample P
      double shape1 = 0;
      double shape2 = 0;
      P(ns-1, ns-1) = 1;

      for (int j =0; j< (ns-1); ++j){
	shape1 =  A0(j,j) + nstate[j] - 1;
	shape2 =  A0(j,j+1) + 1; //
	P(j,j) =  stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
      }
      Matrix<> alpha(ns, 1);
      for (int j = 0; j < ns ; ++j){
	for (int i = 2; i< ncat; ++i){
	  if (i==(ncat-1)){
	    gamma_p(j, i) = stream.rtbnorm_combo(gamma_st(j, i), ::pow(tune[j], 2.0),
						 gamma_p(j, i-1));
	  }
	  else {
	    gamma_p(j, i) = stream.rtnorm_combo(gamma_st(j, i), ::pow(tune[j], 2.0),
						gamma_p(j, i-1), gamma_st(j, i+1));
	  }
	}
	Matrix<int> Yj = Y((beta_count_storage[j] - nstate[j]), 0, (beta_count_storage[j] - 1), 0);
	Matrix<> Xj = X((beta_count_storage[j] - nstate[j]), 0, (beta_count_storage[j] - 1), k-1);
	alpha[j] = oprobit_log_postLX(j, ncat, gamma_p, gamma_st, Yj, Xj, beta, tune, gammafixed);
      }
      if (sum(alpha) > 0){
	pdf_alpha_store(iter) = 0;
      }
      else{
	pdf_alpha_store(iter) = sum(alpha);
      }

    }
    double denominator = mean(pdf_alpha_store);
    double pdf_gamma = numerator - denominator;

    // 2. beta
    Matrix<> density_beta(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      Matrix<> Sout = gaussian_ordinal_state_sampler_fixedsigma(stream, m, Y, X, beta_linear, Sigma, P);
      Matrix<int> s = Sout(_, 0);
      for ( int i = 0; i<N; ++i) {
	const Matrix<> mu = X(i,_)*t(beta(s[i]-1,_));
	Z[i] = stream.rtnorm_combo(mu[0], 1.0, gamma_st(s[i]-1, Y[i]-1), gamma_st(s[i]-1, Y[i]));
      }
      int beta_count = 0;
      Matrix<int> beta_count_storage(ns, 1);
      Matrix<int> nstate(ns, 1);
      for (int j = 0; j <ns ; ++j){
	for (int i = 0; i<N; ++i){
	  if (s[i] == (j+1)) {
	    nstate[j] = nstate[j] + 1;
	  }
	}
	beta_count = beta_count + nstate[j];
	Matrix<> yj = Y((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> Xj = X((beta_count - nstate[j]), 0, (beta_count - 1), k-1);
	Matrix<> Zj = Z((beta_count - nstate[j]), 0, (beta_count - 1), 0);
	Matrix<> XpX = t(Xj)*Xj;
	Matrix<> XpZ = t(Xj)*Zj;
	Matrix<> XpY = t(Xj)*yj;
	Matrix<> Bn = invpd(B0 + XpX/Sigma[0]);
	Matrix<> bn = Bn*(B0*b0 + XpY/Sigma[0]);
	beta_linear(j,_) = stream.rmvnorm(bn, Bn);
	Matrix<> Bn2 = invpd(B0 + XpX);
	Matrix<> bn2 = Bn2*(B0*b0 + XpZ);
	beta(j,_) = stream.rmvnorm(bn2, Bn2);
	beta_count_storage[j] = beta_count;
	density_beta(iter, j) = exp(lndmvn(t(beta_st(j,_)), bn2, Bn2));
      }

      // Sample P
      double shape1 = 0;
      double shape2 = 0;
      P(ns-1, ns-1) = 1;

      for (int j =0; j< (ns-1); ++j){
	shape1 =  A0(j,j) + nstate[j] - 1;
	shape2 =  A0(j,j+1) + 1; //
	P(j,j) =  stream.rbeta(shape1, shape2);
	P(j,j+1) = 1 - P(j,j);
      }
    }

    double pdf_beta = log(prod(meanc(density_beta)));

    // 3. P
    Matrix<> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
      Matrix<> Sout = gaussian_ordinal_state_sampler_fixedsigma(stream, m, Y, X, beta_linear_st, Sigma, P);
      Matrix <double> s = Sout(_, 0);
      double shape1 = 0;
      double shape2 = 0;
      P(ns-1, ns-1) = 1;
      Matrix <double> P_addN(ns, 1);
      for (int j = 0; j<ns; ++j){
	for (int i = 0; i<N; ++i){
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

    // likelihood
    Matrix<> F = Matrix<>(N, ns);
    Matrix<> like(N, 1);
    Matrix<> pr1 = Matrix<>(ns, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);

    for (int t=0; t< N ; ++t){
      Matrix<> mu = X(t,_)*::t(beta_st);
      for (int j=0; j<ns ; ++j){
	int yt = Y[t];
	py[j]  =  oprobit_pdfLX(ncat, yt, mu[j], gamma_st(j,_));
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

    loglike = sum(log(like));

    // log prior ordinate
    Matrix<> density_beta_prior(ns, 1);
    Matrix<> density_P_prior(ns, 1);
    density_P[ns-1] = 1; //

    for (int j=0; j<ns ; ++j){
      density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv);
    }

    for (int j =0; j< (ns-1); ++j){
      density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1)));
    }
    double density_gamma_prior = ns*(ncat-2)*log(dunif(1, 0, 10));

    double logprior = sum(density_beta_prior) + sum(density_P_prior) + density_gamma_prior;
    logmarglike = (loglike + logprior) - (pdf_beta + pdf_P + pdf_gamma);

    if (verbose >0 ){
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("log_prior = %10.5f\n", logprior);
      Rprintf("log_beta = %10.5f\n", pdf_beta);
      Rprintf("log_P = %10.5f\n", pdf_P);
      Rprintf("log_gamma = %10.5f\n", pdf_gamma);
    }
  } // end of marginal likelihood

}//end

extern "C"{
  void cMCMCoprobitChange(double *betaout, double *betalinearout,
			 double *gammaout, double *Pout, double *psout, double *sout,
			 const double *Ydata,
			 const double *Xdata, const int *Xrow, const int *Xcol,
			 const int *m, const int *ncat,
			 const int *burnin, const int *mcmc, const int *thin, const int *verbose,
			 const double *tunedata, // const int *tdfdata,
			 const int *uselecuyer, const int *seedarray, const int *lecuyerstream,
			 const double *betastart,  const double *betalinearstart,
			 const double *gammastart, const double *Pstart,
			 const double *sigmastart, const double *a, const double *b,
			 const double *b0data, const double *B0data,
			 const double *A0data,
			 double *logmarglikeholder, double *loglikeholder,
			 const int *chib, const int *gammafixed)
  {

    // pull together Matrix objects
    const Matrix<> Y(*Xrow, 1, Ydata);
    const Matrix<> X(*Xrow, *Xcol, Xdata);
    const  int nstore = *mcmc / *thin;
    const  int N = *Xrow;
    const  int k = *Xcol;
    const  int gk = *ncat + 1;
    const  int ns = *m + 1;

    // generate starting values
    Matrix<> beta(ns, k, betastart);
    Matrix<> beta_linear(ns, k, betalinearstart);
    Matrix<> Sigma(1, 1, sigmastart);
    Matrix<> P(ns, ns, Pstart);
    Matrix<> b0(k, 1, b0data);
    Matrix<> B0(k, k, B0data);
    Matrix<> tune(ns, 1, tunedata);
    Matrix<> A0(ns, ns, A0data);
    double logmarglike;
    double loglike;

    // storage matrices
    Matrix<> beta_store(nstore, ns*k);
    Matrix<> beta_linear_store(nstore, ns*k);
    Matrix<> Z_store(nstore, N);
    Matrix<> P_store(nstore, ns*ns);
    Matrix<> ps_store(N, ns);
    Matrix<int> s_store(nstore, N);

    Matrix<> gamma(ns, gk, gammastart);
    Matrix<> gamma_store(nstore, ns*gk);
    MCMCPACK_PASSRNG2MODEL(MCMCoprobitChange_impl, *m, *ncat,
			   Y, X, beta, beta_linear, gamma, P,  Sigma,
			   b0, B0, A0,
			   *burnin, *mcmc, *thin, *verbose, tune,
			   *chib,  *gammafixed,
			   beta_store, beta_linear_store, gamma_store, Z_store,
			   P_store, ps_store, s_store,
			   logmarglike, loglike);

    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;

    // return output
    for (int i = 0; i<(nstore*ns*k); ++i){
      betaout[i] = beta_store[i];
      betalinearout[i] = beta_linear_store[i];
    }
    for (int i = 0; i<(nstore*ns*gk); ++i){
      gammaout[i] = gamma_store[i];
      }
    for (int i = 0; i<(nstore*ns*ns); ++i){
      Pout[i] = P_store[i];
    }
    for (int i = 0; i<(N*ns); ++i){
      psout[i] = ps_store[i];
    }
    for (int i = 0; i<(nstore*N); ++i){
      sout[i] = s_store[i];
    }
  }
}
#endif


