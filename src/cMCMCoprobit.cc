//////////////////////////////////////////////////////////////////////////
// cMCMCoprobit.cc is C++ code to estimate a ordinalprobit regression
//   model with a multivariate normal prior
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
// updated to the new version of Scythe 7/26/2004 KQ
// fixed a bug pointed out by Alexander Raach 1/16/2005 KQ
// updated to Scythe 1.0.X 7/10/2007 ADM
// Albert and Chib method added 9/20/2007 JHP
// fixed a bug pointed out by Shawn Treier 5/7/2018 KQ
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef CMCMCOPROBIT_CC
#define CMCMCOPROBIT_CC

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
#include "optimize.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;

static inline double lnmulttdens(const Matrix<>& theta,
				 const Matrix<>& mu,
				 const Matrix<>& C,
				 const double df){

  const int d = theta.size();
  //const Matrix<> z = t(theta - mu) * C;
  // C is now C' if VC mat is C C'
  const Matrix<> z = C * (theta - mu);
  double zsumsq = 0;
  for (int i=0; i<d; ++i){
    zsumsq += std::pow(z[i], 2);
  }

  return ( (-(df + d)/2) * std::log(1 + (1/df) * zsumsq)  );
}

// function that transforms alpha to gamma
Matrix<> gamma2alpha(const Matrix<>& gamma){
  const int m = gamma.rows() - 2;
  Matrix<> alpha(m, 1);
  alpha[0] = std::log(gamma[1]);
  for (int j=1; j< m ; ++j){
    alpha[j] = std::log(gamma[j+1] - gamma[j]);
  }
  return alpha;
}

// function that transforms alpha to gamma
Matrix<> alpha2gamma(const Matrix<>& alpha){
  const int m = alpha.rows();
  Matrix<> gamma(m+2, 1);
  gamma[0] = -300;
  gamma[m+1] = 300;
  for (int j=1; j<m+1 ; ++j){
    double gamma_sum = 0.0;
    for(int i=0;i<j; ++i){
      gamma_sum +=exp(alpha[i]);
    }
    gamma[j] = gamma_sum;
  }
  return gamma;
}

double lndmvn_jhp(const Matrix<double> &x, const Matrix<double> &mu,
		  const Matrix<double> &Sigma)
{
  int k = Sigma.cols();
  double first = (-k/2.0) * ::log(2 * M_PI) -0.5 * ::log(det(Sigma));
  Matrix< > second = ::t(x-mu)*invpd(Sigma)*(x-mu);
  return (first - 0.5*second[0]);
}



// orpobit_logpost
static inline double oprobit_logpost(const Matrix<double>& nY, const Matrix<double>& X,
				     const Matrix<double>& alpha,
				     const Matrix<double>& alpha_prior_mean,
				     const Matrix<double>& alpha_prior_var,
				     const Matrix<double>& beta){

  //  likelihood
  double loglike = 0.0;
  const int n = nY.rows();
  const int ncat = alpha.rows() + 1;

  // the linear predictor
  Matrix<> mu = X*beta;
  Matrix<> gamma = alpha2gamma(alpha);

  // compute prob
  Matrix<> cat_prob(n, ncat-1);
  //cat_prob: CATegory specific PROBability
  //the first col of cat_prob = pnorm(gamma_1 - mu)
  //thus, the col number is ncat - 1
  Matrix<> prob(n, ncat);
  for (int j=0; j< ncat-1; ++j){
    for (int i=0; i<n ; ++i){
      cat_prob(i, j) = pnorm(gamma[j+1] - mu[i], 0.0, 1.0);
    }
  }
  prob(_, ncat-1)  = 1 - cat_prob(_, ncat-2);    // top category
  prob(_, 0)       = cat_prob(_, 0);             // bottom category
  for (int j=1; j<(ncat-1); ++j){               // middle categories are actually differences of cumulatives
    prob(_, j) = cat_prob(_,j) - cat_prob(_, j-1);
  }

  // the loglikelihood
  for (int i=0; i<n; ++i){
    int ind = (int) nY[i];
    loglike += std::log(prob(i,ind-1));
  }

  // prior
  double logprior = 0.0;
  logprior = lndmvn_jhp(alpha, alpha_prior_mean, alpha_prior_var);

  return (loglike + logprior);

}// end of statis oprobit_logpost

// The oprobitModel class stores additional data to be passed to BFGS
class oprobitModel {
public:
  double operator() (const Matrix<double> alpha){
    const int n = y_.rows();
    const int ncat = alpha.rows() + 1;

    // the linear predictor
    Matrix<> mu = X_ * beta_;
    Matrix<> gamma = alpha2gamma(alpha);

    Matrix<> cat_prob(n, ncat-1);
    //cat_prob: CATegory specific PROBability
    //the first col of cat_prob = pnorm(gamma_1 - mu)
    //thus, the col number is ncat - 1
    Matrix<> prob(n, ncat);
    for (int j=0; j< ncat-1; ++j){
      for (int i=0; i<n ; ++i){
	cat_prob(i, j) = pnorm(gamma[j+1] - mu[i], 0.0, 1.0);
      }
    }
    prob(_, ncat-1) = 1 - cat_prob(_, ncat-2);    // top category
    prob(_, 0) = cat_prob(_, 0);               // bottom category
    for (int j=1; j<(ncat-1); ++j){               // middle categories are actually differences of cumulatives
      prob(_, j) = cat_prob(_,j) - cat_prob(_, j-1);
    }

    // the loglikelihood
    double loglike = 0.0;
    for (int i=0; i<n; ++i){
      int ind = y_[i];
      loglike += std::log(prob(i,ind-1));
    }
    return -1 * loglike;
  }
  Matrix<double> y_;
  Matrix<double> X_;
  Matrix<double> beta_;
};

/* MCMCoprobit implementation.  Takes Matrix<> reference which it
 * fills with the posterior.
 */
template <typename RNGTYPE>
void MCMCoprobit_impl (rng<RNGTYPE>& stream, const int * Y,
		       const Matrix<>& nY, const Matrix<>& X, Matrix<>& beta, Matrix<>& gamma,
		       const Matrix<>& b0, const Matrix<>& B0,
		       const Matrix<>& alpha_prior_mean, const Matrix<>& alpha_prior_var,
		       const unsigned int burnin, const unsigned int mcmc, const unsigned int thin,
		       const unsigned int verbose, const Matrix<>& tune, const double tdf,
		       const unsigned int cowles, Matrix<>& result) {

  // define constants and from cross-product matrices
  const unsigned int tot_iter = burnin + mcmc;  // total number of mcmc iterations
  const unsigned int nstore = mcmc / thin;      // number of draws to store
  const unsigned int k = X.cols();
  const unsigned int N = X.rows();
  const int ncat = gamma.rows() - 1;
  const Matrix<> XpX = crossprod(X);

  // storage matrix or matrices
  Matrix<> storemat(nstore, k+ncat+1);

  // initialize Z
  Matrix<> Z(N,1);
  Matrix<> Xbeta = X * beta;

  // Gibbs loop
  int count = 0;
  int accepts = 0;

  Matrix<> gamma_p = gamma;
  Matrix<> gamma_new = gamma + 1;
  Matrix<> alpha = gamma2alpha(gamma_new);
  Matrix<> alpha_hat = alpha;

  // initialize current value stuff
// JHP  Matrix<> propV = tune * alpha_prior_var * tune;
// JHP  Matrix<> propCinvT = t(cholesky(invpd(propV)));
// JHP  double logpost_cur = oprobit_logpost(nY, X, alpha, alpha_prior_mean, alpha_prior_var, beta);
// JHP  double logjump_cur = lnmulttdens(alpha_prior_mean, alpha_hat, propCinvT, tdf);

  double tune_double = tune[0];
  for (unsigned int iter = 0; iter < tot_iter; ++iter) {

    //////////////////
    if (cowles == 1){
      //////////////////
      // Cowles method [gamma | Z, beta]
      for (int i=2; i<(ncat); ++i){
	if (i==(ncat-1)){
	  gamma_p[i] = stream.rtbnorm_combo(gamma[i], ::pow(tune_double, 2.0),
					    gamma_p[i-1]);
	}
	else {
	  gamma_p[i] = stream.rtnorm_combo(gamma[i], ::pow(tune_double, 2.0),
					   gamma_p[i-1],
					   gamma[i+1]);
	}
      }
      double loglikerat = 0.0;
      double loggendenrat = 0.0;

      // loop over observations and construct the acceptance ratio
      for (unsigned int i=0; i<N; ++i){
	if (Y[i] == ncat){
	  loglikerat = loglikerat
	    + log(1.0  - pnorm(gamma_p[Y[i]-1] - Xbeta[i], 0.0, 1.0) )
	    - log(1.0 - pnorm(gamma[Y[i]-1] - Xbeta[i], 0.0, 1.0) );
	}
	else if (Y[i] == 1){
	  loglikerat = loglikerat + log(pnorm(gamma_p[Y[i]] - Xbeta[i], 0.0, 1.0)  )
	    - log(pnorm(gamma[Y[i]] - Xbeta[i], 0.0, 1.0) );
	}
	else{
	  loglikerat = loglikerat
	    + log(pnorm(gamma_p[Y[i]] - Xbeta[i], 0.0, 1.0) -
		  pnorm(gamma_p[Y[i]-1] - Xbeta[i], 0.0, 1.0) )
	    - log(pnorm(gamma[Y[i]] - Xbeta[i], 0.0, 1.0) -
		  pnorm(gamma[Y[i]-1] - Xbeta[i], 0.0, 1.0) );
	}
      }
      for (int j=2; j<ncat; ++j){	   
	loggendenrat = loggendenrat
	  + log(pnorm(gamma[j+1], gamma[j], tune[0]) - 
		pnorm(gamma_p[j-1], gamma[j], tune[0]) )  
	  - log(pnorm(gamma_p[j+1], gamma_p[j], tune[0]) - 
		pnorm(gamma[j-1], gamma_p[j], tune[0]) );
	
      }
      double logacceptrat = loglikerat + loggendenrat;
      if (stream.runif() <= exp(logacceptrat)){
	gamma = gamma_p;
	++accepts;
      }

    }// end of Cowles

    //////////////////
    else {
      //////////////////
      // Albert and Chib (2001) method
      // Step 1: [gamma | Z, beta]
      oprobitModel oprobit_model;
      oprobit_model.y_ = nY;
      oprobit_model.X_ = X;
      oprobit_model.beta_ = beta;

      // tailored proposal MH
      Matrix<> alpha_hat = BFGS(oprobit_model, alpha, stream, 100, 1e-5, false);
      Matrix<> alpha_V = invpd(hesscdif(oprobit_model, alpha_hat));
      //note that oprobit_model contains the multiplication by -1
      Matrix<> propV = tune * alpha_V * tune;
      Matrix<> propCinvT = ::t(cholesky(invpd(propV)));

      // Draw alpha_can from multivariate t
      Matrix<> alpha_can = alpha_hat + stream.rmvt(propV, tdf);

      // compute components of transition kernel
      double logpost_can = oprobit_logpost(nY, X, alpha_can, alpha_prior_mean, alpha_prior_var, beta);
      double logjump_can = lnmulttdens(alpha_can, alpha_hat, propCinvT, tdf);
      double logpost_cur = oprobit_logpost(nY, X, alpha, alpha_prior_mean, alpha_prior_var, beta);
      double logjump_cur = lnmulttdens(alpha, alpha_hat, propCinvT, tdf);
      double ratio = exp(logpost_can - logjump_can - logpost_cur + logjump_cur);
      const double u = stream();
      if (u < ratio) {
	alpha = alpha_can;
	gamma = alpha2gamma(alpha);
	logpost_cur = logpost_can;
	logjump_cur = logjump_can;
	++accepts;
      }
    }// end of AC method

    // Step 2: [Z| gamma, beta, y]
    Xbeta = X * beta;
    // Matrix<> Z_noconst(N,1);
    for (unsigned int i=0; i<N; ++i){
      Z[i] = stream.rtnorm_combo(Xbeta[i], 1.0, gamma[Y[i]-1], gamma[Y[i]]);
    }

    // Step 3: [beta|Z, gamma]
    const Matrix<> XpZ = t(X) * Z;
    beta = NormNormregress_beta_draw(XpX, XpZ, b0, B0, 1.0, stream);

    // store values in matrices
    if (iter >= burnin && ((iter % thin)==0)){
      for (unsigned int j=0; j<k; ++j)
	storemat(count, j) = beta[j];
      for (int j=0; j<(ncat+1); ++j)
	storemat(count, j+k) = gamma[j];
      ++count;
    }

    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n\nMCMCoprobit iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("beta = \n");
      Rprintf("%10.5f\n", beta[0]-gamma[1]);
      for (unsigned int j=1; j<k; ++j)
	Rprintf("%10.5f\n", beta[j]);
      Rprintf("Metropolis acceptance rate for gamma = %3.5f\n\n",
	      static_cast<double>(accepts)/static_cast<double>(iter+1));
    }

    R_CheckUserInterrupt(); // allow user interrupts
  }// end of MCMC
  if (verbose > 0){
    Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    Rprintf("The Metropolis acceptance rate for gamma was %3.5f",
	    static_cast<double>(accepts) / static_cast<double>(tot_iter));
    Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  }
  result = storemat;

}


extern "C"{

  void cMCMCoprobit(double *sampledata, const int *samplerow,
		   const int *samplecol, const int *Y,
		   const double *nYdata, const int *nYrow, const int *nYcol,
		   const double *Xdata,
		   const int *Xrow, const int *Xcol, const int *burnin,
		   const int *mcmc, const int *thin, const double *tunedata,
		   const int *tunerow, const int *tunecol,
		   const double* tdf,
		   const int *uselecuyer, const int *seedarray,
		   const int *lecuyerstream, const int *verbose,
		   const double *betadata, const int *betarow,
		   const int *betacol, const double* gammadata,
		   const int* gammarow, const int* gammacol,
		   const double *b0data, const int *b0row, const int *b0col,
		   const double *B0data, const int *B0row, const int *B0col,
		   const double *a0data, const int *a0row, const int *a0col,
		   const double *A0data, const int *A0row, const int *A0col,
		   const int *cowles) {

    // pull together Matrix objects
    const Matrix <> nY(*nYrow, *nYcol, nYdata);
    const Matrix <> X(*Xrow, *Xcol, Xdata);
    Matrix <> beta(*betarow, *betacol, betadata);
    Matrix <> gamma(*gammarow, *gammacol, gammadata);
    const Matrix <> b0(*b0row, *b0col, b0data);
    const Matrix <> B0(*B0row, *B0col, B0data);
    const Matrix <> alpha_prior_mean(*a0row, *a0col, a0data);
    const Matrix <> alpha_prior_prec(*A0row, *A0col, A0data);
    const Matrix <> alpha_prior_var = invpd(alpha_prior_prec);
    const Matrix<> tune(*tunerow, *tunecol, tunedata);

    Matrix<> storagematrix;
    MCMCPACK_PASSRNG2MODEL(MCMCoprobit_impl, Y, nY, X, beta, gamma, b0, B0,
			   alpha_prior_mean, alpha_prior_var,
			   *burnin, *mcmc, *thin, *verbose, tune, *tdf,
			   *cowles, storagematrix);

    const unsigned int size = *samplerow * *samplecol;
    for (unsigned int i = 0; i < size; ++i)
      sampledata[i] = storagematrix(i);
  }
}

#endif
