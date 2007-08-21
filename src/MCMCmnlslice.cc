// MCMCmnlslice.cc DESCRIPTION HERE
//
// The initial version of this file was generated by the
// auto.Scythe.call() function in the MCMCpack R package
// written by:
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
// Copyright (C) 2004 Andrew D. Martin and Kevin M. Quinn
// 
// This file was initially generated on Wed Dec 29 15:27:40 2004
// REVISION HISTORY
//
// 7/28/07 DBP ported to scythe 1.0

#ifndef MCMCMNLSLICE_CC
#define MCMCMNLSLICE_CC

#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"
#include "MCMCmnl.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;

// eventually all of the slice sampling functions should be made more 
// general and put in MCMCfcds.{h cc}
//
// Radford Neal's (2000) doubling procedure coded for a logdensity
template <typename RNGTYPE>
static void
doubling(double (*logfun)(const Matrix<>&, const Matrix<>&,
			  const Matrix<>&, const Matrix<>&,
			  const Matrix<>&), const Matrix<>& beta, 
	 int index, double z, double w, int p, const Matrix<>& Y,
	 const Matrix<>& X, const Matrix<>& beta_prior_mean,
	 const Matrix<>& beta_prior_prec, rng<RNGTYPE>& stream, 
	 double& L, double& R){
  
  const double U = stream();
  const double x0 = beta(index);
  Matrix<> beta_L = beta;
  Matrix<> beta_R = beta;
  L = x0 - w * U;
  beta_L(index) = L;
  R = L + w;
  beta_R(index) = R;
  int K = p;
  while (K > 0 && 
	 (z < logfun(Y, X, beta_L, beta_prior_mean, beta_prior_prec) ||
	  z < logfun(Y, X, beta_R, beta_prior_mean, beta_prior_prec))) {
    double V = stream();
    if (V < 0.5){
      L = L - (R - L);
      beta_L(index) = L;
    }
    else {
      R = R + (R - L);
      beta_R(index) = R;
    }
    --K;
  }  
}

// Radford Neal's (2000) Accept procedure coded for a logdensity
static const bool 
Accept(double (*logfun)(const Matrix<>&, const Matrix<>&, 
			const Matrix<>&, const Matrix<>&, 
			const Matrix<>&), const Matrix<>& beta,
       int index, double x0, double z, double w, const Matrix<>& Y,
       const Matrix<>& X, const Matrix<>& beta_prior_mean,
       const Matrix<>& beta_prior_prec, double L, double R)
{
  double Lhat = L;
  double Rhat = R;
  bool D = false;
  while ((Rhat - Lhat ) > 1.1 * w) {
    double M = (Lhat + Rhat) / 2.0;
    if ( (x0 < M && beta(index) >= M) || (x0 >= M && beta(index) < M)){
      D = true;
    }
    if (beta(index) < M){
      Rhat = M;
    }
    else {
      Lhat = M;
    }

    Matrix<> beta_L = beta;
    Matrix<> beta_R = beta;
    beta_L[index] = Lhat;
    beta_R[index] = Rhat;

    if (D && z >= logfun(Y, X, beta_L, beta_prior_mean, beta_prior_prec)
	&& 
	z >=  logfun(Y, X, beta_R, beta_prior_mean, beta_prior_prec)) {
      return(false);
    }    
  }

  return(true);
}


// Radford Neal's (2000) shrinkage procedure coded for a log density
template <typename RNGTYPE>
static double
shrinkage(double (*logfun)(const Matrix<>&, const Matrix<>&, 
			   const Matrix<>&, const Matrix<>&,
			   const Matrix<>&), const Matrix<>& beta,
	  int index, double z, double w, const Matrix<>& Y,
	  const Matrix<>& X, const Matrix<>& beta_prior_mean,
	  const Matrix<>& beta_prior_prec, rng<RNGTYPE>& stream,
	  double L, double R)
{
  double Lbar = L;
  double Rbar = R;
  Matrix<> beta_x1 = beta;
  const double x0 = beta[index]; 
  for (;;) {
    const double U = stream();
    const double x1 = Lbar + U*(Rbar - Lbar);
    beta_x1(index) = x1;
    if (z <= logfun(Y, X, beta_x1, beta_prior_mean, beta_prior_prec) &&
	Accept(logfun, beta_x1, index, x0, z, w, Y, X, beta_prior_mean,
	       beta_prior_prec, L, R)) {
      return(x1);
    }
    if (x1 < x0) {
      Lbar = x1;
    } else {
      Rbar = x1;
    }
  } // end infinite loop
}

template <typename RNGTYPE>
void MCMCmnlslice_impl(rng<RNGTYPE>& stream, const Matrix<>& Y,
		       const Matrix<>& X, const Matrix<>& b0,
		       const Matrix<>& B0, const Matrix<>& V,
		       Matrix<>& beta, unsigned int burnin,
		       unsigned int mcmc, unsigned int thin,
		       unsigned int verbose, Matrix<>& storemat){
  // DEFINE constants
  const unsigned int tot_iter = burnin + mcmc;  // total iterations
  const unsigned int nstore = mcmc / thin;      // # of draws to store
  const unsigned int k = X.cols();
  
  // Initialize storage matrix
  storemat = Matrix<>(nstore, k, false);
  
  // proposal parameters
  const Matrix<> propV = invpd(B0 + invpd(V));
  const Matrix<> w_init = ones<double>(k, 1);
  for (unsigned int i = 0; i < k; ++i)
    w_init(i) = sqrt(propV(i,i)) *0.05;
  
  // starting values
  double L = -1.0;
  double R = 1.0;
  
  const unsigned int warmup_iter = 100;
  const unsigned int warmup_burnin = 10;
  const unsigned int p_init = 15;
  const Matrix<> widthmat(warmup_iter - warmup_burnin, k);
  // warm up sampling to choose the slice sampling parameters
  for (unsigned int iter = 0; iter < warmup_iter; ++iter) {
    for (unsigned int index = 0; index < k; ++index) {
      double funval = mnl_logpost(Y, X, beta, b0, B0);
      double z = funval - stream.rexp(1.0);
      doubling(&mnl_logpost, beta, index, z, w_init[index], p_init, Y,
	       X, b0, B0, stream, L, R);
      beta(index) = shrinkage(&mnl_logpost, beta, index, z, 
			      w_init(index), Y, X, b0, B0, stream, L, R);
      if (iter >= warmup_burnin)
	widthmat(iter-warmup_burnin, index) =  R - L;
    }
  }
  const Matrix<> w = meanc(widthmat);
  Matrix<int> p = ones<int>(k,1);
  for (unsigned int index = 0; index < k; ++index) {
    int p_temp = 2;
    while ((w(index) * pow(2.0, p_temp) ) < max(widthmat(_,index))) { 
      ++p_temp;
    } 
    p(index) = p_temp + 1;       
  }
	 
  unsigned int count = 0;
  ///// REAL MCMC SAMPLING OCCURS IN THIS FOR LOOP
  for(unsigned int iter = 0; iter < tot_iter; ++iter) {
    for (unsigned int index = 0; index < k; ++index) {
      double funval = mnl_logpost(Y, X, beta, b0, B0);
      double z = funval - stream.rexp(1.0);
      doubling(&mnl_logpost, beta, index, z, w(index), p(index), Y, X,
	       b0, B0, stream, L, R);
      beta(index) = shrinkage(&mnl_logpost, beta, index, z, w(index),
			      Y, X, b0, B0, stream, L, R);
    }

    // store draws in storage matrix (or matrices)
    if(iter >= burnin && (iter % thin == 0)) {
      for (unsigned int j = 0; j < k; j++)
	storemat(count, j) = beta(j);
      ++count;
    }

    // print output to stdout
    if(verbose > 0 && iter % verbose == 0) {
      Rprintf("\n\nMCMCmnl slice iteration %i of %i \n", (iter+1),
	      tot_iter);
      Rprintf("beta = \n");
      for (unsigned int j = 0; j < k; ++j)
	Rprintf("%10.5f\n", beta[j]);
    }

    R_CheckUserInterrupt(); // allow user interrupts
		 
  } // end MCMC loop
}

extern "C" {

  // MCMC sampling for MNL model via slice sampling
  void MCMCmnlslice(double *sampledata, const int *samplerow, 
		    const int *samplecol, const double *Ydata, 
		    const int *Yrow, const int *Ycol, 
		    const double *Xdata, const int *Xrow, 
		    const int *Xcol, const int *burnin, 
		    const int *mcmc, const int *thin, 
		    const int *uselecuyer, 
		    const int *seedarray, const int *lecuyerstream, 
		    const int *verbose, const double *betastartdata, 
		    const int *betastartrow, const int *betastartcol,
		    const double *b0data, const int *b0row, 
		    const int *b0col, const double *B0data, 
		    const int *B0row, const int *B0col,
		    const double *Vdata, const int *Vrow, 
		    const int *Vcol) {
     
    // pull together Matrix objects
    // REMEMBER TO ACCESS PASSED ints AND doubles PROPERLY
    const Matrix<> Y(*Yrow, *Ycol, Ydata);
    const Matrix<> X(*Xrow, *Xcol, Xdata);
    Matrix<> beta(*betastartrow, *betastartcol, betastartdata);     
    const Matrix<> b0(*b0row, *b0col, b0data);
    const Matrix<> B0(*B0row, *B0col, B0data);
    const Matrix<> V(*Vrow, *Vcol, Vdata);
 
    // storage matrix or matrices
    Matrix<double> storemat;
    MCMCPACK_PASSRNG2MODEL(MCMCmnlslice_impl, Y, X, b0, 
			   B0, V, beta,
			   *burnin, *mcmc, *thin, 
			   *verbose, storemat);

    // load draws into sample array
    for(unsigned int i = 0; i < storemat.size(); ++i)
      sampledata[i] = storemat(i);

  } // end MCMCmnlslice 
} // end extern "C"

#endif
