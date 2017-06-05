//////////////////////////////////////////////////////////////////////////
// cMCMCSVDreg.R samples from the posterior distribution of a Gaussian
// linear regression model in which the X matrix has been decomposed
// with an SVD. Useful for prediction when number of columns of X
// is (possibly much) greater than the number of rows of X.
//
// See West, Mike. 2000. "Bayesian Regression in the 'Large p, Small n'
//      Paradigm." Duke ISDS Discussion Paper 2000-22.
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
// KQ 9/9/2005
// KQ 8/1/2007 ported to Scythe 1.0.2
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef CMCMCSVDREG_CC
#define CMCMCSVDREG_CC

#include "MCMCrng.h"
#include "MCMCfcds.h"
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "rng.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;



template<typename RNGTYPE>
void MCMCSVDreg_impl(rng<RNGTYPE>& stream,
		     double *sampledata, const int *samplerow,
		     const int *samplecol,
		     const double *Ydata, const int *Yrow, const int *Ycol,
		     const int *Ymiss,
		     const double *Adata, const int *Arow, const int *Acol,
		     const double *Ddata, const int *Drow, const int *Dcol,
		     const double *Fdata, const int *Frow, const int *Fcol,
		     const int *burnin, const int *mcmc,
		     const int *thin, const int *uselecuyer,
		     const int *seedarray,
		     const int *lecuyerstream, const int *verbose,
		     const double *taustartdata, const int *taustartrow,
		     const int *taustartcol,
		     const double *g0data,
		     const int *g0row, const int *g0col,
		     const double *a0, const double *b0,
		     const double* c0, const double* d0,
		     const double* w0,
		     const int* betasamp
		     ){


    // pull together Matrix objects
    Matrix <> y(*Yrow, *Ycol, Ydata);
    Matrix <> A(*Arow, *Acol, Adata);
    Matrix <> D(*Drow, *Dcol, Ddata);
    Matrix <> F(*Frow, *Fcol, Fdata);
    Matrix <> g0(*g0row, *g0col, g0data);

    // define constants
    const int tot_iter = *burnin + *mcmc;
    // total number of mcmc iterations
    const int nstore = *mcmc / *thin; // number of draws to store
    const int k = *Arow;
    const int n = *Yrow;
    Matrix<> FtD = t(F) * D;
    Matrix<> DinvF = invpd(D) * F;
    Matrix<> Dg0 = D * g0;
    //double dsquared[n]; OLD (NEW BELOW)
    double* dsquared = new double[n];

    for (int i=0; i<n; ++i){
      dsquared[i] = std::pow(D(i,i), 2);
    }
    int holder = 0;
    for (int i=0; i<n; ++i){
      if (Ymiss[i] == 1){
	++holder;
      }
    }
    const int nYmiss = holder;

    // storage matrices
    Matrix<> beta_store;
    if (*betasamp == 1){
      beta_store = Matrix<double>(k, nstore);
    }
    Matrix<> gamma_store(n, nstore);
    Matrix<> Y_store(nYmiss, nstore);
    Matrix<> sigma2_store(1, nstore);
    Matrix<> tau2_store(n, nstore);


    // set starting values
    // double tau2[n]; OLD WAY (NEW BELOW)
    double* tau2 = new double[n];

    for (int i=0; i<n; ++i){
      tau2[i] = taustartdata[i];
    }
    Matrix <> gamma(n, 1);
    double sigma2 = 1;
    Matrix <double> beta;
    if (*betasamp == 1){
      beta = Matrix<double>(k, 1);
    }



    /////////////////// Gibbs sampler ///////////////////
    int count = 0;
    for (int iter = 0; iter < tot_iter; ++iter) {

      // sample [sigma2 | Y, A, D, F, tau2]
      Matrix <double> Fy = F * y;
      Fy = Fy - Dg0;
      double q = 0.0;
      for (int i=0; i<n; ++i){
	q += (std::pow(Fy[i], 2) / (1 + tau2[i]));
      }
      sigma2 = stream.rigamma((*a0+n)*0.5, (*b0 + q)*0.5);


      // sample [gamma | Y, A, D, F, sigma2, tau2]
      // w0[i] is assumed to be the prior prob that gamma[i] == 0
      Matrix <> gammahat = DinvF * y;
      for (int i=0; i<n; ++i){
	double mstar = (g0[i] + tau2[i]*gammahat[i]) / (1 + tau2[i]);
	double vstar = (sigma2 * tau2[i]) / (dsquared[i] * (tau2[i] + 1));
	//double wstar = 1.0 - (1.0 - w0[i]) /
	//  (1.0 - w0[i] + w0[i] * dnorm(0.0, mstar, std::sqrt(vstar)));

	Matrix<> gammanoti = gamma;
	gammanoti[i] = 0.0;

	Matrix<> residvec = y - FtD * gammanoti;
	Matrix<> margmeanvec = FtD(_,i) * g0[i];
	Matrix<> margVarmat = FtD(_,i) * t(FtD(_,i)) *
	  (sigma2 * tau2[i]) / (dsquared[i]) + eye<double>(n) * sigma2;

	double logw0 = std::log(w0[i]);
	double log1minusw0 = std::log(1.0 - w0[i]);
	double logf0 = 0.0;
	for (int j=0; j<n; ++j){
	  logf0 += lndnorm(residvec[j], 0.0, std::sqrt(sigma2));
	}
	double logfnot0 = lndmvn(residvec, margmeanvec, margVarmat);

	double logdenom;
	if ( (logw0 + logf0) > (log1minusw0 + logfnot0)){
	  logdenom = logw0 + logf0 + std::log(1.0 +
					      std::exp(log1minusw0 -
						       logw0 + logfnot0 -
						       logf0));
	}
	else{
	  logdenom = log1minusw0 + logfnot0 + std::log(1.0 +
						       std::exp(logw0 -
								log1minusw0 + logf0 -
								logfnot0));
	}
	double wstar = std::exp(logw0 + logf0 - logdenom);


	if (stream.runif() < wstar){
	  gamma[i] = 0.0;
	}
	else {
	  gamma[i] = stream.rnorm(mstar, std::sqrt(vstar));
	}
      }


      // sample [tau2 | Y, A, D, F, gamma, sigma2]
      for (int i=0; i<n; ++i){
	double gamg2 = std::pow((gamma[i] - g0[i]), 2);
	tau2[i] = stream.rigamma((1+c0[i])*0.5,
				  (gamg2 * (dsquared[i] / sigma2) + d0[i])
				  *0.5);
      }


      // sample [y[miss] | Y[obs], A, D, F, gamma, sigma2, tau2]
      Matrix <> eta = FtD * gamma;
      for (int i=0; i<n; ++i){
	if (Ymiss[i] == 1){
	  y[i] = stream.rnorm(eta[i], std::sqrt(sigma2));
	}
      }


      // optional sample [beta | Y, A, D, F, gamma, sigma2, tau2]
      if (*betasamp == 1){
	beta = A * gamma;
      }





      // store draws in storage matrix (or matrices)
      if (iter >= *burnin && (iter % *thin == 0)) {
	sigma2_store[count] = sigma2;
	int Ymiss_count = 0;
	for (int i = 0; i<n; ++i){
	  gamma_store(i, count) = gamma[i];
	  tau2_store(i, count) = tau2[i];
	  if (Ymiss[i] == 1){
	    Y_store(Ymiss_count, count) = y[i];
	    ++Ymiss_count;
	  }
	}
	if (*betasamp == 1){
	  for (int i=0; i<k; ++i){
	    beta_store(i, count) = beta[i];
	  }
	}
	++count;
      }


      // print output to stdout
      if(*verbose > 0 && iter % *verbose == 0) {
	Rprintf("\n\nMCMCSVDreg iteration %i of %i \n",
		(iter+1), tot_iter);
	Rprintf("gamma = \n");
	for (int j=0; j<n; ++j)
	  Rprintf("%10.5f\n", gamma[j]);
	Rprintf("tau2 = \n");
	for (int j=0; j<n; ++j)
	  Rprintf("%10.5f\n", tau2[j]);
	Rprintf("sigma2 = %10.5f\n", sigma2);
      }

      R_CheckUserInterrupt(); // allow user interrupts

    } // end MCMC loop


    // cleanup
    delete [] dsquared;
    delete [] tau2;


    // load draws into sample array
    Matrix <> storagematrix = cbind(t(Y_store), t(gamma_store));
    storagematrix = cbind(storagematrix, t(tau2_store));
    storagematrix = cbind(storagematrix, t(sigma2_store));
    if (*betasamp == 1){
      storagematrix = cbind(storagematrix, t(beta_store));
    }

    const int size = *samplerow * *samplecol;
    for(int i = 0; i < size; ++i)
      sampledata[i] = storagematrix[i];



}



extern "C" {

  // simulate from posterior distribution and return an mcmc by parameters
  // matrix of the posterior sample
  void cMCMCSVDreg(double *sampledata, const int *samplerow,
		  const int *samplecol,
		  const double *Ydata, const int *Yrow, const int *Ycol,
		  const int *Ymiss,
		  const double *Adata, const int *Arow, const int *Acol,
		  const double *Ddata, const int *Drow, const int *Dcol,
		  const double *Fdata, const int *Frow, const int *Fcol,
		  const int *burnin, const int *mcmc,
		  const int *thin, const int *uselecuyer,
		  const int *seedarray,
		  const int *lecuyerstream, const int *verbose,
		  const double *taustartdata, const int *taustartrow,
		  const int *taustartcol,
		  const double *g0data,
		  const int *g0row, const int *g0col,
		  const double *a0, const double *b0,
		  const double* c0, const double* d0,
		  const double* w0,
		  const int* betasamp) {



    MCMCPACK_PASSRNG2MODEL(MCMCSVDreg_impl,
			   sampledata, samplerow, samplecol,
			   Ydata, Yrow, Ycol,
			   Ymiss,
			   Adata, Arow, Acol,
			   Ddata, Drow, Dcol,
			   Fdata, Frow, Fcol,
			   burnin, mcmc, thin,
			   uselecuyer, seedarray, lecuyerstream, verbose,
			   taustartdata, taustartrow, taustartcol,
			   g0data, g0row, g0col,
			   a0, b0, c0, d0, w0,
			   betasamp);


  } // end cMCMCSVDreg
} // end extern "C"



#endif
