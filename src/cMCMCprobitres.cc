//////////////////////////////////////////////////////////////////////////
// MCMCprobitres.cc is a program that simulates draws from the posterior
// density of a probit regression model and returns latent residuals.  
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
// updated to Scythe 1.0.X 7/28/2007 KQ
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef MCMCPROBITRES_CC
#define MCMCPROBITRES_CC

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

/* MCMCprobitres implementation.  Takes Matrix<> reference and fills with the
 * posterior.
 */
template <typename RNGTYPE>
void MCMCprobitres_impl (rng<RNGTYPE>& stream, const Matrix<>& Y,
			 const Matrix<>& X, Matrix<>& beta,
			 Matrix<>& resvec,
			 const Matrix<>& b0,
			 const Matrix<>& B0,  unsigned int burnin, 
			 unsigned int mcmc,
			 unsigned int thin, unsigned int verbose,  bool chib,  
			 Matrix<>& result, 
			 double& logmarglike) {
  
  // define constants and from cross-product matrices
  const unsigned int tot_iter = burnin + mcmc;  // total number of mcmc iterations
  const unsigned int nstore = mcmc / thin;      // number of draws to store
  const unsigned int k = X.cols();
  const unsigned int N = X.rows();
  const Matrix<> XpX = crossprod(X);
  const Matrix<> B0inv = invpd(B0);
  
  // storage matrix or matrices
  Matrix<> beta_store(nstore, k);
  Matrix<> Z_store(nstore, N);
   
  // initialize Z
  Matrix<> Z(N,1);
  
  // MCMC sampling starts here
  unsigned int count = 0;
  for (unsigned int iter = 0; iter < tot_iter; ++iter){
    
    // [Z| beta, y]
    const Matrix<> Z_mean = X * beta;
    for (unsigned int i=0; i<N; ++i){
      if (Y[i] == 1.0)
	Z[i] = stream.rtbnorm_combo(Z_mean[i], 1.0, 0);
      if (Y[i] == 0.0)
	Z[i] = stream.rtanorm_combo(Z_mean[i], 1.0, 0);
    }
    
    // [beta|Z]
    const Matrix<> XpZ = t(X) * Z;
    beta = NormNormregress_beta_draw(XpX, XpZ, b0, B0, 1.0, stream);
    
    // store values in matrices
    if (iter >= burnin && ((iter % thin)==0)){ 
      for (unsigned int j = 0; j < k; j++){
	beta_store(count, j) = beta[j];
      }
      Z_store(count,_) = Z;
      for (unsigned int j=0; j<(resvec.rows()); ++j){
	const int i = static_cast<int>(resvec[j]) - 1;
	beta_store(count, j+k) = Z[i] - Z_mean[i];
      }
      ++count;
    }
    
    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n\nMCMCprobit iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("beta = \n");
      for (unsigned int j=0; j<k; ++j)
	Rprintf("%10.5f\n", beta[j]);
    }
    
    R_CheckUserInterrupt(); // allow user interrupts   
    
  } // end MCMC loop
  
  
  if(chib==1){
    Rprintf("\n Marginal Likelihood Computation Starts!\n"); 

    Matrix<double> beta_star = meanc(beta_store);
    Matrix<double> density_beta(nstore, 1);
    for (unsigned int iter = 0; iter<nstore; ++iter){
      const Matrix<> Z_reduced = Z_store(iter,_);
      const Matrix<double> XpZ = (::t(X)*Z_reduced);
      const Matrix<double> Bn = invpd(B0inv + XpX);
      const Matrix<double> bn = Bn*gaxpy(B0inv, b0, XpZ);
      density_beta(iter) = exp(lndmvn(beta_star, bn, Bn));	
    }
    double logbeta = log(prod(meanc(density_beta)));
    
    double loglike = 0.0;
    Matrix<> eta = X * beta_star;
    for (unsigned int i = 0; i < X.rows(); ++i) {
      double phi = pnorm(eta(i), 0, 1);
      loglike += log(dbinom(Y(i), 1, phi));
    }
     
    // calculate log prior ordinate
    double logprior = lndmvn(beta_star, b0, B0inv);
    
    logmarglike = loglike + logprior - logbeta;
    
    Rprintf("\n logmarglike %10.5f", logmarglike, "\n"); 
    Rprintf("\n loglike %10.5f", loglike, "\n"); 
    
  }// end of marginal likelihood computation
 
  result = beta_store;  
}




extern "C"{
  
  void MCMCprobitres(double *sampledata, const int *samplerow, 
		     const int *samplecol, const double *Ydata, 
		     const int *Yrow, const int *Ycol, const double *Xdata, 
		     const int *Xrow, const int *Xcol,  
		     const double *resvecdata, const int *resvecrow,
		     const int *resveccol, const int *burnin,
		     const int *mcmc, const int *thin, const int *uselecuyer, 
		     const int *seedarray, const int *lecuyerstream, 
		     const int *verbose, const double *betastartdata, 
		     const int *betastartrow, const int *betastartcol, 
		     const double *b0data, const int *b0row, 
		     const int *b0col, const double *B0data, 
		     const int *B0row, const int *B0col, 
		     double *logmarglikeholder, // double *loglikeholder, 
		     const int *chib) {  
    
    
    // pull together Matrix objects
    const Matrix <> Y(*Yrow, *Ycol, Ydata);
    const Matrix <> X(*Xrow, *Xcol, Xdata);
    Matrix <> resvec(*resvecrow, *resveccol, 
			   resvecdata);
    Matrix <> beta(*betastartrow, *betastartcol, 
		   betastartdata);
    const Matrix <> b0(*b0row, *b0col, b0data);
    const Matrix <> B0(*B0row, *B0col, B0data);
    double logmarglike;
   
    Matrix<> storagematrix;
    MCMCPACK_PASSRNG2MODEL(MCMCprobitres_impl, Y, X, beta, resvec, 
			   b0, B0, *burnin,
			   *mcmc, *thin, *verbose, *chib, 
			   storagematrix, 
 			   logmarglike);
    
    // return output
    const unsigned int size = *samplerow * *samplecol;
    for (unsigned int i=0; i<size; ++i)
      sampledata[i] = storagematrix(i);
  }
  
}

#endif
