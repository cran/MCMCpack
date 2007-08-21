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
// Copyright (C) 2004 Andrew D. Martin and Kevin M. Quinn
// 
// updated to the new version of Scythe 7/26/2004 KQ
// updated to Scythe 1.0.X 7/28/2007 KQ

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
			 unsigned int thin, unsigned int verbose,  
			 Matrix<>& result) {
  
  // define constants and from cross-product matrices
  const unsigned int tot_iter = burnin + mcmc;  // total number of mcmc iterations
  const unsigned int nstore = mcmc / thin;      // number of draws to store
  const unsigned int k = X.cols();
  const unsigned int N = X.rows();
  const Matrix<> XpX = crossprod(X);
  
  // holding matrices
  Matrix<> storemat(nstore, k+resvec.rows());
  
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
      for (unsigned int j = 0; j < k; j++)
	storemat(count, j) = beta[j];
      for (unsigned int j=0; j<(resvec.rows()); ++j){
	const int i = static_cast<int>(resvec[j]) - 1;
	storemat(count, j+k) = Z[i] - Z_mean[i];
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
  
  result = storemat;
  
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
		     const int *B0row, const int *B0col) {  
    
    // pull together Matrix objects
    const Matrix <> Y(*Yrow, *Ycol, Ydata);
    const Matrix <> X(*Xrow, *Xcol, Xdata);
    Matrix <> resvec(*resvecrow, *resveccol, 
			   resvecdata);
    Matrix <> beta(*betastartrow, *betastartcol, 
		   betastartdata);
    const Matrix <> b0(*b0row, *b0col, b0data);
    const Matrix <> B0(*B0row, *B0col, B0data);

    Matrix<> storagematrix;
    MCMCPACK_PASSRNG2MODEL(MCMCprobitres_impl, Y, X, beta, resvec, 
			   b0, B0, *burnin,
			   *mcmc, *thin, *verbose, 
			   storagematrix);
    
    // return output
    const unsigned int size = *samplerow * *samplecol;
    for (unsigned int i=0; i<size; ++i)
      sampledata[i] = storagematrix(i);
  }
  
}

#endif
