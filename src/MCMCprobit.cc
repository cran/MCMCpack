// MCMCprobit.cc is C++ code to estimate a probit regression model with
//   a multivariate normal prior
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

#include <iostream>
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace SCYTHE;
using namespace std;

extern "C"{
  
  void MCMCprobit(double *sampledata, const int *samplerow, 
		  const int *samplecol, const double *Ydata, 
		  const int *Yrow, const int *Ycol, const double *Xdata, 
		  const int *Xrow, const int *Xcol, const int *burnin, 
		  const int *mcmc, const int *thin, const int *lecuyer, 
		  const int *seedarray, const int *lecuyerstream, 
		  const int *verbose, const double *betastartdata, 
		  const int *betastartrow, const int *betastartcol, 
		  const double *b0data, const int *b0row, const int *b0col, 
		  const double *B0data, const int *B0row, const int *B0col) {  
    
    // pull together Matrix objects
    const Matrix <double> Y = r2scythe(*Yrow, *Ycol, Ydata);
    const Matrix <double> X = r2scythe(*Xrow, *Xcol, Xdata);
    Matrix <double> beta = r2scythe(*betastartrow, *betastartcol, 
				    betastartdata);
    const Matrix <double> b0 = r2scythe(*b0row, *b0col, b0data);
    const Matrix <double> B0 = r2scythe(*B0row, *B0col, B0data);

    // define constants and from cross-product matrices
    const int tot_iter = *burnin + *mcmc;  // total number of mcmc iterations
    const int nstore = *mcmc / *thin;      // number of draws to store
    const int k = X.cols();
    const int N = X.rows();
    const Matrix<double> XpX = crossprod(X);
    
    // storage matrix or matrices
    Matrix<double> storemat(nstore, k);
    
    // initialize rng stream
    rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
    
    // initialize Z
    Matrix<double> Z(N,1);

    // MCMC sampling starts here
    int count = 0;
    for (int iter = 0; iter < tot_iter; ++iter)
      {

	// [Z| beta, y]
	const Matrix<double> Z_mean = X * beta;
	for (int i=0; i<N; ++i){
	  if (Y[i] == 1.0)
	    Z[i] = stream->rtbnorm_combo(Z_mean[i], 1.0, 0);
	  if (Y[i] == 0.0)
	    Z[i] = stream->rtanorm_combo(Z_mean[i], 1.0, 0);
	}
	
	// [beta|Z]
	const Matrix<double> XpZ = t(X) * Z;
	beta = NormNormregress_beta_draw(XpX, XpZ, b0, B0, 1.0, stream);
	
	// store values in matrices
	if (iter >= *burnin && ((iter % *thin)==0)){ 
	  for (int j = 0; j < k; ++j)
	    storemat(count, j) = beta[j];
	  ++count;
	}
	
	// print output to stdout
	if(*verbose > 0 && iter % *verbose == 0){
	  Rprintf("\n\nMCMCprobit iteration %i of %i \n", (iter+1), tot_iter);
	  Rprintf("beta = \n");
	  for (int j=0; j<k; ++j)
	    Rprintf("%10.5f\n", beta[j]);
	}
	
	R_CheckUserInterrupt(); // allow user interrupts    	
	
      } // end MCMC loop
     delete stream; // clean up random number stream

    // return output
    const int size = *samplerow * *samplecol;
    for (int i=0; i<size; ++i)
      sampledata[i] = storemat[i];
    
  }
  
}
