// MCMCmetrop1R.cc samples from a user-written posterior code in R using a
// random walk Metropolis algorithm
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
// KQ 6/24/2004
// updated to work with new Scythe and RNGs ADM 7/24/2004

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

extern "C" {

#include <Rdefines.h>
#include <Rinternals.h>


  // function that evaluatees the user supplied R function
  static double user_fun_eval(SEXP fun, SEXP theta, SEXP myframe){
    
    SEXP R_fcall;
    if(!isFunction(fun)) error("`fun' must be a function");
    if(!isEnvironment(myframe)) error("myframe must be an environment");
    PROTECT(R_fcall = lang2(fun, R_NilValue));
    SETCADR(R_fcall, theta);
    SEXP funval = eval(R_fcall, myframe);
    if (!isReal(funval)) error("`fun' must return a double");
    double fv = REAL(funval)[0];
    UNPROTECT(1);
    return fv;
  }


  // the function that actually does the sampling and returns a value to R
  SEXP MCMCmetrop1R_cc(SEXP fun, SEXP theta, SEXP myframe, SEXP burnin_R,
		       SEXP mcmc_R, SEXP thin_R, 
		       SEXP verbose, SEXP lecuyer_R, SEXP seedarray_R,
		       SEXP lecuyerstream_R, SEXP logfun, 
		       SEXP propvar_R){

    // put burnin_R, mcmc_R, and thin_R into ints
    const int burnin = INTEGER(burnin_R)[0];
    const int mcmc = INTEGER(mcmc_R)[0];
    const int thin = INTEGER(thin_R)[0];
    
    // put rng stuff together and initiate stream
    int seedarray[6];
    for(int i=0; i<6; ++i) seedarray[i] = INTEGER(seedarray_R)[i];
    rng *stream = MCMCpack_get_rng(INTEGER(lecuyer_R)[0],
       seedarray, INTEGER(lecuyerstream_R)[0]);    
       
    // put propvar_R into a Matrix
    double* propvar_data = REAL(propvar_R);
    const int propvar_nr = nrows(propvar_R);
    const int propvar_nc = ncols(propvar_R);
    Matrix <double> propvar (propvar_nc, propvar_nr, propvar_data);
    propvar = t(propvar);
    
    // define constants
    const int npar = length(theta);
    const int tot_iter = burnin + mcmc;
    const int nsamp = mcmc/thin;
    const Matrix <double> propc  = cholesky(propvar);
    
    // define matrix to hold the sample
    Matrix <double> sample (nsamp, npar);

    // put theta into a Scythe Matrix 
    double* theta_data = REAL(theta);
    const int theta_nr = length(theta);
    const int theta_nc = 1;
    Matrix <double> theta_M (theta_nc, theta_nr, theta_data);
    theta_M = t(theta_M);

    // evaluate userfun at starting value
    double userfun_cur =  user_fun_eval(fun, theta, myframe);
    if (INTEGER(logfun)[0]==0) userfun_cur = ::log(userfun_cur);
    
    
    // THE METROPOLIS SAMPLING
    int count = 0;
    int accepts = 0;
    for (int iter=0; iter<tot_iter; ++iter){

      // generate candidate value of theta 
      Matrix <double> theta_can_M = theta_M + propc * stream->rnorm(npar,1);
      
      // put theta_can_M into a SEXP
      SEXP theta_can;
      theta_can = PROTECT(allocVector(REALSXP, npar));
      for (int i=0; i<npar; ++i){
	REAL(theta_can)[i] = theta_can_M[i];
      }
      
      // evaluate user function fun at candidate theta
      double userfun_can = user_fun_eval(fun, theta_can, myframe);
      if (INTEGER(logfun)[0]==0) userfun_can = ::log(userfun_can);
      const double ratio = ::exp(userfun_can - userfun_cur);
      
      if (stream->runif() < ratio){
	theta = theta_can;
	theta_M = theta_can_M;
	userfun_cur = userfun_can;
	++accepts;
      }

      // store values in matrices
      if ((iter%thin)==0 && iter >= burnin){ 
	for (int j = 0; j < npar; j++)
	  sample(count, j) = REAL(theta)[j];
	++count;
      }
            
      if (iter % INTEGER(verbose)[0] == 0 && INTEGER(verbose)[0] > 0 ) {
	Rprintf("MCMCmetrop1R iteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("theta = \n");
	for (int i=0; i<npar; ++i)
	  Rprintf("%10.5f\n", REAL(theta)[i]);
	Rprintf("function value = %10.5f\n", userfun_cur);
	Rprintf("Metropolis acceptance rate = %3.5f\n\n", 
		static_cast<double>(accepts) / 
		static_cast<double>(iter+1));	
      } 

      
      UNPROTECT(1);      
      R_CheckUserInterrupt(); // allow user interrupts
    }

    // put the sample into a SEXP and return it   
    SEXP sample_SEXP;
    sample_SEXP = PROTECT(allocMatrix(REALSXP, nsamp, npar));
    for (int i=0; i<nsamp; ++i){
      for (int j=0; j<npar; ++j){
	REAL(sample_SEXP)[i + nsamp*j] = sample(i,j);
      }
    }
    UNPROTECT(1);

     delete stream; // clean up random number stream

    // print the the acceptance rate to the console in a way that 
    // everyone (even Windows users) can see
    Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    Rprintf("The Metropolis acceptance rate was %3.5f", 
	    static_cast<double>(accepts) / static_cast<double>(tot_iter));
    Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

    // return the sample
    return sample_SEXP;
    
  }
}
