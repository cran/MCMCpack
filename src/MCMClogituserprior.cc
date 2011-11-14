//////////////////////////////////////////////////////////////////////////
// MCMClogituserprior.cc samples from the posterior distribution of a 
// logistic regression model with user-written prior coded in R using a
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
// KQ 8/17/2005 (based on current version of MCMCmetrop1R.cc)
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////



#ifndef MCMCLOGITUSERPRIOR_CC
#define MCMCLOGITUSERPRIOR_CC

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

#include <Rdefines.h>
#include <Rinternals.h>


using namespace std;
using namespace scythe;



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

static double logit_loglike(const Matrix<double>& Y, 
			    const Matrix<double>& X, 
			    const Matrix<double>& beta){
    
  // likelihood
  const Matrix<double> eta = X * beta;
  const Matrix<double> p = 1.0/(1.0 + exp(-eta));
  double loglike = 0.0;
  for (unsigned int i=0; i<Y.rows(); ++i)
    loglike += Y[i]*::log(p[i]) + (1-Y[i])*::log(1-p[i]);

  return loglike;
}

template <typename RNGTYPE>
void MCMClogituserprior_impl(rng<RNGTYPE>& stream, SEXP fun, 
			     SEXP theta, SEXP myframe,
			     unsigned int burnin,
			     unsigned int mcmc, unsigned int thin,
			     unsigned int verbose, bool logfun,
			     const Matrix<>& propvar, const Matrix<>& Y_M, 
			     const Matrix<>& X_M, SEXP& sample_SEXP){

    // define constants
    const unsigned int npar = length(theta);
    const unsigned int tot_iter = burnin + mcmc;
    const unsigned int nsamp = mcmc/thin;
    const Matrix <> propc  = cholesky(propvar);
    
    // define matrix to hold the sample
    Matrix <> sample (nsamp, npar, false);

    // put theta into a Scythe Matrix 
    double* theta_data = REAL(theta);
    const int theta_nr = length(theta);
    const int theta_nc = 1;
    Matrix <> theta_M (theta_nc, theta_nr, theta_data);
    theta_M = t(theta_M);

    // evaluate userfun at starting value
    double loglike_val = logit_loglike(Y_M, X_M, theta_M);
    double userfun_cur =  user_fun_eval(fun, theta, myframe);
    if (! logfun) 
      userfun_cur = ::log(userfun_cur);
    userfun_cur += loglike_val;
    

    // THE METROPOLIS SAMPLING
    int count = 0;
    int accepts = 0;
    for (unsigned int iter=0; iter<tot_iter; ++iter){
      
      // generate candidate value of theta 
      Matrix <> theta_can_M = theta_M + propc * stream.rnorm(npar,1, 0, 1);
      
      // put theta_can_M into a SEXP
      SEXP theta_can;
      theta_can = PROTECT(allocVector(REALSXP, npar));
      for (unsigned int i=0; i<npar; ++i){
	REAL(theta_can)[i] = theta_can_M(i);
      }
      
      // evaluate user function fun at candidate theta
      loglike_val = logit_loglike(Y_M, X_M, theta_can_M);
      double userfun_can = user_fun_eval(fun, theta_can, myframe);
      if (! logfun) 
	userfun_can = std::log(userfun_can);
      userfun_can += loglike_val;
      
      const double ratio = std::exp(userfun_can - userfun_cur);
      
      if (stream() < ratio){
	theta = theta_can;
	theta_M = theta_can_M;
	userfun_cur = userfun_can;
	++accepts;
      }

      // store values in matrices
      if ((iter%thin)==0 && iter >= burnin){ 
	for (unsigned int j = 0; j < npar; j++)
	  sample(count, j) = REAL(theta)[j];
	++count;
      }
            
      if (verbose && iter % verbose == 0) {
	Rprintf("MCMClogit iteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("beta = \n");
	for (unsigned int i=0; i<npar; ++i)
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
    sample_SEXP = PROTECT(allocMatrix(REALSXP, nsamp, npar));
    for (unsigned int i=0; i<nsamp; ++i){
      for (unsigned int j=0; j<npar; ++j){
	REAL(sample_SEXP)[i + nsamp*j] = sample(i,j);
      }
    }
    UNPROTECT(1);


    // print the the acceptance rate to the console in a way that 
    // everyone (even Windows users) can see
    Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    Rprintf("The Metropolis acceptance rate was %3.5f", 
	    static_cast<double>(accepts) / static_cast<double>(tot_iter));
    Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");


}








extern "C" {  

  // the function that actually does the sampling and returns a value to R
  SEXP MCMClogituserprior_cc(SEXP fun, SEXP Y_R, SEXP X_R, 
			     SEXP theta, SEXP myframe, SEXP burnin_R,
			     SEXP mcmc_R, SEXP thin_R, 
			     SEXP verbose, SEXP lecuyer_R, SEXP seedarray_R,
			     SEXP lecuyerstream_R, SEXP logfun, 
			     SEXP propvar_R){
    
    
    // put rng stuff together and initiate stream
    int seedarray[6];
    for(int i=0; i<6; ++i) seedarray[i] = INTEGER(seedarray_R)[i];
    int uselecuyer_cc = INTEGER(lecuyer_R)[0];
    int lecuyerstream_cc = INTEGER(lecuyerstream_R)[0];
    int* uselecuyer = &uselecuyer_cc;
    int* lecuyerstream = &lecuyerstream_cc;

       
    // put propvar_R into a Matrix
    double* propvar_data = REAL(propvar_R);
    const int propvar_nr = nrows(propvar_R);
    const int propvar_nc = ncols(propvar_R);
    Matrix <double,Row> propvarpre (propvar_nc, propvar_nr, propvar_data);
    //propvar = Matrix<>(propvar.rows(), propvar.cols(),
    //    propvar.begin_f<Row>());
    Matrix<> propvar = t(propvarpre);

    // put Y_R into a Scythe Matrix
    int* Y_data = INTEGER(Y_R);
    const int Y_nr = length(Y_R);
    const int Y_nc = 1;
    Matrix <int,Row> Y_Mpre (Y_nc, Y_nr, Y_data);
    //Y_M = Matrix<int>(Y_M.rows(), Y_M.cols(), Y_M.begin_f<Row>());
    Matrix<int> Y_M = t(Y_Mpre);

    // put X_R into a Scythe Matrix
    double* X_data = REAL(X_R);
    const int X_nr = nrows(X_R);
    const int X_nc = ncols(X_R);
    Matrix <double,Row> X_Mpre (X_nc, X_nr, X_data);
    //X_M = Matrix<>(X_M.rows(), X_M.cols(), X_M.begin_f<Row>());
    Matrix<> X_M = t(X_Mpre);


    SEXP sample_SEXP;
    MCMCPACK_PASSRNG2MODEL(MCMClogituserprior_impl, fun, theta, myframe, 
			   INTEGER(burnin_R)[0], INTEGER(mcmc_R)[0], 
			   INTEGER(thin_R)[0],
			   INTEGER(verbose)[0], INTEGER(logfun)[0], 
			   propvar, Y_M, X_M, 
			   sample_SEXP)
    

    // return the sample
    return sample_SEXP;
    
  }
}



#endif
