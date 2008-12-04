//////////////////////////////////////////////////////////////////////////
// MCMCmnl.h contains multinomial logit functions called by both the
// metropolis hastings and slice sampling implemenetation of MCMCmnl,
// such as a function that returns the log posterior.
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
// DBP 7/27/2007
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef MCMCMNL_H
#define MCMCMNL_H

#include "matrix.h"
#include "algorithm.h"
#include "distributions.h"
#include "la.h"
#include "smath.h"

using namespace std;
using namespace scythe;

inline double 
mnl_logpost(const Matrix<>& Y, const Matrix<>& X, const Matrix<>& beta,
	    const Matrix<>& beta_prior_mean, 
	    const Matrix<>& beta_prior_prec)
{
  
  //  likelihood
  double loglike = 0.0;
  const Matrix<double,Row> numera = exp(X * beta);
  //numer = reshape(numer, Y.rows(), Y.cols());
  //numer.resize(Y.rows(), Y.cols(), true);
  Matrix<double,Row> numer(Y.rows(), Y.cols(), false);
  copy<Row,Row>(numera, numer);
  double *denom = new double[Y.rows()];
  for (unsigned int i = 0; i < Y.rows(); ++i) {
    denom[i] = 0.0;
    for (unsigned int j = 0; j < Y.cols(); ++j) {
      if (Y(i,j) != -999){
	denom[i] += numer(i,j);
      }
    }
    for (unsigned int j = 0; j < Y.cols(); ++j) {
      if (Y(i,j) == 1.0){
	loglike += std::log(numer(i,j) / denom[i]);
      }
    }
  }
  
  delete [] denom;
  
  // prior
  //  double logprior = 0.0;
  //if (beta_prior_prec(0,0) != 0) {
  //  logprior = lndmvn(beta, beta_prior_mean, invpd(beta_prior_prec));
  // }
  //
  // the following is only up to proportionality
  const double logprior = -0.5 *(t(beta - beta_prior_mean) * 
                          beta_prior_prec * 
  			   (beta - beta_prior_mean))(0);
   
  

  return (loglike + logprior);
}

#endif
