/* 
 * Scythe Statistical Library
 * Copyright (C) 2000-2002 Andrew D. Martin and Kevin M. Quinn;
 * 2002-2004 Andrew D. Martin, Kevin M. Quinn, and Daniel
 * Pemstein.  All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * under the terms of the GNU General Public License as published by
 * Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.  See the text files COPYING
 * and LICENSE, distributed with this source code, for further
 * information.
 * --------------------------------------------------------------------
 * scythestat/rng.h
 *
 * Provides the implementation for the rng class.  This abstract
 * class forms the foundation of random number generation in Scythe.
 *
 */

#ifndef SCYTHE_RNG_CC
#define SCYTHE_RNG_CC

#include <iostream>
#include <cmath>

#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#ifdef SCYTHE_COMPILE_DIRECT
#include "rng.h"
#include "distributions.h"
#include "util.h"
#include "ide.h"
#include "stat.h"
#else
#include "scythestat/rng.h"
#include "scythestat/distributions.h"
#include "scythestat/util.h"
#include "scythestat/ide.h"
#include "scythestat/stat.h"
#endif

namespace SCYTHE {
  
  /* Default constructor */
  rng::rng ()
  {
  }
  
  rng::~rng()
  {
  }
  
  /* Random Numbers */
  Matrix<double>
  rng::runif (const int &rows, const int &cols)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Attempted to create Matrix of size <= 0");
    }
    
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = runif();
    
    return temp;
  }
  
  double
  rng::rbeta (const double& alpha, const double& beta)
  {
    static double report;
    double xalpha, xbeta;
    
    // Check for allowable parameters
    if (alpha <= 0) {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
			       __LINE__, "alpha <= 0");
    }
    if (beta <= 0) {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
			       __LINE__, "beta <= 0");
    }
    
    xalpha = rchisq (2 * alpha);
    xbeta = rchisq (2 * beta);
    report = xalpha / (xalpha + xbeta);
    
    return (report);
  }
  
  Matrix<double>
  rng::rbeta (const int& rows, const int& cols, const double& alpha,
	      const double& beta)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Attempted to create Matrix of size <= 0");
    }
    
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rbeta (alpha, beta);
    
    return temp;
  }
  
  /* Return a pseudo-random deviate from a non-cental hypergeometric
   * distribution
   */
  double 
  rng::rnchypgeom(const double& m1, const double& n1, 
		  const double& n2, const double& psi, 
		  const double& delta)
  {
    // Calculate mode of mass function
    double a = psi - 1;
    double b = -1 * ((n1+m1+2)*psi + n2 - m1);
    double c = psi * (n1+1) * (m1+1);
    double q = -0.5 * ( b + sgn(b)* ::sqrt(::pow(b,2) - 4*a*c));
    double root1 = c/q;
    double root2 = q/a;
    double el = max(0.0, m1-n2);
    double u = min(n1,m1);
    double mode = floor(root1);
    int exactcheck = 0;
    if (u<mode || mode<el) {
      mode = floor(root2);
      exactcheck = 1;
    }
 

    int size = static_cast<int>(u+1);

    double *fvec = new double[size];
    fvec[static_cast<int>(mode)] = 1.0;
    double s;
    // compute the mass function at y
    if (delta <= 0 || exactcheck==1){  //exact evaluation 
      // sum from mode to u
      double f = 1.0;
      s = 1.0;
      for (double i=(mode+1); i<=u; ++i){
        double r = ((n1-i+1)*(m1-i+1))/(i*(n2-m1+i)) * psi;
        f = f*r;
        s += f;
        fvec[static_cast<int>(i)] = f;
      }
     
      // sum from mode to el
      f = 1.0;
      for (double i=(mode-1); i>=el; --i){
        double r = ((n1-i)*(m1-i))/((i+1)*(n2-m1+i+1)) * psi;
        f = f/r;
        s += f;
        fvec[static_cast<int>(i)] = f;
      }
    } else { // approximation
    	double epsilon = delta/10.0;
      // sum from mode to ustar
      double f = 1.0;
      s = 1.0;
      double i = mode+1;
      double r;
      do {
        if (i>u) break;
        r = ((n1-i+1)*(m1-i+1))/(i*(n2-m1+i)) * psi;
        f = f*r;
        s += f;
        fvec[static_cast<int>(i)] = f;
        ++i;
      } while(f>=epsilon || r>=5.0/6.0);
     
      // sum from mode to elstar
      f = 1.0;
      i = mode-1;
      do {
        if (i<el) break;
        r = ((n1-i)*(m1-i))/((i+1)*(n2-m1+i+1)) * psi;
        f = f/r;
        s += f;
        fvec[static_cast<int>(i)] = f;
        --i;
      } while(f>=epsilon || r <=6.0/5.0);         
    }

    double udraw = runif();
    double psum = fvec[static_cast<int>(mode)]/s;
    if (udraw<=psum)
      return mode;
    double lower = mode-1;
    double upper = mode+1;

    do{
      double fl;
      double fu;
      if (lower >= el)
        fl = fvec[static_cast<int>(lower)];
      else 
        fl = 0.0;

      if (upper <= u)
        fu = fvec[static_cast<int>(upper)];
      else
        fu = 0.0;

      if (fl > fu) {
        psum += fl/s;
        if (udraw<=psum)
          return lower;
        --lower;
      } else {
        psum += fu/s;
        if (udraw<=psum)
          return upper;
        ++upper;
      }
    } while(udraw>psum);
   
    delete [] fvec;
    exit(500000);
  }
	
  Matrix<double>
  rng::rnchypgeom(const int &rows, const int &cols, const double &m1,
		  const double &n1, const double &n2, 
		  const double &psi, const double &delta)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Attempted to create Matrix of size <= 0");
    }
    
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rnchypgeom (m1, n1, n2, psi, delta);
    
    return temp;
  }
  
  /* Random Numbers */
  int
  rng::rbinom (const int& n, const double& p)
  {
    static int report;
    int count = 0;
    double hold;
      
    // Check for allowable parameters
    if (n <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "n <= 0");
    }
    if (p < 0 || p > 1) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "p not in [0,1]");
    }
      
    // Loop and count successes
    for (int i = 0; i < n; i++) {
      hold = runif ();
      if (hold < p)
        count++;
    }
    report = count;
    
    return (report);
  }

  Matrix<double>
  rng::rbinom (const int& rows, const int& cols, const int& n,
    					const double& p)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rbinom (n, p);
    
    return temp;
  }

  double
  rng::rchisq (const double &nu)
  {
    static double report;
      
    // Check for allowable paramter
    if (nu <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Degrees of freedom <= 0");
    }
  
    // Return Gamma(nu/2, 1/2) deviate
    report = rgamma (nu / 2, .5);
    
    return (report);
  }

  Matrix<double>
  rng::rchisq (const int &rows, const int &cols, const double &nu)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rchisq (nu);
    
    return temp;
  }

  double
  rng::rexp (const double &beta)
  {
    static double report;
    
    // Check for allowable parameter
    if (beta <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Inverse scale parameter beta <= 0");
    }
    
    report = -std::log (runif ()) / beta;
    
    return (report);
  }

  Matrix<double>
  rng::rexp (const int &rows, const int &cols, const double &beta)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rexp (beta);
    
    return temp;
  }

  double
  rng::rf (const double &n1, const double &n2)
  {
    if (n1 <= 0.0 || n2 <= 0.0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "n1 or n2 <= 0");

    return ((rchisq(n1) / n1) / (rchisq(n2) / n2));
  }

  Matrix<double>
  rng::rf(const int &rows, const int &cols, const double &n1,
					const double &n2)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rf (n1, n2);
    
    return temp;
  }

  double
  rng::rgamma (const double &alpha, const double &beta)
  {
    static double report;

    // Check for allowable parameters
    if (alpha <= 0) {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
             __LINE__, "alpha <= 0");
    }
    if (beta <= 0) {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
             __LINE__, "beta <= 0");
    }

    if (alpha > 1)
      report = rgamma1 (alpha) / beta;
    else if (alpha == 1)
      report = -::log (runif ()) / beta;
    else if (alpha < 1)
      report = rgamma1 (alpha + 1) * ::pow (runif (), 1 / alpha) / beta;

    return (report);
  }

  Matrix<double>
  rng::rgamma(const int& rows, const int& cols, const double& alpha, 
    					const double& beta)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rgamma (alpha, beta);
  
    return temp;
  }

  double
  rng::rlogis (const double& alpha, const double& beta)
  {
    static double report;
    double unif;
      
    // Check for allowable paramters
    if (beta <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "beta <= 0");
    }
    
    unif = runif ();
    report = alpha + beta * std::log (unif / (1 - unif));
    
    return (report);
  }

  Matrix<double>
  rng::rlogis(const int& rows, const int& cols, const double& alpha, 
    					const double& beta)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rlogis (alpha, beta);
    
    return temp;
  }

  double
  rng::rlnorm (const double &logmean, const double &logsd)
  {
    if (logsd < 0.0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "standard deviation < 0");

    return std::exp(rnorm(logmean, logsd));
  }
  
  Matrix<double>
  rng::rlnorm(const int &rows, const int &cols, const double &logmean,
    					const double &logsd)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rlnorm (logmean, logsd);

    return temp;
  }

  double
  rng::rnbinom (const double &n, const double &p)
  {
    if (n <= 0 || p <= 0 || p > 1)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "n <= 0, p <= 0, or p > 1");

    return rpois(rgamma(n, (1 - p) / p));
  }

  Matrix<double>
  rng::rnbinom (const int &rows, const int &cols, const double &n,
    						const double &p)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = rnbinom(n, p);
    
    return temp;
  }

  double
  rng::rnorm (const double &mu, const double &sigma)
  {
    if (sigma <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Negative variance");
    }
    
    return (mu + rnorm1 () * sigma);
  }

  Matrix<double>
  rng::rnorm (const int &rows, const int &cols, const double &mu,
    					const double &sigma)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rnorm (mu, sigma);

    return temp;
  }

  int
  rng::rpois(const double &lambda)
  {
    if (lambda <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "lambda <= 0");
    }
    int n;
    
    if (lambda < 33){
      double cutoff = ::exp(-lambda);
      n = -1;
      double t = 1.0;
      do {
        ++n;
        t *= runif();
      } while (t > cutoff);    
    } else{
      int accept = 0;
      double c = 0.767 - 3.36/lambda;
      double beta = M_PI/::sqrt(3*lambda);
      double alpha = lambda*beta;
      double k = ::log(c) - lambda - ::log(beta);
        
      while (accept == 0){
        double u1 = runif();
        double x = (alpha - ::log((1-u1)/u1))/beta;
        while (x <= -0.5){
          u1 = runif();
          x = (alpha - ::log((1-u1)/u1))/beta;
        } 
        n = static_cast<int>(x + 0.5);
        double u2 = runif();
        double lhs = alpha - beta*x +
          ::log(u2/::pow(1+::exp(alpha-beta*x),2));
        double rhs = k + n*::log(lambda) - lnfactorial(n);
        if (lhs <= rhs)
          accept = 1;
      }
    }
    
    return n;
  }

  Matrix<int>
  rng::rpois (const int &rows, const int &cols, const double &lambda)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = rpois(lambda);
    
    return temp;
  }

  double
  rng::rt (const double& mu, const double& sigma2,
      const double& nu)
  {
    static double report;
    double x, z;
      
    // Check for allowable paramters
    if (sigma2 <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Variance parameter sigma2 <= 0");
    }
    if (nu <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "D.O.F parameter nu <= 0");
    }
    
    z = rnorm1 ();
    x = rchisq (nu);
    report = mu + ::sqrt (sigma2) * z * ::sqrt (nu) / ::sqrt (x);
    
    return (report);
  }

  Matrix<double>
  rng::rt (const int& rows, const int& cols, const double& mu,
      const double& sigma2, const double& nu)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rt (mu, sigma2, nu);
    
    return temp;
  }


  double
  rng::rweibull (const double &shape, const double &scale)
  {
    if (shape <= 0 || scale <= 0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "shape or scale <= 0");

    return scale * std::pow(-std::log(runif()), 1.0 / shape);
  }

  Matrix<double>
  rng::rweibull(const int& rows, const int& cols, const double& shape,
      					const double& scale)
  {
    int size = rows * cols;
    if (size <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = rweibull(shape,scale);
    
    return temp;
  }

  double
  rng::richisq (const double &nu)
  {
    static double report;
      
    // Check for allowable parameter
    if (nu <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Degrees of freedom <= 0");
    }
      
    // Return Inverse-Gamma(nu/2, 1/2) deviate
    report = rigamma (nu / 2, .5);
    return (report);
  }

  Matrix<double>
  rng::richisq (const int& rows, const int& cols, const double& nu)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = richisq (nu);
    
    return temp;
  }


  double
  rng::rigamma (const double& alpha, const double& beta)
  {
    static double report;
    
    // Check for allowable parameters
    if (alpha <= 0) {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
             __LINE__, "alpha <= 0");
    }
    if (beta <= 0) {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
             __LINE__, "beta <= 0");
    }
    // Return reciprocal of gamma deviate
    report = ::pow (rgamma (alpha, beta), -1);

    return (report);
  }

  Matrix<double>
  rng::rigamma (  const int &rows, const int &cols, const double &alpha,
    const double &beta)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rigamma (alpha, beta);
    
    return temp;
  }

  Matrix<double>
  rng::rwish(const int &v, const Matrix<double> &S)
  {
    if (! S.isSquare()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "S not square");
    }
    if (v < S.rows()) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "v < S.rows()");
    }
      
    Matrix<double> A(S.rows(), S.rows()); //XXX init to 0?
    Matrix<double> C = cholesky(S);
    Matrix<double> alpha;
      
    for (int i = 0; i < v; ++i) {
      alpha = C * rnorm(S.rows(), 1);
      A = A + (alpha * (!alpha));
    }

    return(A);
  }

  /* Dirichlet generator */
  Matrix<double> 
	rng::rdirich(const Matrix<double> &alpha) 
	{ 
    // Check for allowable parameters
    if (min(alpha) <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "alpha has elements < 0");
    }
    if (alpha.cols() > 1) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "alpha not column vector");
    }     
 
    int dim = alpha.rows();
    Matrix<double> y(dim, 1);
    double ysum = 0;
    for (int i=0; i<dim; ++i){
      y[i] = rgamma(alpha[i], 1);
      ysum += y[i];
    }

    Matrix<double> report = y;
    for (int i=0; i<dim; ++i)
      report[i] = y[i]/ysum;
    return(report);
  }
   

  /* Multivariate Normal */
  Matrix<double>
  rng::rmvnorm(const Matrix<double> &mu, const Matrix<double> &sigma)
	{  
    int dim = mu.rows();
    if (mu.cols() != 1) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "mu not column vector");
    }
    if (! sigma.isSquare()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "sigma not square");
    }
    if (sigma.rows() != dim) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "mu and sigma not conformable");
    }       
            
    Matrix<double> A = mu + cholesky(sigma) * rnorm(dim,1);
    return(A);
  }

  /* Multivariate t */
  Matrix<double>
  rng::rmvt (const Matrix<double> &sigma, const double &nu) {
    Matrix<double> result;
    if (nu <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "D.O.F parameter nu <= 0");
    }
    result = rmvnorm(Matrix<double>(sigma.rows(), 1, true, 0), sigma);
    return result / std::sqrt(rchisq(nu) / nu);
  }

  /* Bernoulli */
  int
  rng::rbern (const double &p)
  {
    static int report;
    double unif;
      
    // Check for allowable paramters
    if (p < 0 || p > 1) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "p parameter not in[0,1]");
    }
    
    unif = runif ();
    if (unif < p)
      report = 1;
    else
      report = 0;
    
    return (report);
  }

  Matrix<double>
  rng::rbern (const int& rows, const int& cols, const double& p)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Attempted to create Matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; ++i)
      temp[i] = rbern (p);
    
    return temp;
  }

  double 
  rng::rtnorm(const double& m, const double& v, const double& below, 
	      const double& above)
  {  
    if (below > above) {
      std::cout << "mean = " << m << " and var = " << v << std::endl << 
	"below = " << below << "  and above = " << above << std::endl;
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, \
				"Truncation bound not logically consistent");
    }
		
    double FA = 0.0;
    double FB = 0.0;
    double sd = ::sqrt(v);
    if ((::fabs((above-m)/sd) < 8.2) && (::fabs((below-m)/sd) < 8.2)){
      FA = pnorm2((above-m)/sd, true, false);
      FB = pnorm2((below-m)/sd, true, false);
    }
    if ((((above-m)/sd) < 8.2)  && (((below-m)/sd) <= -8.2) ){ 
      FA = pnorm2((above-m)/sd, true, false);
      FB = 0.0;
    }
    if ( (((above-m)/sd) >= 8.2)  && (((below-m)/sd) > -8.2) ){ 
      FA = 1.0;
      FB = FB = pnorm2((below-m)/sd, true, false);
    } 
    if ( (((above-m)/sd) >= 8.2) && (((below-m)/sd) <= -8.2)){
      FA = 1.0;
      FB = 0.0;
    }
    double term = runif()*(FA-FB)+FB;
    if (term < 5.6e-17)
      term = 5.6e-17;
    if (term > (1 - 5.6e-17))
      term = 1 - 5.6e-17;
    double draw = m + sd * qnorm1(term);
    if (draw > above)
      draw = above;
    if (draw < below)
      draw = below;
 		 
    return draw;
  }

  double
  rng::rtnorm_combo(const double& m, const double& v,
		    const double& below, const double& above)
  {
    if (below >= above) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Truncation bound not logically consistent");
    }
    double s = ::sqrt(v);

    if (( ((above-m)/s > 0.5) && ((m-below)/s > 0.5)) ||
    	( ((above-m)/s > 2.0) && ((below-m)/s < 0.25)) ||
    	( ((m-below)/s > 2.0) && ((above-m)/s > -0.25)) ){ 
      double x = rnorm(m, s);
      while ((x > above) || (x < below))
	x = rnorm(m,s);
      return x;
    } else {
      // use the inverse cdf method
      double FA = 0.0;
      double FB = 0.0;
      if ((::fabs((above-m)/s) < 8.2) && (::fabs((below-m)/s) < 8.2)){
	FA = pnorm2((above-m)/s, true, false);
	FB = pnorm2((below-m)/s, true, false);
      }
      if ((((above-m)/s) < 8.2)  && (((below-m)/s) <= -8.2) ){ 
	FA = pnorm2((above-m)/s, true, false);
	FB = 0.0;
      }
      if ( (((above-m)/s) >= 8.2)  && (((below-m)/s) > -8.2) ){ 
	FA = 1.0;
	FB = FB = pnorm2((below-m)/s, true, false);
      } 
      if ( (((above-m)/s) >= 8.2) && (((below-m)/s) <= -8.2)){
	FA = 1.0;
	FB = 0.0;
      }
      double term = runif()*(FA-FB)+FB;
      if (term < 5.6e-17)
	term = 5.6e-17;
      if (term > (1 - 5.6e-17))
	term = 1 - 5.6e-17;
      double x = m + s * qnorm1(term);
      if (x > above)
	x = above;
      if (x < below)
	x = below;
      return x;
    }    
  }

  double
  rng::rtbnorm_slice (const double& m, const double& v,
		      const double& below, const int& iter)
  {
    if (below < m) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Truncation point < mean");
    }
    if (v <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Variance non-positive");
    }
 		 
    double z = 0;
    double x = below + .00001;
 		 
    for (int i=0; i<iter; ++i){
      z = runif()*::exp(-1*::pow((x-m),2)/(2*v));
      x = runif()*( (m + ::sqrt(-2*v*::log(z))) - below) + below;
    }
    if (! finite(x)) {
      std::cerr << "WARNING in "
		<< __FILE__ << ", " << __PRETTY_FUNCTION__ << ", "
		<< __LINE__ << ": Mean extremely far from truncation point. "
		<< "Returning truncation point" << std::endl;
      return below; 
    }
    return x;
  }

  double
  rng::rtanorm_slice (const double& m, const double& v,
		      const double& above, const int& iter)
  {
    if (above > m) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Truncation point > mean");
    }
    if (v <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Variance non-positive");
    }
  
    double below = -1*above;
    double newmu = -1*m;
    double z = 0;
    double x = below + .00001;
 		 
    for (int i=0; i<iter; ++i){
      z = runif()*::exp(-1*::pow((x-newmu),2)/(2*v));
      x = runif()*( (newmu + ::sqrt(-2*v*::log(z))) - below) + below;
    }
    if (! finite(x)) {
      std::cerr << "WARNING in "
		<< __FILE__ << ", " << __PRETTY_FUNCTION__ << ", "
		<< __LINE__ << ": Mean extremely far from truncation point. "
		<< "Returning truncation point" << std::endl;
      return above; 
    }
		
    return -1*x;
  }

  double
  rng::rtbnorm_combo (const double& m, const double& v,
		      const double& below, const int& iter)
  {
    if (v <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Variance non-positive");
    }
    
    double s = ::sqrt(v);
    // do rejection sampling and return value
    //if (m >= below){
    if ((m/s - below/s ) > -0.5){
      double x = rnorm(m, s);
      while (x < below)
	x = rnorm(m,s);
      return x; 
    } else if ((m/s - below/s ) > -5.0 ){
      // use the inverse cdf method
      double above =  std::numeric_limits<double>::infinity();
      double x = rtnorm(m, v, below, above);
      return x;
    } else {
      // do slice sampling and return value
      double z = 0;
      double x = below + .00001;
      for (int i=0; i<iter; ++i){
	z = runif()*::exp(-1*::pow((x-m),2)/(2*v));
	x = runif()*( (m + ::sqrt(-2*v*::log(z))) - below) + below;
      }
    if (! finite(x)) {
	std::cerr << "WARNING in "
		  << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", "
		  << __LINE__ << ": Mean extremely far from truncation point. "
		  << "Returning truncation point" << std::endl;
	return below; 
      }
      return x;
    }
  }

  double
  rng::rtanorm_combo (const double& m, const double& v,
		      const double& above, const int& iter)
  {
    if (v <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
				__LINE__, "Variance non-positive");
    }
    double s = ::sqrt(v);
    // do rejection sampling and return value
    if ((m/s - above/s ) < 0.5){ 
      double x = rnorm(m, s);
      while (x > above)
	x = rnorm(m,s);
      return x;
    } else if ((m/s - above/s ) < 5.0 ){
      // use the inverse cdf method
      double below =  -std::numeric_limits<double>::infinity();
      double x = rtnorm(m, v, below, above);
      return x;
    } else {
      // do slice sampling and return value
      double below = -1*above;
      double newmu = -1*m;
      double z = 0;
      double x = below + .00001;
   			 
      for (int i=0; i<iter; ++i){
	z = runif()*::exp(-1*::pow((x-newmu),2)/(2*v));
	x = runif()*( (newmu + ::sqrt(-2*v*::log(z))) - below) + below;
      }
    if (! finite(x)) {
	std::cerr << "WARNING in "
		  << __FILE__ << ", " << __PRETTY_FUNCTION__ << ", "
		  << __LINE__ << ": Mean extremely far from truncation point. "
		  << "Returning truncation point" << std::endl;
	return above; 
      }
      return -1*x;
    }
  }
  
	double
  rng::rnorm1 ()
  {
    static int rnorm_count = 1;
    static double x2;
    double nu1, nu2, rsquared, sqrt_term;
    if (rnorm_count == 1){ // odd numbered passses
      do {
        nu1 = -1 +2*runif();
        nu2 = -1 +2*runif();
        rsquared = ::pow(nu1,2) + ::pow(nu2,2);
      } while (rsquared >= 1 || rsquared == 0.0);
      sqrt_term = ::sqrt(-2*::log(rsquared)/rsquared);
      x2 = nu2*sqrt_term;
      rnorm_count = 2;
      return nu1*sqrt_term;
    } else { // even numbered passes
      rnorm_count = 1;
      return x2;
    } 
  }

  double
  rng::rgamma1 (const double &alpha)
  {
    int test;
    double u, v, w, x, y, z, b, c;
    static double accept;

    // Check for allowable parameters
    if (alpha <= 1) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "alpha < 1");
    }

    // Implement Best's (1978) simulator
    b = alpha - 1;
    c = 3 * alpha - 0.75;
    test = 0;
    while (test == 0) {
      u = runif ();
      v = runif ();

      w = u * (1 - u);
      y = ::sqrt (c / w) * (u - .5);
      x = b + y;

      if (x > 0) {
        z = 64 * ::pow (v, 2) * ::pow (w, 3);
        if (z <= (1 - (2 * ::pow (y, 2) / x))) {
          test = 1;
          accept = x;
        } else if ((2 * (b * ::log (x / b) - y)) >= ::log (z)) {
          test = 1;
          accept = x;
        } else {
          test = 0;
        }
      }
    }
    
    return (accept);
  }
}

#endif /* SCYTHE_RNG_CC */
