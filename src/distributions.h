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
 * scythestat/distributions.h
 *
 * Provides definitions for PDFs, CDFs, and some common functions
 * (gamma, beta, etc).
 *
 */


#ifndef SCYTHE_DISTRIBUTIONS_H
#define SCYTHE_DISTRIBUTIONS_H

#include <cfloat>

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#else
#include "scythestat/matrix.h"
#endif

/* Fill in some defs from R that aren't in math.h */
#ifndef M_PI
#define M_PI 3.141592653589793238462643383280
#endif
#define M_LN_SQRT_2PI 0.918938533204672741780329736406
#define M_LN_SQRT_PId2  0.225791352644727432363097614947
#define M_1_SQRT_2PI  0.39894228040143267793994605993
#define M_2PI   6.28318530717958647692528676655
#define M_SQRT_32 5.656854249492380195206754896838

namespace SCYTHE {

  void throw_on_nonconv (bool);

  /*************
   * Functions *
   *************/

  /* The gamma function */
  double gammafn (const double &);

  /* The natural log of the absolute value of the gamma function */
  double lngammafn (const double &);

  /* The beta function */
  double betafn (const double &, const double &);

  /* The natrual log of the beta function */
  double lnbetafn(const double &, const double &);

  /* factorial */
  int factorial (const int &n);
  
  /* The natural log of the factorial */
  double lnfactorial (const int &n);

  /*********************************
   * Fully Specified Distributions *
   *********************************/
  
  /**** beta distribution ****/
  
  /* CDFs */
  double pbeta (const double &, const double &,
                const double &);

  Matrix<double> pbeta( const int &, const int &,
                        const double &, const double &,
                        const double &);
  /* PDFs */
  double dbeta (const double &, const double &,
                const double &);

  Matrix<double> dbeta( const int &, const int &,
                        const double &, const double &,
                        const double &b);

  /* Other */

  /* Returns the natural log of the ordinate of the Beta density
   * evaluated at x with Shape1 a, and Shape2 b
   */
  double lndbeta1(const double &, const double &,
                  const double &b);
  
  /**** binomial distribution ****/
  
	/* CDFs */
  double pbinom(const double &, const double &,const double &);

  Matrix<double> pbinom ( const int &, const int &,
                          const double &, const double &,
                          const double &);

  /* PDFs */
  double dbinom(const double &, const double &, const double &);

  Matrix<double> dbinom(const int &, const int &,
                        const double &, const double &,
                        const double &);


  /**** The Chi Square Distribution ****/
  
  /* CDFs */
  double pchisq(const double &, const double &);

  Matrix<double> pchisq(const int &, const int &,
                        const double &, const double &);
  /* PDFs */
  double dchisq(const double &, const double &);

  Matrix<double> dchisq(const int &, const int &,
                        const double &, const double &);
  
  /**** The Exponential Distribution ****/

  /* CDFs */
  double pexp(const double &, const double &);

  Matrix<double> pexp(const int &, const int &, const double &,
                      const double &);

  /* PDFs */
  double dexp(const double &, const double &);

  Matrix<double> dexp(const int &, const int &, const double &,
                      const double &);

  /**** The f Distribution ****/

  /* CDFs */
  double pf(const double &, const double &, const double &);

  Matrix<double> pf(const int &, const int &, const double &,
                    const double &, const double &);
  /* PDFs */
  double df(const double &, const double &, const double &);

  Matrix<double> df(const int &, const int &, const double &,
                    const double &, const double &);

  /**** The Gamma Distribution ****/

  /* CDFs */
  // note that this takes scale (or 1/rate for R-users)
  double pgamma (double, const double &, const double &);

  Matrix<double> pgamma(const int &, const int &, const double &, 
                        const double &scale = 1);

  /* PDFs */
  // note that this takes scale (or 1/rate for R-users)
  double dgamma(const double &, const double &,
                const double &scale = 1);

  Matrix<double> dgamma(const int &, const int &,
                        const double &, const double &scale = 1);

  /**** The Logistic Distribution ****/

  /* CDFs */
  double plogis(const double &, const double &location = 0.0,
                const double &scale = 1.0);

  Matrix<double> plogis(const int &, const int &,
                        const double &, const double &location = 0.0,
                        const double &scale = 1.0);

  /* PDFs */
  double dlogis(const double &, const double &location = 0.0,
                const double &scale = 1.0);

  Matrix<double> dlogis(const int &, const int &,
                        const double &, const double &location = 0.0,
                        const double &scale = 1.0);


  /**** The Log Normal Distribution ****/

  /* CDFs */
  double plnorm(const double& x, const double& logmean = 0.0,
                const double& logsd = 1.0);

  Matrix<double> plnorm(const int& rows, const int& cols,
                        const double& x, const double& logmean = 0.0,
                        const double& logsd = 1.0);

  /* PDFs */
  double dlnorm(const double& x, const double& logmean = 0.0,
                const double& logsd = 1.0);

  Matrix<double> dlnorm(const int& rows, const int& cols,
                        const double& x, const double& logmean = 0.0,
                        const double& logsd = 1.0);

  /**** The Negative Binomial Distribution ****/

  /* CDFs */
  double pnbinom(const double &, const double &, const double &);

  Matrix<double> pnbinom( const int &, const int &,
                          const double &, const double &,
                          const double &);

  /* PDFs */
  double dnbinom(const double &, const double &, const double &);

  Matrix<double> dnbinom( const int &, const int &,
                          const double &, const double &,
                          const double &);

  /**** The Normal Distribution ****/
  
  /* CDFs */
  double pnorm (const double &x, const double &mu = 0.0,
                const double &sigma = 1.0);

  /* PDFs */
  double dnorm( const double& x, const double& mu = 0.0,
                const double& sigma = 1.0);

  Matrix<double> dnorm( const int& rows, const int& cols,
                        const double& x, const double& mu = 0.0,
                        const double& sigma = 1.0);
  /* Other */
  
  /* Returns the univariate standard normal cumulative distribution
   * function (CDF)
   */
  double pnorm2(const double &x, const bool &lower_tail,
                const bool &log_p);

  void pnorm_both(double x, double *cum, double *ccum, int i_tail,
                  bool log_p);

  /* Returns the quantile of the standard normal distribution 
   * associated with a given probability p
   */
  double qnorm1 (const double& in_p);

  /* Returns the log of the univariate normal density ordinate 
   * evaluated  at x
   */
  double lndnorm (const double& x, const double& mu = 0.0,
                  const double& sigma = 1.0);

  /**** The Poison Distribution ****/
  
  /* CDFs */

  double ppois(const double &, const double &);

  Matrix<double> ppois( const int &, const int &,
                        const double &, const double &);
  /* PDFs */
  double dpois(const int &, const double &);
  
  Matrix<double> dpois( const int &, const int &,
                        const double &, const double &);

  /**** The t Distribution ****/

  /* CDFs */

  double pt(const double &, const double &);

  Matrix<double> pt(const int &, const int &, const double &,
                    const double &);
  /* PDFs */
  double dt(const double &, const double &);

  Matrix<double> dt(const int &, const int &, const double &,
                    const double &);
  
  /* Others */

  /* Returns the univariate Student-t density evaluated at x 
   * with mean mu, scale sigma^2, and nu degrees of freedom
   */
  double dt1( const double &, const double &, const double &, 
              const double &);
  
  /* Returns the natural log of the univariate Student-t density 
   * evaluated at x with mean mu, scale sigma^2, and nu 
   * degrees of freedom
   */
  double lndt1(const double &, const double &,
                const double &, const double &);

  /**** The Uniform Distribution ****/

  /* CDFs */
  double punif (const double &, const double &a = 0.0,
                const double &b = 1.0);
  
  Matrix<double> punif( const int &, const int &,
                        const double &, const double &a = 0.0,
                        const double &b = 1.0);

  /* PDFs */
  double dunif( const double &, const double &a = 0.0,
                const double &b = 1.0);

  Matrix<double> dunif( const int &, const int &,
                        const double &, const double &a = 0.0,
                        const double &b = 1.0);

  /**** The Weibull Distribution ****/
  
  /* CDFs */
  double pweibull(const double &, const double &,
                  const double &scale = 1.0);

  Matrix<double> pweibull(const int &, const int &,
                          const double &, const double &,
                          const double &scale = 1.0);

  /* PDFs */
  double dweibull(const double &, const double &, 
                  const double &scale = 1.0);

  Matrix<double> dweibull(const int &, const int &,
                          const double &, const double &,
                          const double &scale = 1.0);

  /************************************
   * Partially Finished Distributions *
   ************************************/
  
  /* Multivariate Normal */
  double lndmvn ( const Matrix<double> &, const Matrix<double> &, 
                  const Matrix<double> &);


  /********************
  * Helper Functions *
  ********************/
  namespace INTERNAL {
  
    /* Evaluates an Chebyshev series at a given point */
    double chebyshev_eval (const double &, const double *, const int &);

    /* Computes the log gamma correction factor for x >= 10 */
    double lngammacor (const double &);

    /* Helper for dpois and dgamma */
    double dpois_raw (const double &, const double &);

    /* Evaluates the "deviance part" */
    double bd0 (const double &, const double &);

    /* Computes the log of the error term in Stirling's formula */
    double stirlerr (const double &);

    double pbeta_raw(const double &, const double &, const double &);

    double dbinom_raw(const double &, const double &, const double &,
                      const double &);
  }

}

#endif
