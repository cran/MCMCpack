/* Scythe_Simulate.h
 *
 * This header provides definitions for random number generators,
 * PDFs, CDFs, and some common functions (gamma, beta, etc) for the
 * Scythe Statistical Library.
 *
 * Scythe Statistical Library
 * Copyright (C) 2003, Andrew D. Martin, Kevin M. Quinn, and Daniel
 * Pemstein.  All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.  A copy of this license is included
 * with this library (LICENSE.GPL).
 *
 * This library utilizes code from a number of other open source
 * projects.  Specific copyright information is provided with the
 * applicable code.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 * USA.
 *
 * This code written by:
 *
 * Andrew D. Martin
 * Assistant Professor
 * Deptartment of Political Science
 * Campus Box 1063
 * Washington University
 * One Brookings Drive
 * St. Louis, MO 63130
 * admartin@artsci.wustl.edu
 *
 * Kevin M. Quinn
 * Assistant Professor
 * Department of Government and
 * Center for Basic Research in the Social Sciences
 * 34 Kirkland Street
 * Harvard University
 * Cambridge, MA 02138
 * kquinn@fas.harvard.edu
 *
 * Daniel Pemstein
 * Deptartment of Poltical Science
 * 702 South Wright Street
 * University of Illinois at Urbana-Champaign
 * Urbana, IL 61801
 * dbp@uiuc.edu
 */


#ifndef SCYTHE_SIMULATE_H
#define SCYTHE_SIMULATE_H

#include "Scythe_Matrix.h"
#include <cfloat>

namespace SCYTHE {

  /*******************
   * Underlying RNGs *
   *******************/
  
  /* Marsaglia's SWB + KISS */
  // XXX get cite

  /* Set the random seeds */
  void set_ranmars_seed(const unsigned long &z, const unsigned long &w,
                        const unsigned long &jsr,
                        const unsigned long &jcong);

  /* Generates a random uniform number */
  double ranmars();

  /* Mersenne Twister by Nishimura and Matsumoto:  Default seed is
   * 5489UL
   *
   * mersenne_* code based on:
   * A C-program for MT19937, with initialization improved 2002/1/26.
   * Coded by Takuji Nishimura and Makoto Matsumoto.
   * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   * All rights reserved.
   *
   * The bulk of the mersenne code is in Scythe_Simluate.cc,
   * along with a full copy of the license (BSD) under which the
   * original code was distributed.
   */
  
  void set_mersenne_seed(const unsigned long &s);

  /* Generate a random unsigned 32bit integer */
  unsigned long mersenne_int32();

  /* Generate a random double on [0,1] */
  inline double mersenne_real1()
  {
    return mersenne_int32()*(1.0/4294967295.0);
  }

  /* Generate a random double on [0,1) */
  inline double mersenne_real2()
  {
    return mersenne_int32()*(1.0/4294967296.0);
  }

  /* Generate a random double on (0,1) */
  inline double mersenne_real3()
  {
    return (((double)mersenne_int32()) + 0.5)*(1.0/4294967296.0);
  }

  /* "Static" var specifying rng to use */
#ifdef __MINGW32__
  static double (*rng)() = mersenne_real3;
#else
  namespace {
    double (*rng)() = mersenne_real3;
  }
#endif

  /* Set the generator to use.  Default is mersenne_real3 */
  inline void set_rng (double (*fun)())
  {
    rng = fun;
  }

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
  
  /* Random Numbers */
  double rbeta (const double &, const double &);

  Matrix<double> rbeta (const int &, const int &,
                        const double &, const double &);
  
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
  
  /* returns a pseudo-random deviate from a non-central hypergeometric
   * distribution
   * References: Uses algorithm of Liang and Rosen. 2001.
   * "Fast and Stable Algorithms for Computing and Sampling From the
   * Noncentral Hypergeometric Distribution", The American Statistician,
   * 55: 366-369.
   * 
   * Consider the following contingency table:
   * --------------
   * | y1  |   | n1
   * --------------
   * | y2  |   | n2
   * --------------
   * | m1  |   | 
   *
   * with Yi ~ Binom(ni, pi_i)
   *
   * rnchypgeom() samples Y1|Y1+Y2=m1
   */
  double rnchypgeom(const double& m1, const double& n1, 
                    const double& n2, const double& psi, 
                    const double& delta=1e-14); 

  /**** binomial distribution ****/
  
  /* Random Numbers */
  int rbinom (const int &, const double &);
  
  Matrix<double> rbinom ( const int &, const int &,
                          const int &, const double &);
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
  
  /* Random Numbers */
  double rchisq (const double &);

  Matrix<double> rchisq ( const int &, const int &,
                          const double &);

  /* CDFs */
  double pchisq(const double &, const double &);

  Matrix<double> pchisq(const int &, const int &,
                        const double &, const double &);
  /* PDFs */
  double dchisq(const double &, const double &);

  Matrix<double> dchisq(const int &, const int &,
                        const double &, const double &);
  
  /**** The Exponential Distribution ****/

  /* Random Numbers */
  double rexp (const double &);

  Matrix<double> rexp ( const int &rows, const int &cols,
                        const double &);

  /* CDFs */
  double pexp(const double &, const double &);

  Matrix<double> pexp(const int &, const int &, const double &,
                      const double &);

  /* PDFs */
  double dexp(const double &, const double &);

  Matrix<double> dexp(const int &, const int &, const double &,
                      const double &);

  /**** The f Distribution ****/

  /* Random Numbers */
  double rf(const double &, const double &);
  
  /* CDFs */
  double pf(const double &, const double &, const double &);

  Matrix<double> pf(const int &, const int &, const double &,
                    const double &, const double &);
  /* PDFs */
  double df(const double &, const double &, const double &);

  Matrix<double> df(const int &, const int &, const double &,
                    const double &, const double &);

  /**** The Gamma Distribution ****/

  /* Random Numbers */
  double rgamma (const double &, const double &);

  Matrix<double> rgamma ( const int &, const int &,
                          const double &, const double &);

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

  /* Others */
  double rgamma1 (const double &);
  
  /**** The Logistic Distribution ****/

  /* Random Numbers */
  double rlogis (const double &, const double &);

  Matrix<double> rlogis ( const int &, const int &,
                          const double &, const double &);
  
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

  /* Random Numbers */
  double rlnorm(const double &logmean = 0.0,
                const double &logsd = 1.0);

  Matrix<double> rlnorm(const int &, const int &,
                        const double &logmean = 0.0,
                        const double &logsd = 1.0);

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

  /* Random Numbers */
  double rlnorm(const double & , const double &);

  Matrix<double> rlnorm(const int &, const int &, const double &,
                        const double &);


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
  
  /* Random Numbers */
  double rnorm (const double &mu=0.0, const double &sigma=1.0);

  Matrix<double> rnorm (const int &rows, const int &cols,
                        const double &mu=0.0, const double &sigma=1.0);

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
  double rnorm1();
  
  /* Returns the univariate standard normal cumulative distribution
   * function (CDF)
   */
  double pnorm1(double, const double& eps=DBL_EPSILON,
                const double& maxit=500, const double& minit=10);

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
  
  /* Random Numbers */
  int rpois(const double &);

  Matrix<int> rpois (const int &, const int &, const double &);

  /* CDFs */

  double ppois(const double &, const double &);

  Matrix<double> ppois( const int &, const int &,
                        const double &, const double &);
  /* PDFs */
  double dpois(const int &, const double &);
  
  Matrix<double> dpois( const int &, const int &,
                        const double &, const double &);

  /**** The t Distribution ****/

  /* Random Numbers */
  double rt (const double &, const double &, const double &);

  Matrix<double> rt ( const int &, const int &,
                      const double &, const double &,
                      const double &);

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

  /* Random Numbers */
  inline double runif () {
    return rng();
  }
  
  Matrix<double> runif (const int &, const int &);

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
  /* Random Numbers */
  double rweibull (const double &, const double &scale = 1.0);

  Matrix<double> rweibull(const int &, const int &, const double &,
                          const double &scale = 1.0);
  
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
  
  /* Inverse Chi Squared */
  double richisq (const double &);

  Matrix<double> richisq (const int &, const int &,
                          const double &);


  /* Inverse Gamma */
  double rigamma (const double &, const double &);

  Matrix<double> rigamma( const int &, const int &,
                          const double &, const double &);

  /* Wishart (Only for Matrix) */ 
  Matrix<double> rwish(const int &, const Matrix<double> &);

  /* Dirichlet generator */
  Matrix<double> rdirich(const Matrix<double> &);
    
  /* Multivariate Normal generator */
  Matrix<double> rmvnorm(const Matrix<double> &,const Matrix<double> &);

  /* Multivariate t */
  Matrix<double> rmvt (const Matrix<double> &, const double &);
  
  /* Bernoulli */
  
  int rbern (const double &);

  Matrix<double> rbern (const int &, const int &,
                        const double &);

  /* Beta-Binomial */
  int rbetabin (const int &, const double &, const double &);

  Matrix<double> rbetabin ( const int &, const int &,
                            const int &, const double &,
                            const double &);

  /* Truncated Normal */
  
  /* Simulates from a truncated Normal distribution */
  double rtnorm(const double &, const double &, const double &, 
                const double &);

  /* Simulate from a truncated Normal distribution */
  double rtnorm_combo(const double &, const double &, 
              const double &, const double &);

  /* Sample from a truncated Normal distribution */
  double rtbnorm_slice( const double &, const double &,
                        const double &, const int &iter = 10);

  /* Sample from a truncated Normal distribution */
  double rtanorm_slice( const double &, const double &,
                        const double &, const int &iter = 10);

  /* Sample from a truncated Normal distribution */
  double rtbnorm_combo( const double &, const double &,
                        const double &, const int &iter=10);

  /* Sample from a truncated Normal distribution */
  double rtanorm_combo( const double &, const double &,
                        const double &, const int &iter=10);

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
