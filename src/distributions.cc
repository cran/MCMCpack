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
 * scythestat/distributions.cc
 *
 * Provides implementations of PDFs, CDFs, and some common functions
 * (gamma, beta, etc).
 *
 */

#ifndef SCYTHE_DISTRIBUTIONS_CC
#define SCYTHE_DISTRIBUTIONS_CC

#include <iostream>
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>

#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#ifdef SCYTHE_COMPILE_DIRECT
#include "distributions.h"
#include "error.h"
#include "util.h"
#include "ide.h"
#include "stat.h"
#include "la.h"
#else
#include "scythestat/distributions.h"
#include "scythestat/error.h"
#include "scythestat/util.h"
#include "scythestat/ide.h"
#include "scythestat/stat.h"
#include "scythestat/la.h"
#endif

#ifndef HAVE_TRUNC
inline double trunc(double x) throw () 
{
	if (x >= 0)
		return floor(x);
	else
		return ceil(x);
}
#endif

/* Many random number generators, pdfs, cdfs, and functions (gamma,
 * etc) in this file are based on code from the R Project, version
 * 1.6.0-1.7.1.  This code is available under the terms of the GNU
 * GPL.  Original copyright:
 * 
 * Copyright (C) 1998      Ross Ihaka
 * Copyright (C) 2000-2002 The R Development Core Team
 * Copyright (C) 2003      The R Foundation
 */

namespace SCYTHE {
  
  /* Non-convergence flag setup */
#ifdef __MINGW32__
  static bool THROW_ON_NONCONV = false;
#else
  namespace {
    bool THROW_ON_NONCONV = false;
  }
#endif

  void throw_on_nonconv (bool t) {
    THROW_ON_NONCONV = t;
  }

  /*************
   * Functions *
   *************/
  
  /* The gamma function */
  double 
  gammafn (const double &x)
  {
    const double gamcs[22] = {
      +.8571195590989331421920062399942e-2,
      +.4415381324841006757191315771652e-2,
      +.5685043681599363378632664588789e-1,
      -.4219835396418560501012500186624e-2,
      +.1326808181212460220584006796352e-2,
      -.1893024529798880432523947023886e-3,
      +.3606925327441245256578082217225e-4,
      -.6056761904460864218485548290365e-5,
      +.1055829546302283344731823509093e-5,
      -.1811967365542384048291855891166e-6,
      +.3117724964715322277790254593169e-7,
      -.5354219639019687140874081024347e-8,
      +.9193275519859588946887786825940e-9,
      -.1577941280288339761767423273953e-9,
      +.2707980622934954543266540433089e-10,
      -.4646818653825730144081661058933e-11,
      +.7973350192007419656460767175359e-12,
      -.1368078209830916025799499172309e-12,
      +.2347319486563800657233471771688e-13,
      -.4027432614949066932766570534699e-14,
      +.6910051747372100912138336975257e-15,
      -.1185584500221992907052387126192e-15,
    };


    double y = fabs(x);

    if (y <= 10) {

      /* Compute gamma(x) for -10 <= x <= 10
       * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
       * first of all. */

      int n = (int) x;
      if (x < 0)
        --n;
      
      y = x - n;/* n = floor(x)  ==>  y in [ 0, 1 ) */
      --n;
      double value = INTERNAL::chebyshev_eval(y * 2 - 1, gamcs, 22)
        + .9375;
      
      if (n == 0)
        return value;/* x = 1.dddd = 1+y */

      if (n < 0) {
        /* compute gamma(x) for -10 <= x < 1 */

        /* If the argument is exactly zero or a negative integer */
        /* then return NaN. */
        if (x == 0 || (x < 0 && x == n + 2))
          throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
              __LINE__, "x is 0 or a negative integer");

        /* The answer is less than half precision */
        /* because x too near a negative integer. */
        if (x < -0.5 && std::fabs(x - (int)(x - 0.5) / x) < 67108864.0)
          throw scythe_precision_error(__FILE__, __PRETTY_FUNCTION__,
              __LINE__,
              std::string("Answer < 1/2 precision because x is ")
              & "too near a negative integer");

        /* The argument is so close to 0 that the result
         * * would overflow. */
        if (y < 2.2474362225598545e-308)
          throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
              __LINE__, "x too close to 0");

        n = -n;

        for (int i = 0; i < n; i++)
          value /= (x + i);
        
        return value;
      } else {
        /* gamma(x) for 2 <= x <= 10 */

        for (int i = 1; i <= n; i++) {
          value *= (y + i);
        }
        return value;
      }
    } else {
      /* gamma(x) for   y = |x| > 10. */

      if (x > 171.61447887182298)    /* Overflow */
        throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Overflow");

      if (x < -170.5674972726612)       /* Underflow */
        throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Underflow");

      double value = std::exp((y - 0.5) * std::log(y) - y 
          + M_LN_SQRT_2PI + INTERNAL::lngammacor(y));

      if (x > 0)
        return value;

      if (std::fabs((x - (int)(x - 0.5))/x) < 67108864.0)
        throw scythe_precision_error(__FILE__, __PRETTY_FUNCTION__,
            __LINE__,
            std::string("Answer < 1/2 precision because x is ")
            & "too near a negative integer");

      double sinpiy = std::sin(M_PI * y);
      
      if (sinpiy == 0)     /* Negative integer arg - overflow */
        throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Overflow");

      return -M_PI / (y * sinpiy * value);
    }
  }

  /* The natural log of the absolute value of the gamma function */
  double
  lngammafn(const double &x)
  {
    if (x <= 0 && x == (int)x)
      throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
             __LINE__, "x is 0 or a negative integer");

    double y = fabs(x);

    if (y <= 10)
      return log(fabs(gammafn(x)));

    if (y > 2.5327372760800758e+305)
      throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
             "Overflow");

    if (x > 0) /* i.e. y = x > 10 */
      return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x
        + INTERNAL::lngammacor(x);
    
    /* else: x < -10; y = -x */
    double sinpiy = fabs(sin(M_PI * y));

    if (sinpiy == 0) /* Negative integer argument */
      throw scythe_exception("UNEXPECTED ERROR",
           __FILE__, __PRETTY_FUNCTION__, __LINE__,
           "ERROR:  Should never happen!");

    double ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy)
      - INTERNAL::lngammacor(y);

    if(fabs((x - (int)(x - 0.5)) * ans / x) < 1.490116119384765696e-8)
      throw scythe_precision_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, std::string("Answer < 1/2 precision because x is ")
           & "too near a negative integer");
    
    return ans;
  }

  /* The beta function */
  double
  betafn(const double &a, const double &b)
  {
    if (a <= 0 || b <= 0)
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
             "a or b < 0");

    if (a + b < 171.61447887182298) /* ~= 171.61 for IEEE */
      return gammafn(a) * gammafn(b) / gammafn(a+b);

    double val = lnbetafn(a, b);
    if (val < -708.39641853226412)
      throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
             "Underflow");
    
    return std::exp(val);
  }

  /* The natural log of the beta function */
  double
  lnbetafn (const double &a, const double &b)
  {
    double p, q;

    p = q = a;
    if(b < p) p = b;/* := min(a,b) */
    if(b > q) q = b;/* := max(a,b) */

    if (p <= 0 || q <= 0)
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
             "a or b <= 0");

    if (p >= 10) {
      /* p and q are big. */
      double corr = INTERNAL::lngammacor(p) + INTERNAL::lngammacor(q)
        - INTERNAL::lngammacor(p + q);
      return std::log(q) * -0.5 + M_LN_SQRT_2PI + corr
        + (p - 0.5) * log(p / (p + q)) + q * log(1 + (-p / (p + q)));
    } else if (q >= 10) {
      /* p is small, but q is big. */
      double corr = INTERNAL::lngammacor(q)
        - INTERNAL::lngammacor(p + q);
      return lngammafn(p) + corr + p - p * std::log(p + q)
        + (q - 0.5) * log(1 + (-p / (p + q)));
    }
    
    /* p and q are small: p <= q > 10. */
    return std::log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
  }

  /* Compute the factorial of a non-negative integer */
  int
  factorial (const int &n)
  {
    if (n < 0)
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
             "n < 0");

    if (n == 0)
      return 1;

    return n * factorial(n - 1);
  }

  /* Compute the natural log of the factorial of a non-negative
   * integer
   */
  double
  lnfactorial (const int &n)
  {
    if (n < 0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "n < 0");
      
    double x = n+1;
    double cof[6] = {
      76.18009172947146, -86.50532032941677,
      24.01409824083091, -1.231739572450155,
      0.1208650973866179e-2, -0.5395239384953e-5
    };
    double y = x;
    double tmp = x + 5.5 - (x + 0.5) * ::log(x + 5.5);
    double ser = 1.000000000190015;
    for (int j = 0; j <= 5; j++) {
      ser += (cof[j] / ++y);
    }
    return(std::log(2.5066282746310005 * ser / x) - tmp);
  }
  
  /*********************************
   * Fully Specified Distributions * 
   *********************************/
  
  /**** The Beta Distribution ****/

  /* CDFs */
  double
  pbeta(const double& x, const double& pin, const double& qin)
  {
    if (pin <= 0 || qin <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "pin or qin <= 0");
    }
    
    if (x <= 0)
      return 0.;
    if (x >= 1)
      return 1.;
    
    return INTERNAL::pbeta_raw(x,pin,qin);
  }

  Matrix<double>
  pbeta(const int& rows, const int& cols, const double& x,
  const double& pin, const double& qin)
  {
    int size = rows*cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);  
    for (int i=0; i<size; i++)
      temp[i] = pbeta(x,pin,qin);
    
    return temp;
  }
  
  /* PDFs */
  double
  dbeta(const double& x, const double& a, const double& b)
  {
    if ((x < 0.0) || (x > 1.0)) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "x not in [0,1]");
    }
    if (a < 0.0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "a < 0");
    }
    if (b < 0.0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "b < 0");
    }

    return (std::pow(x, (a-1.0)) * std::pow((1.0-x), (b-1.0)) )
      / betafn(a,b);
  }

  Matrix<double>
  dbeta(const int& rows, const int& cols, const double& x,
  const double& a, const double& b)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dbeta(x,a,b);
    
    return temp;  
  }
  
  /* Returns the natural log of the ordinate of the Beta density
   * evaluated at x with Shape1 a, and Shape2 b
   */
  double
  lndbeta1(const double& x, const double& a, const double& b)
  { 
    if ((x < 0.0) || (x > 1.0)) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "x not in [0,1]");
    }
    if (a < 0.0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "a < 0");
    }
    if (b < 0.0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "b < 0");
    }
      
    return (a-1.0) * std::log(x) + (b-1) * std::log(1.0-x)
      - lnbetafn(a,b);
  }

  /**** The Binomial Distribution ****/

  /* CDFs */
  double
  pbinom(const double &x, const double &n,const double &p) 
  {
    double N = std::floor(n + 0.5);
      
    if (N <= 0 || p < 0 || p > 1){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "floor(n + 0.5) <= 0 or p < 0 or p > 1");
    }
    double X = std::floor(x);
      
    if (X < 0.0)
      return 0;
    
    if (N <= X)
      return 1;
      
    return pbeta(1 - p, N - X, X + 1);
  }

  Matrix<double>
  pbinom (const int& rows, const int& cols, const double& x,
    const double& n, const double& p)
  {
    int size = rows*cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for(int i = 0; i < size; i++)
      temp[i] = pbinom(x,n,p);
      
    return temp;
  }    
  
  /* PDFs */
  double
  dbinom(const double& x, const double& n, const double& p)
  {
    if (p < 0 || p > 1) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "p not in [0,1]");
    }
    
    double N = floor(n + 0.5);
    double X = floor(x + 0.5);
    
    return INTERNAL::dbinom_raw(X, N, p, 1 - p);
  }

  Matrix<double>
  dbinom( const int& rows, const int& cols, const double& x,
    const double& n, const double& p)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dbinom(x,n,p);
    
    return temp;  
  }

  /**** The Chi Squared Distribution ****/
  
  /* CDFs */
  double
  pchisq(const double& x, const double& df)
  {
    return pgamma(x, df/2.0, 2.0);
  }


  Matrix<double>
  pchisq (const int& rows, const int& cols, const double& x,
    const double& df)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    
    for (int i = 0; i < size; i++)
      temp[i] = pchisq(x,df);
    
    return temp;
  }    

  /* PDFs */
  double
  dchisq(const double& x, const double& df)
  {
    return dgamma(x, df / 2.0, 2.0);
  }

  Matrix<double>
  dchisq( const int& rows, const int& cols, const double& x,
    const double& df)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dchisq(x,df);
      
    return temp;
  }

  /**** The Exponential Distribution ****/

  /* CDFs */
  double
  pexp(const double& x, const double& scale)
  {
    if (scale <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "scale <= 0");
    }
    
    if (x <= 0)
      return 0;
    
    return (1 - std::exp(-x*scale));
  }

  Matrix<double>
  pexp( const int& rows, const int& cols, const double& x,
  const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = pexp(x,scale);
    
    return temp;
  }
  
  /* PDFs */
  double
  dexp(const double& x, const double& scale)
  {
    if (scale <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "scale <= 0");
    }
    
    if (x < 0)
      return 0;
      
    return std::exp(-x * scale) * scale;
  }

  Matrix<double>
  dexp( const int& rows, const int& cols, const double& x,
  const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dexp(x,scale);
    
    return temp;
  }

  /**** The f Distribution ****/

  /* CDFs */
  double
  pf(const double& x, const double& n1, const double& n2)
  {
    if (n1 <= 0 || n2 <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "n1 or n2 <= 0");
    }
  
    if (x <= 0)
      return 0;
    
    if (n2 > 4e5)
      return pchisq(x*n1,n1);
    if (n1 > 4e5)
      return 1-pchisq(n2/x,n2);
    
    return (1-pbeta(n2/(n2+n1*x),n2/2.0,n1/2.0));
  }

  Matrix<double>
  pf (const int& rows, const int& cols, const double& x,
      const double& n1, const double& n2)
  {
    int size = rows * cols;
    if (size <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for(int i = 0; i < size; i++)
      temp [i] = pf(x,n1,n2);
    
    return temp;
  }
  
  /* PDFs */
  double
  df(const double& x, const double& m, const double& n)
  {
    double dens;
    
    if (m <= 0 || n <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "m or n <= 0");
    }
    
    if (x <= 0)
      return 0;
      
    double f = 1 / (n + x * m);
    double q = n * f;
    double p = x * m * f;
    
    if (m >= 2) {
      f = m * q / 2;
      dens = INTERNAL::dbinom_raw((m - 2) / 2,(m + n - 2) / 2, p, q);
    } else {
      f = (m * m * q) /(2 * p * (m + n));
      dens = INTERNAL::dbinom_raw(m / 2,(m + n)/ 2, p, q);
    }
    
    return f*dens;
  }

  Matrix<double>
  df( const int& rows, const int& cols, const double& x,
      const double& m, const double& n)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = df(x, m, n);
    
    return temp;
  }

  /**** The Gamma Distribution ****/

  /* CDFs */
  double
  pgamma (double x, const double &alph, const double &scale)
  {
    const double xbig = 1.0e+8, xlarge = 1.0e+37, 
      alphlimit = 1000.;/* normal approx. for alph > alphlimit */
      
    int lower_tail = 1;

    double pn1, pn2, pn3, pn4, pn5, pn6, arg, a, b, c, an, osum, sum;
    long n;
    int pearson;

    /* check that we have valid values for x and alph */

    if(alph <= 0. || scale <= 0.)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "alph or scale <= 0");

    x /= scale;
    
    if (x <= 0.)
      return 0.0;

    /* use a normal approximation if alph > alphlimit */

    if (alph > alphlimit) {
      pn1 = std::sqrt(alph) * 3. * (std::pow(x/alph, 1./3.) + 1.
            / (9. * alph) - 1.);
      return pnorm(pn1, 0., 1.);
    }

    /* if x is extremely large __compared to alph__ then return 1 */

    if (x > xbig * alph)
      return 1.0;

    if (x <= 1. || x < alph) {
      pearson = 1;/* use pearson's series expansion. */
      arg = alph * std::log(x) - x - lngammafn(alph + 1.);
      c = 1.;
      sum = 1.;
      a = alph;
      do {
        a += 1.;
        c *= x / a;
        sum += c;
      } while (c > DBL_EPSILON);
      arg += std::log(sum);
    }
    else { /* x >= max( 1, alph) */
      pearson = 0;/* use a continued fraction expansion */
      arg = alph * std::log(x) - x - lngammafn(alph);
      a = 1. - alph;
      b = a + x + 1.;
      pn1 = 1.;
      pn2 = x;
      pn3 = x + 1.;
      pn4 = x * b;
      sum = pn3 / pn4;
      for (n = 1; ; n++) {
        a += 1.;/* =   n+1 -alph */
        b += 2.;/* = 2(n+1)-alph+x */
        an = a * n;
        pn5 = b * pn3 - an * pn1;
        pn6 = b * pn4 - an * pn2;
        if (std::fabs(pn6) > 0.) {
          osum = sum;
          sum = pn5 / pn6;
          if (std::fabs(osum - sum) <= DBL_EPSILON * min(1., sum))
            break;
        }
        pn1 = pn3;
        pn2 = pn4;
        pn3 = pn5;
        pn4 = pn6;
        if (std::fabs(pn5) >= xlarge) {
          /* re-scale terms in continued fraction if they are large */
          pn1 /= xlarge;
          pn2 /= xlarge;
          pn3 /= xlarge;
          pn4 /= xlarge;
        }
      }
      arg += std::log(sum);
    }

    lower_tail = (lower_tail == pearson);

    sum = exp(arg);

    return (lower_tail) ? sum : 1 - sum;
  }
  
  Matrix<double>
  pgamma( const int& rows, const int& cols, const double& x,
    const double& alph, const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = pgamma(x, alph, scale);
    
    return temp;
  }
  

  /* PDFs */
  double
  dgamma(const double& x, const double& shape, const double& scale)
  {
    if (shape <= 0 || scale <= 0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "shape or scale <= 0");


    if (x < 0)
      return 0.0;
    
    if (x == 0) {
      if (shape < 1)
        throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "x == 0 and shape < 1");
      
      if (shape > 1)
        return 0.0;
      
      return 1 / scale;
    }
    
    if (shape < 1) { 
      double pr = INTERNAL::dpois_raw(shape, x/scale);
      return pr * shape / x;
    }
    
    /* else  shape >= 1 */
    double pr = INTERNAL::dpois_raw(shape - 1, x / scale);
    return pr / scale;
  }

  Matrix<double>
  dgamma( const int& rows, const int& cols, const double& x,
    const double& shape, const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dgamma(x, shape, scale);
    
    return temp;
  }
  
  /**** The Logistic Distribution ****/

  /* CDFs */
  double
  plogis (const double& x, const double& location,
    const double& scale)
  {
    if (scale <= 0.0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "scale <= 0");
    }
    
    double X = (x-location) / scale;
      
    X = std::exp(-X);
      
    return 1 / (1+X);
  }

  Matrix<double>
  plogis (const int& rows, const int& cols, const double& x,
    const double& location, const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = plogis(x,location,scale);
    
    return temp;
  }    
  
  /* PDFs */
  double
  dlogis( const double& x, const double& location,
    const double& scale)
  {
    if (scale <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "scale <= 0");
    }
    
    double X = (x - location) / scale;
    double e = std::exp(-X);
    double f = 1.0 + e;
      
    return e / (scale * f * f);
  }

  Matrix<double>
  dlogis( const int& rows, const int& cols, const double& x,
    const double& location, const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dlogis(x,location,scale);
    
    return temp;
  }

  double
  lndlogis( const double& x, const double& location,
    const double& scale)
  {
    if (scale <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "scale <= 0");
    }
    
    double X = (x - location) / scale;
    double e = std::exp(-X);
    double f = 1.0 + e;
    return std::log(e) - std::log(scale) - 2.0*std::log(f);
  }

  Matrix<double>
  lndlogis( const int& rows, const int& cols, const double& x,
    const double& location, const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = lndlogis(x,location,scale);
    
    return temp;
  }

  /**** The Log Normal Distribution ****/

  /* CDFs */

  double
  plnorm (const double& x, const double &logmean,
    const double & logsd)
  {
    if (logsd <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "logsd <= 0");
    }
    
    if (x > 0)
      return pnorm(std::log(x), logmean, logsd);
    
    return 0;
  }

  Matrix<double>
  plnorm (const int& rows, const int& cols, const double& x,
    const double& logmean, const double& logsd)
  {
    int size = rows * cols;
    if (size <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for(int i=0; i<size; i++)
      temp[i] = plnorm(x,logmean,logsd);
    
    return temp;
  }  

  /* PDFs */
  double
  dlnorm( const double& x, const double& logmean,
    const double& logsd)
  {
    if (logsd <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "logsd <= 0");
    }
    
    if (x == 0)
      return 0;
    
    double y = (::log(x) - logmean) / logsd;
    
    return (1 / (::sqrt(2 * M_PI))) * ::exp(-0.5 * y * y) / (x * logsd);
  }

  Matrix<double>
  dlnorm( const int& rows, const int& cols, const double& x,
    const double& logmean, const double& logsd)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dlnorm(x, logmean, logsd);
    
    return temp;
  }

  /**** The Negative Binomial Distribution ****/

  /* CDFs */
  double
  pnbinom(const double& x, const double& n, const double& p)
  {
    if (n <= 0 || p <= 0 || p >= 1){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "n <= 0 or p not in (0,1)");
    }
    
    double X = floor(x + 1e-7);
      
    if (X < 0)
      return 0.0;
    
    return pbeta(p, n, X + 1);
  }

  Matrix<double>
  pnbinom(const int& rows, const int& cols, const double& x,
    const double& n, const double& p)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    
    for (int i = 0; i < size; i++)
      temp[i] = pnbinom(x,n,p);
    
    return temp;
  }    

  /* PDFs */
  double
  dnbinom(const double& x, const double& n, const double& p)
  {
    if (p < 0 || p > 1 || n <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "p not in [0,1] or n <= 0");
    }
    
    if (x < 0)
      return 0;
    
    double X = floor(x + 0.5);
    
    double prob = INTERNAL::dbinom_raw(n, X + n, p, 1 - p);
    double P = (double) n / (n + x);
    
    return P * prob;
  }

  Matrix<double>
  dnbinom(const int& rows, const int& cols, const double& x,
    const double& n, const double& p)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dnbinom(x,n,p);
    
    return temp;
  }

  
  /**** The Normal Distribution ****/

  /* CDFs */
  double
  pnorm (const double &x, const double &mu, const double &sigma)
  
  {
    if (sigma <= 0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "negative standard deviation");

    return pnorm2((x - mu) / sigma, true, false);
  }
  
  Matrix<double>
  pnorm(const int& rows, const int& cols, const double& x,
  const double& mu, const double& sigma)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = pnorm(x,mu,sigma);
    
    return temp;
  }
  
  /* PDFs */
  double
  dnorm(const double& x, const double& mu,
  const double& sigma)
  {
    if (sigma <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "negative standard deviation");
    }
    
    double X = (x - mu) / sigma;
    
    return (M_1_SQRT_2PI * std::exp(-0.5 * X * X) / sigma);
  }

  Matrix<double>
  dnorm(const int& rows, const int& cols, const double& x,
  const double& mu, const double& sigma)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dnorm(x,mu,sigma);
    
    return temp;
  }
  
  /* A newer version of pnorm for 0.4+.  The pnorm wrapper has been
   * updated to use this function, as have library calls that
   * previously used pnorm1.
   */

	 // Many original comments left in for reference.

#define SIXTEN 16
#define do_del(X)              \
  xsq = trunc(X * SIXTEN) / SIXTEN;        \
  del = (X - xsq) * (X + xsq);          \
  if(log_p) {              \
      *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);  \
      if((lower && x > 0.) || (upper && x <= 0.))      \
      *ccum = log1p(-exp(-xsq * xsq * 0.5) *     \
        exp(-del * 0.5) * temp);    \
  }                \
  else {                \
      *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;  \
      *ccum = 1.0 - *cum;            \
  }

#define swap_tail            \
  if (x > 0.) {/* swap  ccum <--> cum */      \
      temp = *cum; if(lower) *cum = *ccum; *ccum = temp;  \
  }

  void
  pnorm_both(  double x, double *cum, double *ccum, int i_tail,
              bool log_p)
  {
    const double a[5] = {
      2.2352520354606839287,
      161.02823106855587881,
      1067.6894854603709582,
      18154.981253343561249,
      0.065682337918207449113
    };
    const double b[4] = {
      47.20258190468824187,
      976.09855173777669322,
      10260.932208618978205,
      45507.789335026729956
    };
    const double c[9] = {
      0.39894151208813466764,
      8.8831497943883759412,
      93.506656132177855979,
      597.27027639480026226,
      2494.5375852903726711,
      6848.1904505362823326,
      11602.651437647350124,
      9842.7148383839780218,
      1.0765576773720192317e-8
    };
    const double d[8] = {
      22.266688044328115691,
      235.38790178262499861,
      1519.377599407554805,
      6485.558298266760755,
      18615.571640885098091,
      34900.952721145977266,
      38912.003286093271411,
      19685.429676859990727
    };
    const double p[6] = {
      0.21589853405795699,
      0.1274011611602473639,
      0.022235277870649807,
      0.001421619193227893466,
      2.9112874951168792e-5,
      0.02307344176494017303
    };
    const double q[5] = {
      1.28426009614491121,
      0.468238212480865118,
      0.0659881378689285515,
      0.00378239633202758244,
      7.29751555083966205e-5
    };
    
    double xden, xnum, temp, del, eps, xsq, y;
    int i, lower, upper;

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = std::fabs(x);
    if (y <= 0.67448975) {
      /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
      if (y > eps) {
        xsq = x * x;
        xnum = a[4] * xsq;
        xden = xsq;
        for (i = 0; i < 3; ++i) {
          xnum = (xnum + a[i]) * xsq;
          xden = (xden + b[i]) * xsq;
        }
      } else xnum = xden = 0.0;
      
      temp = x * (xnum + a[3]) / (xden + b[3]);
      if(lower)  *cum = 0.5 + temp;
      if(upper) *ccum = 0.5 - temp;
      if(log_p) {
        if(lower)  *cum = log(*cum);
        if(upper) *ccum = log(*ccum);
      }
    } else if (y <= M_SQRT_32) {
      /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) 
       * ~= 5.657 */

      xnum = c[8] * y;
      xden = y;
      for (i = 0; i < 7; ++i) {
        xnum = (xnum + c[i]) * y;
        xden = (xden + d[i]) * y;
      }
      temp = (xnum + c[7]) / (xden + d[7]);
      do_del(y);
      swap_tail;
    } else if (log_p
              || (lower && -37.5193 < x && x < 8.2924)
              || (upper && -8.2929 < x && x < 37.5193)
        ) {
      /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
      xsq = 1.0 / (x * x);
      xnum = p[5] * xsq;
      xden = xsq;
      for (i = 0; i < 4; ++i) {
        xnum = (xnum + p[i]) * xsq;
        xden = (xden + q[i]) * xsq;
      }
      temp = xsq * (xnum + p[4]) / (xden + q[4]);
      temp = (M_1_SQRT_2PI - temp) / y;
      do_del(x);
      swap_tail;
    } else {
      if (x > 0) {
        *cum = 1.;
        *ccum = 0.;
      } else {
        *cum = 0.;
        *ccum = 1.;
      }
      if (THROW_ON_NONCONV)
        throw scythe_convergence_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, std::string("x (") & x & "did not converge");
    }

    return;
  }
#undef SIXTEN
#undef do_del
#undef swap_tail

  double
  pnorm2 (const double &x, const bool &lower_tail, const bool &log_p)
  {
    if (! finite(x))
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Quantile x is inifinte (+/-Inf) or NaN");

    double p, cp;
    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

    return (lower_tail ? p : cp);
  }
  
  /* Returns the quantile of the standard normal distribution 
   * associated with a given probability p
   */
  double
  qnorm1 (const double& in_p)
  {
    double lim = 10e-20;
    double p0 = -0.322232431088;
    double q0 = 0.0993484626060;
    double p1 = -1.0;
    double q1 = 0.588581570495;
    double p2 = -0.342242088547;
    double q2 = 0.531103462366;
    double p3 = -0.0204231210245;
    double q3 = 0.103537752850;
    double p4 = -0.453642210148e-4;
    double q4 = 0.38560700634e-2;
    double xp = 0.0;
    double p = in_p;
      
    if (p > 0.5)
      p = 1 - p;
        
    if (p < lim) {
      throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
             __LINE__, "p outside accuracy limit");
    }
      
    if (p == 0.5)
      return xp;
      
    double y = ::sqrt (::log (1.0 / ::pow (p, 2)));
    xp = y + ((((y * p4 + p3) * y + p2) * y + p1) * y + p0) /
      ((((y * q4 + q3) * y + q2) * y + q1) * y + q0);
      
    if (in_p < 0.5)
      xp = -1 * xp;
    
    return xp;
  }

  /* Return the natrual log of the normal PDF */
  double
  lndnorm (const double& x, const double& mu, const double& sigma)
  {
    if (sigma < 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "negative standard deviation");
    }
    
    if (sigma == 0.0){
      if (x != mu){
	return -std::numeric_limits<double>::infinity();
      }
      else {
	return std::numeric_limits<double>::infinity();
      }
    }
    
    double X = (x - mu) / sigma;
    
    return -(M_LN_SQRT_2PI  +  0.5 * X * X + std::log(sigma));
  }
  
  /**** The Poison Distribution ****/

  /* CDFs */
  double
  ppois(const double& x, const double& lambda)
  {
    if(lambda<=0.0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "lambda <= 0");
    }
    
    double X = floor(x + 1e-7);
    
    if (X < 0)
      return 0;
    if (lambda == 1)
      return 1;
    
    return 1 - pgamma(lambda, X + 1, 1.0);
  }

  Matrix<double>
  ppois(const int& rows, const int& cols, const double& x,
  const double& lambda)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = ppois(x,lambda);
      
    return temp;
  }
  
  /* PDFs */
  double
  dpois(const int &x, const double &lambda)
  {
    if (x < 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "x < 0");
    }
    if (lambda <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "lambda <= 0");
    }
    
    // compute log(x!)
    double xx = x+1;
    double cof[6] = {
      76.18009172947146, -86.50532032941677,
      24.01409824083091, -1.231739572450155,
      0.1208650973866179e-2, -0.5395239384953e-5
    };
    double y = xx;
    double tmp = xx + 5.5 - (xx + 0.5) * ::log(xx + 5.5);
    double ser = 1.000000000190015;
    for (int j = 0; j <= 5; j++) {
      ser += (cof[j] / ++y);
    }
    double lnfactx = ::log(2.5066282746310005 * ser / xx) - tmp;
      
    return (::exp( -1*lnfactx + x * ::log(lambda) - lambda));
  }

  Matrix<double>
  dpois(const int& rows, const int& cols, const double& x,
  const double& lambda)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dpois((int)x,lambda);
    
    return temp;
  }

  /**** The t Distribution ****/

  /* CDFs */
  double
  pt(const double& x, const double& n)
  {
    double val;
    
    if (n <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "n <= 0");
    }
    
    if (n > 4e5) {
      val = 1/(4*n);
      return pnorm2(x * (1 - val) / ::sqrt(1 + x * x * 2. * val), 
          true, false);
    }
    
    val = pbeta(n / (n + x * x), n / 2.0, 0.5);
    
    val /= 2;
    
    if (x <= 0)
      return val;
    else
      return 1 - val;
  }

  Matrix<double>
  pt(const int& rows, const int& cols, const double& x, const double& n)
  { 
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = pt(x,n);
    
    return temp;
  }
  
  /* PDFs */
  double
  dt(const double& x, const double& n)
  {
    double u;
    if (n <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "n <= 0");
    }
    
    double t = -INTERNAL::bd0(n/2., (n + 1) / 2.)
      + INTERNAL::stirlerr((n + 1) / 2.)
      - INTERNAL::stirlerr(n / 2.);
    if(x*x > 0.2*n)
      u = std::log(1+x*x/n)*n/2;
    else
      u = -INTERNAL::bd0(n/2., (n+x*x)/2.) + x*x/2;
    
    return std::exp(t-u)/std::sqrt(2*M_PI*(1+x*x/n));
  }

  Matrix<double>
  dt(const int& rows, const int& cols, const double& x, const double& n)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dt(x,n);
    
    return temp;
  }
  
  /* Other */
  
  /* Returns the univariate Student-t density evaluated at x 
   * with mean mu, scale sigma^2, and nu degrees of freedom
   */
  double
  dt1(const double& x, const double& mu, const double& sigma2, 
      const double& nu)
  {
    double logdens =   lngammafn((nu + 1.0) /2.0)
      - std::log(std::sqrt(nu * M_PI))
      - lngammafn(nu / 2.0) - std::log(std::sqrt(sigma2))
      - (nu + 1.0) / 2.0 * std::log(1.0 + (std::pow((x - mu), 2.0))
            / (nu * sigma2));
    
    return(std::exp(logdens));
  }

  /* Returns the natural log of the univariate Student-t density 
   * evaluated at x with mean mu, scale sigma^2, and nu 
   * degrees of freedom
   */
  double 
  lndt1(  const double& x, const double& mu, const double& sigma2, 
    const double& nu)
  {
    double logdens = lngammafn((nu+1.0)/2.0)
      - std::log(std::sqrt(nu*M_PI))
      - lngammafn(nu/2.0) - std::log(std::sqrt(sigma2))
      - (nu+1.0)/2.0 * std::log(1.0 + (std::pow((x-mu),2.0))
        /(nu * sigma2));
    
    return(logdens);
  }

  /**** The Uniform Distribution ****/

  /* CDFs */
  double
  punif(const double& x, const double& a, const double& b)
  {
    if (b <= a) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "b <= a");
    }
      
    if (x <= a)
      return 0.0;
        
    if (x >= b)
      return 1.0;
      
    return (x - a) / (b - a);
  }

  Matrix<double>
  punif(const int& rows, const int& cols, const double& x,
  const double& a, const double& b)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = punif(x,a,b);
    
    return temp;
  }

  /* PDFs */
  double
  dunif(const double& x, const double& a, const double& b)
  {
    if (b <= a) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "b <= a");
    }
    
    if (a <= x && x <= b)
      return 1.0 / (b - a);
    
    return 0.0;
  }

  Matrix<double>
  dunif(const int& rows, const int& cols, const double& x,
  const double& a, const double& b)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dunif(x,a,b);
    
    return temp;
  }

  /**** The Weibull Distribution ****/

  /* CDFs */
  double
  pweibull( const double& x, const double& shape,
      const double& scale)
  {
    if (shape <= 0 || scale <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "shape or scale <= 0");
    }
    
    if (x <= 0)
      return 0.0;
    
    return 1 - std::exp(-std::pow(x / scale, shape));
  }

  Matrix<double>
  pweibull( const int& rows, const int& cols, const double& x,
      const double& shape, const double& scale)
  {
    int size = rows * cols;
    if (size <= 0){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = pweibull(x,shape,scale);
    
    return temp;
  }

  /* PDFs */
  double
  dweibull( const double& x, const double& shape,
      const double& scale)
  {
    if (shape <= 0 || scale <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "shape or scale <= 0");
    }
    if (x < 0)
      return 0;
      
    double tmp1 = std::pow(x / scale, shape - 1);
    double tmp2 = tmp1*(x / scale);
      
    return shape * tmp1 * std::exp(-tmp2) / scale;
  }

  Matrix<double>
  dweibull( const int& rows, const int& cols, const double& x,
      const double& shape, const double& scale)
  {
    int size = rows * cols;
    if (size <= 0) {
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "Tried to create matrix of size <= 0");
    }
    Matrix<double> temp(rows, cols, false);
    for (int i = 0; i < size; i++)
      temp[i] = dweibull(x,shape,scale);
    
    return temp;
  }

  /************************************
   * Partially Finished Distributions *
   ************************************/
  
  /* Multivariate Normal */
  double
  lndmvn (const Matrix<double> &x, const Matrix<double> &mu, 
    const Matrix<double> &Sigma)
  {
    if (! x.isColVector()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "x not column vector");
    }
    if (! mu.isColVector()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "mu not column vector");
    }
    if (! Sigma.isSquare()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "Sigma not square");
    }
    if (mu.rows() != Sigma.rows() || x.rows() != Sigma.rows()){
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "mu, x have different number of rows than Sigma");
    }
    int k = mu.rows();
    return ( (-k/2.0)*::log(2*M_PI) -0.5 * ::log(~Sigma) 
       -0.5 * (!(x - mu)) * invpd(Sigma) * (x-mu) )[0];
  }

  /********************
   * Helper Functions *
   ********************/
  namespace INTERNAL {

    /* Evaluate a Chebysheve series at a given point */
    double
    chebyshev_eval (const double &x, const double *a,
        const int &n)
    {
      if (n < 1 || n > 1000)
        throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "n not on [1, 1000]");
  
      if (x < -1.1 || x > 1.1)
        throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "x not on [-1.1, 1.1]");
      
      double b0, b1, b2;
      b0 = b1 = b2 = 0;
  
      double twox = x * 2;
  
      for (int i = 1; i <= n; ++i) {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
      }
  
      return (b0 - b2) * 0.5;
    }

    /* Computes the log gamma correction factor for x >= 10 */
    double
    lngammacor(const double &x)
    {
      const double algmcs[15] = {
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
      };
    
      if (x < 10) {
        throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "This function requires x >= 10");  
      } else if (x >= 3.745194030963158e306) {
        throw scythe_range_error(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Underflow");
      
      } else if (x < 94906265.62425156) {
        double tmp = 10 / x;
        return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, 5) / x;
      }
      
      return 1 / (x * 12);
    }

    /* Helper for dpois and dgamma */
    double
    dpois_raw (const double &x, const double &lambda)
    {
      if (lambda == 0)
        return ( (x == 0) ? 1.0 : 0.0);

      if (x == 0)
        return std::exp(-lambda);

      if (x < 0)
        return 0.0;

      return std::exp(-stirlerr(x) - bd0(x, lambda))
        / std::sqrt(2 * M_PI * x);
    }

    /* Evaluates the "deviance part" */
    double
    bd0(const double &x, const double &np)
    {
      
      if(std::fabs(x - np) < 0.1 * (x + np)) {
        double v = (x - np) / (x + np);
        double s = (x - np) * v;
        double ej = 2 * x * v;
        v = v * v;
        for (int j = 1; ; j++) {
          ej *= v;
          double s1 = s + ej / ((j << 1) + 1);
          if (s1 == s)
            return s1;
          s = s1;
        }
      }
      
      return x * std::log(x / np) + np - x;
    }
  
    /* Computes the log of the error term in Stirling's formula */
    double
    stirlerr(const double &n)
    {
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
      
      /* error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0 */
      const double sferr_halves[31] = {
        0.0, /* n=0 - wrong, place holder only */
        0.1534264097200273452913848,  /* 0.5 */
        0.0810614667953272582196702,  /* 1.0 */
        0.0548141210519176538961390,  /* 1.5 */
        0.0413406959554092940938221,  /* 2.0 */
        0.03316287351993628748511048, /* 2.5 */
        0.02767792568499833914878929, /* 3.0 */
        0.02374616365629749597132920, /* 3.5 */
        0.02079067210376509311152277, /* 4.0 */
        0.01848845053267318523077934, /* 4.5 */
        0.01664469118982119216319487, /* 5.0 */
        0.01513497322191737887351255, /* 5.5 */
        0.01387612882307074799874573, /* 6.0 */
        0.01281046524292022692424986, /* 6.5 */
        0.01189670994589177009505572, /* 7.0 */
        0.01110455975820691732662991, /* 7.5 */
        0.010411265261972096497478567, /* 8.0 */
        0.009799416126158803298389475, /* 8.5 */
        0.009255462182712732917728637, /* 9.0 */
        0.008768700134139385462952823, /* 9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690  /* 15.0 */
      };
      double nn;
      
      if (n <= 15.0) {
        nn = n + n;
        if (nn == (int)nn)
          return(sferr_halves[(int)nn]);
        return (lngammafn(n + 1.) - (n + 0.5) * std::log(n) + n -
            std::log(std::sqrt(2 * M_PI)));
      }
      
      nn = n*n;
      if (n > 500)
        return((S0 - S1 / nn) / n);
      if (n > 80)
        return((S0 - (S1 - S2 / nn) / nn) / n);
      if (n > 35)
        return((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
      /* 15 < n <= 35 : */
      return((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
    }
  
    /* helper for pbeta */
    double
    pbeta_raw(const double& x, const double& pin, const double& qin)
    {
      double ans, c, finsum, p, ps, p1, q, term, xb, xi, y;
      int n, i, ib, swap_tail;
      
      const double eps = .5 * DBL_EPSILON;
      const double sml = DBL_MIN;
      const double lneps = std::log(eps);
      const double lnsml = std::log(eps);
      
      if (pin / (pin + qin) < x) {
        swap_tail = 1;
        y = 1 - x;
        p = qin;
        q = pin;
      } else {
        swap_tail=0;
        y = x;
        p = pin;
        q = qin;
      }
      
      if ((p + q) * y / (p + 1) < eps) {
        ans = 0;
        xb = p * std::log(max(y,sml)) - std::log(p) - lnbetafn(p,q);
        if (xb > lnsml && y != 0)
          ans = std::exp(xb);
        if (swap_tail)
          ans = 1-ans;
      } else {
        ps = q - std::floor(q);
        if (ps == 0)
          ps = 1;
        xb = p * std::log(y) - lnbetafn(ps, p) - std::log(p);
        ans = 0;
        if (xb >= lnsml) {
          ans = std::exp(xb);
          term = ans * p;
          if (ps != 1) {
            n = (int)max(lneps/std::log(y), 4.0);
            for(i = 1; i <= n; i++){
              xi = i;
              term *= (xi-ps)*y/xi;
              ans += term/(p+xi);
            }
          }
        }
        if (q > 1) {
          xb = p * std::log(y) + q * std::log(1 - y)
            - lnbetafn(p, q) - std::log(q);
          ib = (int) max(xb / lnsml, 0.0);
          term = std::exp(xb - ib * lnsml);
          c = 1 / (1 - y);
          p1 = q * c / (p + q - 1);
              
          finsum = 0;
          n = (int) q;
          if(q == n)
            n--;
          for (i = 1; i <= n; i++) {
            if(p1 <= 1 && term / eps <= finsum)
              break;
            xi = i;
            term = (q -xi + 1) * c * term / (p + q - xi);
            if (term > 1) {
              ib--;
              term *= sml;
            }
            if (ib == 0)
              finsum += term;
          }
          ans += finsum;
        }
        
        if(swap_tail)
          ans = 1-ans;
        ans = max(min(ans,1.),0.);
      }
      return ans;
    }
  

    double
    dbinom_raw (const double &x, const double &n, const double &p,
    const double &q)
    { 
      double f, lc;

      if (p == 0)
        return((x == 0) ? 1.0 : 0.0);
      if (q == 0)
        return((x == n) ? 1.0 : 0.0);

      if (x == 0) { 
        if(n == 0)
          return 1.0;
        
        lc = (p < 0.1) ? -bd0(n, n * q) - n * p : n * std::log(q);
        return(std::exp(lc));
      }
      if (x == n) { 
        lc = (q < 0.1) ? -bd0(n,n * p) - n * q : n * std::log(p);
        return(std::exp(lc));
      }

      if (x < 0 || x > n)
        return 0.0;

      lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) -
        bd0(n - x, n * q);
      
      f = (M_2PI * x * (n-x)) / n;

      return (std::exp(lc) / std::sqrt(f));
    }

  
  }  // end namespace INTERNAL
} // end namespace SCYTHE


#endif
