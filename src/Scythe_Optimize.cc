/* Scythe_Optimize.cc
 *
 * This file provides optimization functions for the Scythe
 * Statistical Library.
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

#ifndef SCYTHE_OPTIMIZE_CC
#define SCYTHE_OPTIMIZE_CC

#include <cmath>
#include <iostream>
#include "Scythe_Error.h"
#include "Scythe_Util.h"
#include "Scythe_Simulate.h"
#include "Scythe_Optimize.h"
#include "Scythe_LA.h"

// Avoid NameSpace Pollution
namespace SCYTHE {

  /* Functions (private to this file) that do very little... */
  namespace {
    template <class T>
    T
    donothing(const Matrix<T> &x) {return 0.0;}

    template <class T>
    T
    donothing(const T &x)  { return 0.0; }
  }


  /* Return the machine epsilon 
   * Notes: Algorithm taken from Sedgewick, Robert. 1992. Algorithms
   * in C++. Addison Wesley. pg. 561
   */
  template <class T>
  T
  epsilon()
  {
    T eps, del, neweps;
    del    = (T) 0.5;
    eps    = (T) 0.0;
    neweps = (T) 1.0;
  
    while ( del > 0 ) {
      if ( 1 + neweps > 1 ) {  /* Then the value might be too large */
        eps = neweps;    /* ...save the current value... */
        neweps -= del;    /* ...and decrement a bit */
      } else {      /* Then the value is too small */
        neweps += del;    /* ...so increment it */
      }
      del *= 0.5;      /* Reduce the adjustment by half */
    }

    return eps;
  }

  /* Calculate the definite integral of a function from a to b */
  template <class T>
  T
  intsimp(T (*fun)(const T &), const T &a, const T &b, const int &N)
  {
    if (a > b)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "Lower limit larger than upper");
    
    if (N <= 0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "Number of subintervals negative");

    T I = (T) 0;
    T w = (b - a) / N;
    for (int i = 1; i <= N; i++)
      I += w * (fun(a +(i - 1) *w) + 4 * fun(a - w / 2 + i * w) +
          fun(a + i * w)) / 6;
   
    return I;
  }
  
  /* Calculate the definite integral of a function from a to b
   * Notes: Algorithm taken from Sedgewick, Robert. 1992. Algorithms
   * in C++. Addison Wesley. pg. 562
   */
  template <class T>
  T
  adaptsimp(T (*fun)(const T &), const T &a, const T &b, const int &N,
      const T tol)
  {
    if (a > b)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "Lower limit larger than upper");
    
    if (N <= 0)
      throw scythe_invalid_arg (__FILE__, __PRETTY_FUNCTION__, __LINE__,
        "Number of subintervals negative");

    T I = intsimp(fun, a, b, N);
    if (::fabs(I - intsimp(fun, a, b, N / 2)) > tol)
      return adaptsimp(fun, a, (a + b) / 2, N, tol)
        + adaptsimp(fun, (a + b) / 2, b, N, tol);

    return I;
  }


  /* Numerically calculates the first derivative of a function
   * Notes: Algorithm taken from Nocedal and Wright. 1999. section 7.1
   * with some additional tricks from Press et al. NRC
   */
  template <class T>
  Matrix<T>
  gradfdif (T (*fun)(const Matrix<T> &, const Matrix<T> &,
         const Matrix<T> &), const Matrix<T> &theta,
      const Matrix<T> &y, const Matrix<T> &X)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");

    int k = theta.size();
    // stepsize CAREFUL-- THIS IS MACHINE-SPECIFIC!!!
    T h = std::sqrt(epsilon<T>()); // 2.25e-16
    //T h = std::sqrt(2.25e-16);

   
    Matrix<T> grad(k,1);
  
    for (int i = 0; i < k; ++i){
      Matrix<T> e(k,1);
      e[i] = h;
      Matrix<T> temp = theta + e;
      donothing(temp);
      e = temp - theta;
      grad[i] = (fun(theta + e, y, X) - fun(theta, y, X)) / e[i];
    }

    return grad;
  }


  /* Numerically calculates the first derivative of a function
   * Notes: Algorithm taken from Nocedal and Wright. 1999. section 7.1
   * with some additional tricks from Press et al. NRC
   */
  template <class T>
  T
  gradfdifls (T (*fun)(const Matrix<T> &, const Matrix<T> &,
           const Matrix<T> &),  const T &alpha,
        const Matrix<T> &theta, const Matrix<T> &p, 
        const Matrix<T> &y, const Matrix<T> &X)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
    if (! p.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "p not column vector");

    int k = theta.size();
    // stepsize CAREFUL-- THIS IS MACHINE-SPECIFIC!!!
    T h = std::sqrt(epsilon<T>()); //2.2e-16 
    //T h = std::sqrt(2.2e-16);

    T deriv;

    for (int i = 0; i < k; ++i) {
      T temp = alpha + h;
      donothing(temp);
      T e = temp - alpha;
      deriv = (fun(theta + (alpha + e) * p, y, X)
         - fun(theta + alpha * p, y, X)) / e;
    }
    
    return deriv;
  }



  /* Numerically calculates the gradient of a function
   * Notes: Algorithm taken from Nocedal and Wright. 1999. section 7.1
   * with some additional tricks from Press et al. NRC
   */
  template <class T>
  Matrix<T>
  jacfdif(Matrix<T> (*fun)(const Matrix<T> &), const Matrix<T> &theta)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
   
    Matrix<T> fval = fun(theta);

    int k = theta.rows();
    int n = fval.rows();
    // stepsize CAREFUL -- THIS IS MACHINE-SPECIFIC!!!
    T h = std::sqrt(epsilon<T>()); //2.2e-16
    //T h = std::sqrt(2.2e-16);
    Matrix<T> J(n,k);
    
    for (int i = 0; i < k; ++i) {
      Matrix<T> e(k,1);
      e[i] = h;
      Matrix<T> temp = theta + e;
      donothing(temp);
      e = temp - theta;
      Matrix<T> fthetae = fun(theta + e);
      Matrix<T> ftheta = fun(theta);
      for (int j = 0; j < n; ++j) {
        J(j,i) = (fthetae[j] - ftheta[j]) / e[i];
      }
    }
   
    return J;
  }

  /* Numerically calculates the gradient of a function
   * Notes: Algorithm taken from Nocedal and Wright. 1999. section 7.1
   * with some additional tricks from Press et al. NRC
   */
  template <class T>
  Matrix<T>
  jacfdif(Matrix<T> (*fun)(const Matrix<T> &, const Matrix<T> &), 
    const Matrix<T> &theta, const Matrix<T> &psi)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
   
    Matrix<T> fval = fun(theta, psi);
    
    int k = theta.rows();
    int n = fval.rows();
    // stepsize CAREFUL -- THIS IS MACHINE-SPECIFIC!!!
    //T h = std::sqrt(epsilon<T>()); //2.2e-16
    T h = std::sqrt(2.2e-16);
    Matrix<T> J(n, k);

    for (int i = 0; i < k; ++i) {
      Matrix<T> e(k,1);
      e[i] = h;
      Matrix<T> temp = theta + e;
      donothing(temp);
      e = temp - theta;
      Matrix<T> fthetae = fun(theta + e, psi);
      Matrix<T> ftheta = fun(theta, psi);
      for (int j = 0; j < n; ++j) {
        J(j,i) = (fthetae[j] - ftheta[j]) / e[i];
      }
    }
   
    return J;
  }

  /* Numerically calculates the Hessian of a function */
  template <class T>
  Matrix<T>
  hesscdif (T (*fun)(const Matrix<T> &, const Matrix<T> &,
         const Matrix<T> &), const Matrix<T> &theta,
      const Matrix<T> &y, const Matrix<T> &X)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
    
    T fval = fun(theta,y,X);

    int k = theta.rows();

    // stepsize CAREFUL -- THIS IS MACHINE-SPECIFIC!!!
    T h2 = (T) 1e-10;
    T h = ::sqrt(h2); 
    Matrix<T> H(k,k);

    for (int i=0; i<k; ++i) {
      Matrix<T> ei(k, 1);
      ei[i] = h;
      Matrix<T> temp = theta + ei;
      donothing(temp);
      ei = temp - theta;
      for (int j = 0; j < k; ++j){
        Matrix<T> ej(k,1);
        ej[j] = h;
        temp = theta + ej;
        donothing(temp);
        ej = temp - theta;
        
        if (i==j){
          H(i,i) = ( -fun(theta + 2.0 * ei, y, X) + 16.0 *
              fun(theta+ei,y,X) - 30.0 * fval + 16.0 *
              fun(theta-ei,y,X) -
              fun(theta-2.0 * ei, y, X)) / (12.0 * h2);
        } else {
          H(i,j) = ( fun(theta + ei + ej, y, X) - fun(theta+ei-ej, y, X)
              - fun(theta - ei + ej, y, X) + fun(theta-ei-ej, y, X))
            / (4.0 * h2);
        }
      }
    }
       
    return H;
  }

  /* Performs a linesearch to find the step length (\a alpha)
   * Notes: Algorithm taken from Nocedal and Wright. 1999. Procedure
   * 3.1
   */
  template <class T>
  T
  linesearch1(T (*fun)(const Matrix<T> &, const Matrix<T> &,
           const Matrix<T> &), const Matrix<T> &theta,
        const Matrix<T> &p, const Matrix<T> &y,
        const Matrix<T> &X)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
    if (! p.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "p not column vector");

    T alpha_bar = (T) 1.0;
    T rho = (T) 0.9;
    T c   = (T) 0.5;
    T alpha = alpha_bar;
    Matrix<T> fgrad = gradfdif(fun, theta, y, X);

    while (fun(theta + alpha * p, y, X) > (fun(theta, y, X) + c
             * alpha * t(fgrad) * p)[0]) {
      alpha = rho * alpha;
    }

    return alpha;
  }

  /* Performs a linesearch to find the step length (\a alpha)
   * Notes: Algorithm taken from Nocedal and Wright. 1999. Algorithm
   * 3.2
   */
  template <class T>
  T
  linesearch2(T (*fun)(const Matrix<T> &, const Matrix<T> &,
           const Matrix<T> &), const Matrix<T> &theta,
        const Matrix<T> &p, const Matrix<T> &y,
        const Matrix<T> &X)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
    if (! p.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "p not column vector");

    T alpha_last = (T) 0.0;
    T alpha_cur = (T) 1.0;
    T alpha_max = (T) 10.0;
    T c1 = (T) 1e-4;
    T c2 = (T) 0.5;
    int max_iter = 50;
    T fgradalpha0 = gradfdifls(fun, (T) 0, theta, p, y, X);

    for (int i = 0; i < max_iter; ++i) {
      T phi_cur = fun(theta + alpha_cur * p, y, X);
      T phi_last = fun(theta + alpha_last * p, y, X);
     
      if ((phi_cur > (fun(theta, y, X) + c1 * alpha_cur * fgradalpha0))
          ||
          ((phi_cur >= phi_last) && (i > 0))) {
        T alphastar = zoom(fun, alpha_last, alpha_cur, theta, p, y, X);
        return alphastar;
      }

      T fgradalpha_cur = gradfdifls(fun, alpha_cur, theta, p, y, X);
      if ( ::fabs(fgradalpha_cur) <= -1 * c2 * fgradalpha0)
        return alpha_cur;

      if ( fgradalpha_cur >= (T) 0.0) {
        T alphastar = zoom(fun, alpha_cur, alpha_last, theta, p, y, X);
        return alphastar;
      }
      
      alpha_last = alpha_cur;
      alpha_cur = runif() * (alpha_max - alpha_cur) + alpha_cur;
    }

    return 0.001;
  }


  /* Finds the minimum of a function once bracketed (i.e. over a 
   * closed interval).
   * Notes: Algorithm taken from Nocedal and Wright. 1999. Algorithm
   * 3.3
   */
  template <class T>
  T
  zoom (T (*fun)(const Matrix<T> &, const Matrix<T> &,
        const Matrix<T> &), const T &alo, const T &ahi,
        const Matrix<T> &theta, const Matrix<T> &p,
        const Matrix<T> &y, const Matrix<T> &X )
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
    if (! p.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "p not column vector");

    T alpha_lo = alo;
    T alpha_hi = ahi;
    T alpha_j = (alo + ahi) / 2.0;
    T phi_0 = fun(theta, y, X);
    T c1 = (T) 1e-4;
    T c2 = (T) 0.5;
    T fgrad0 = gradfdifls(fun,(T) 0, theta, p, y, X);

    int count = 0;
    int maxit = 20;
    while(count < maxit) {
      T phi_j = fun(theta + alpha_j * p, y, X);
      T phi_lo = fun(theta + alpha_lo * p, y, X);
     
      if ((phi_j > (phi_0 + c1 * alpha_j * fgrad0))
          || (phi_j >= phi_lo)){
        alpha_hi = alpha_j;
      } else {
        T fgradj = gradfdifls(fun, alpha_j, theta, p, y, X);
        if (::fabs(fgradj) <= -1 * c2 * fgrad0){ 
          return alpha_j;
        }
        if ( fgradj * (alpha_hi - alpha_lo) >= 0){
          alpha_hi = alpha_lo;
        }
        alpha_lo = alpha_j;
      }
      ++count;
    }
   
    return alpha_j;
  }



  /* Find the minimum of a function using the BFGS algorithm
   * Notes: Algorithm taken from Nocedal and Wright. 1999. algorithm
   * 8.1
   */
  template <class T>
  Matrix<T>
  BFGS (T (*fun)(const Matrix<T> &, const Matrix<T> &,
     const Matrix<T> &), const Matrix<T> &theta, 
  const Matrix<T> &y, const Matrix<T> &X, const int &maxit = 1000,
  const T &tolerance = 1e-6)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
    
    int n = theta.size();

    // H is initial inverse hessian
    Matrix<T> H = inv(hesscdif(fun, theta, y, X));

    // gradient at starting values
    Matrix<T> fgrad = gradfdif(fun, theta, y, X);
    Matrix<T> thetamin = theta;
    Matrix<T> fgrad_new = fgrad;
    Matrix<T> I = eye<T>(n); 

    int count = 0;
    while( (t(fgrad_new)*fgrad_new)[0] > tolerance) {
      Matrix<T> p = -1 * H * fgrad;
      T alpha = linesearch2(fun, thetamin, p, y, X);
      Matrix<T> thetamin_new = thetamin + alpha*p;
      fgrad_new = gradfdif(fun, thetamin_new, y, X);
      Matrix<T> s = thetamin_new - thetamin;
      Matrix<T> y = fgrad_new - fgrad;
      T rho = 1.0 / (t(y) * s)[0];
      H = (I - rho * s * t(y)) * H *(I - rho * y * t(s))
        + rho * s * (!s);

      thetamin = thetamin_new;
      fgrad = fgrad_new;
      ++count;

      std::cout << "BFGS iteration = " << count << std::endl;
      std::cout << "thetamin = " << (!thetamin).toString() << std::endl;
      std::cout << "gradient = " << (!fgrad).toString() << std::endl;
      std::cout << "t(gradient) * gradient = " << 
        ((!fgrad) * fgrad).toString() << std::endl;

      if (count > maxit)
        throw scythe_convergence_error (__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Failed to converge.  Try better starting values");
    }
   
    return thetamin;
  }


  /* Zero a function using Broyen's Method
   * Notes: Algorithm taken from Nocedal and Wright. 1999. algorithm 11
   * line search is not used to determine alpha (this should probably
   * be changed at some point.
   */
  template <class T>
  Matrix<T>
  nls_broyden(Matrix<T> (*fun)(const Matrix<T> &),
        const Matrix<T> &theta,  const int &maxit = 5000,
        const T &tolerance = 1e-6)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");

    Matrix<T> thetastar = theta;
    Matrix<T> B = jacfdif(fun, thetastar);

    Matrix<T> fthetastar;
    Matrix<T> p;
    Matrix<T> thetastar_new;
    Matrix<T> fthetastar_new;
    Matrix<T> s;
    Matrix<T> y;

    for (int i = 0; i < maxit; ++i) {
      fthetastar = fun(thetastar);
      p = lu_solve(B, -1 * fthetastar);
      T alpha = (T) 1.0;
      thetastar_new = thetastar + alpha*p;
      fthetastar_new = fun(thetastar_new);
      s = thetastar_new - thetastar;
      y = fthetastar_new - fthetastar;
      B = B + ((y - B * s) * (!s)) / ((!s) * s);
      thetastar = thetastar_new;
      if (max(fabs(fthetastar_new)) < tolerance)
        return thetastar;
    }
    
    throw scythe_convergence_error (__FILE__, __PRETTY_FUNCTION__,
      __LINE__,std::string("Failed to converge.  Try better starting") &
            " values or increase maxit");
  }


  /* Zero a function using Broyen's Method
   * Notes: Algorithm taken from Nocedal and Wright. 1999. algorithm
   * 11.3
   * line search is not used to determine alpha (this should probably
   * be changed at some point.
   */
  template <class T>
  Matrix<T>
  nls_broyden(Matrix<T> (*fun)(const Matrix<T> &, const Matrix<T> &), 
        const Matrix<T> &theta,  const Matrix<T> &psi, 
        const int &maxit=5000, const T &tolerance=1e-6)
  {
    if (! theta.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Theta not column vector");
    if (! psi.isColVector())
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Psi not column vector");

    Matrix<T> thetastar = theta;
    Matrix<T> B = jacfdif(fun, thetastar, psi);

    Matrix<T> fthetastar;
    Matrix<T> p;
    Matrix<T> thetastar_new;
    Matrix<T> fthetastar_new;
    Matrix<T> s;
    Matrix<T> y;

    for (int i = 0; i < maxit; ++i) {
      fthetastar = fun(thetastar, psi);
      p = lu_solve(B, -1 * fthetastar);
      T alpha = (T) 1.0;
      thetastar_new = thetastar + alpha*p;
      fthetastar_new = fun(thetastar_new, psi);
      s = thetastar_new - thetastar;
      y = fthetastar_new - fthetastar;
      B = B + ((y - B * s) * (!s)) / ((!s) * s);
      thetastar = thetastar_new;
      if (max(fabs(fthetastar_new)) < tolerance)
        return thetastar;
    }
    
    throw scythe_convergence_error (__FILE__, __PRETTY_FUNCTION__,
      __LINE__,std::string("Failed to converge.  Try better starting") &
            " values or increase maxit");
  }

}  // namespace dec

#endif /* SCYTHE_OPTIMIZE_CC */
