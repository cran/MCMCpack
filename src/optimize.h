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
 * sythestat/optimization.h
 *
 * Provides definitions for various numerical optimization
 * routines.
 *
 */

#ifndef SCYTHE_OPTIMIZE_H
#define SCYTHE_OPTIMIZE_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#include "rng.h"
#else
#include "scythestat/matrix.h"
#include "scythestat/rng.h"
#endif

// Avoid NameSpace Pollution
namespace SCYTHE {
  
  /* returns the machine epsilon (float, double, or long double) */
  template <class T>
  T epsilon();
 
  /* Caculates the definite integral of a function from a to b, given
   * the function, a, b, and the number of subintervals
   */
  template <class T>
  T  intsimp(T (*fun)(const T &), const T &, const T &, const int &);
  
  /* Calculates the definite integral of a function from a to b, given
   * the function, a, b, the  of subintervals, and a tolerance
   */
  template <class T>
  T adaptsimp(T (*fun)(const T &), const T &, const T &, const int &, 
              const T &tol = 1e-5);
  
  /* Numerically calculates the gradient of a function at theta using
   * a forward difference formula, given the function, theta (col
   * vector), and two matrix args to be sent to the function.  (The
   * function takes three matrices)
   */
  template <class T>
  Matrix<T> gradfdif (T (*fun)(const Matrix<T> &, const Matrix<T> &,
                      const Matrix<T> &), const Matrix<T> &,
                      const Matrix<T> &, const Matrix<T> &);
  
  /* Numerically calculates the first deriv.of a function wrt alpha at
   * (theta + alpha *p) using a forward difference formula, given the
   * function, alpha, theta, p, and two matrices describing the
   * function.  (Primarily useful in linesearches)
   */
  template <class T>
  T gradfdifls (T (*fun)(const Matrix<T> &, const Matrix<T> &,
                const Matrix<T> &), const T &,
                const Matrix<T> &, const Matrix<T> &,
                const Matrix<T> &, const Matrix<T> &);
  
  /* Numerically calculates the jacobian of a function at theta using
   * a forward difference formula, given the function and theta
   */
  template <class T>
  Matrix<T> jacfdif(Matrix<T> (*fun)(const Matrix<T> &),
                    const Matrix<T> &);
  
  /* Numerically calculates the Jacobian of a function a theta using a
   * forward difference formula given the function, theta, and psi ( a
   * column vector of parameter values at which to calculate the
   * jacobian)
   */
  template <class T>
  Matrix<T> jacfdif(Matrix<T> (*fun)(const Matrix<T> &,
                    const Matrix<T> &), const Matrix<T> &,
                    const Matrix<T> &);

  /* Numerically calculates the Hessian of a function at theta using a
   * central difference formula given the function, theta, and two
   * matrix arguments for the function
   */
  template <class T>
  Matrix<T> hesscdif (T (*fun)(const Matrix<T> &, const Matrix<T> &,
                      const Matrix<T> &), const Matrix<T> &,
                      const Matrix<T> &, const Matrix<T> &);
  
  /* Performs a line search to find the step length alpha that
   * approximately minimizes an implied 1d function, given the
   * function to minimize, col-vector theta, and two matrix arguments
   * for the function
   */
  template <class T>
  T linesearch1(T (*fun)(const Matrix<T> &, const Matrix<T> &, 
                const Matrix<T> &), const Matrix<T> &,
                const Matrix<T> &, const Matrix<T> &,
                const Matrix<T> &);
  
  /* Performs a line search to find the step length alpha that
   * approximately minimizes an implied 1d function, given the
   * function, theta, direction vec p, and two matrix args
   */
  template <class T>
  T linesearch2(T (*fun)(const Matrix<T> &, const Matrix<T> &, 
                const Matrix<T> &), const Matrix<T> &,
                const Matrix<T> &, const Matrix<T> &,
                const Matrix<T> &, rng *);
  
  /* Finds minimum of a function once bracketed, given the function,
   * lower bracket, upper bracket, theta, direction vector p, and to
   * matrix arguments
   */
  template <class T>
  T zoom (T (*fun)(const Matrix<T> &, const Matrix<T> &,
          const Matrix<T> &), const T &, const T &, const Matrix<T> &,
          const Matrix<T> &, const Matrix<T> &, const Matrix<T> &);
  
  
  /* Numerically finds the minimum of a function using the BFGS
   * algorithm, given the function, theta, and two matrix args
   */
  template <class T>
  Matrix<T> BFGS (T (*fun)(const Matrix<T> &, const Matrix<T> &,
                  const Matrix<T> &), rng *, const Matrix<T> &,
                  const Matrix<T> &y = Matrix<T>(1,1),
                  const Matrix<T> & X = Matrix<T>(1,1),
                  const int &maxit=1000, const T &tolerance=1e-4);
  
  /* Solves a system of n nonlinear equations in n unknowns of the form
   * fun(thetastar) = 0 for thetastar given the function, starting 
   * value theta, max number of iterations, and tolerance.
   * Uses Broyden's method.
   */
  template <class T>
  Matrix<T> nls_broyden(Matrix<T> (*fun)(const Matrix<T> &),
                        const Matrix<T> &, const int &maxit=5000,
                        const T &tolerance=1e-6);

  /* Nls_broyden  - Solves a system of n nonlinear equations in n
   * unknowns of the form:
   * fun(thetastar) = 0 
   * for thetastar given, the function, the  starting value theta,
   * matrix of fixed parameters psi, max iteration, and tolerance.
   * Uses Broyden's method.
   */
  template <class T>
  Matrix<T> nls_broyden(Matrix<T> (*fun)(const Matrix<T> &,
                        const Matrix<T> &), const Matrix<T> &,
                        const Matrix<T> &, const int &maxit = 5000,
                        const T& tolerance = 1e-6);
  
}  // namespace dec

#if defined (SCYTHE_COMPILE_DIRECT) && \
	  (defined (__GNUG__) || defined (__MWERKS__) || \
		 defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION))
#include "optimize.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_OPTIMIZE_H */
