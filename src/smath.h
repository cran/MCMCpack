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
 * scythestat/math.h
 *
 * Provides definitions for the template wrapper functions
 * that allow common math.h operations to be performed on
 * Scythe matrices.
 *
 */

#ifndef SCYTHE_MATH_H
#define SCYTHE_MATH_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#else
#include "scythestat/matrix.h"
#endif

/* This portion of the library mimics math.h for Matrix<T> objects.  T
 * should be a floating point type of either float, double or long
 * double.  Using ints or other objects will cause warnings or errors.
 * You'll also need a sufficiently modern c++ compiler.  Note that the
 * f and l versions of these functions do not exist in c++ (you can
 * use the c calls if you wish) because type promotion takes care of
 * these functions.
 *
 * NOTE:  When I refer to x and y in the documentation below it means
 * the first and second matrix arguements respectively.  Also, I will
 * typically just write x when I mean forall elements in x
 */

namespace SCYTHE {

  /* acos - inverse cosine function */
  template <class T>
  Matrix<T> acos (Matrix<T>);
  
//  /* acosh - inverse hyperbolic cosine function */
//  template <class T>
//  Matrix<T> acosh (Matrix<T>);

  /* asin - inverse sine function */
  template <class T>
  Matrix <T> asin (Matrix<T>);

//  /* asinh - inverse hyperbolic sine function */
//  template <class T>
//  Matrix<T> asinh (Matrix<T>);

  /* atan - inverse tangent function */
  template <class T>
  Matrix<T> atan (Matrix<T>);
  
//  /* atanh - inverse hyperbolic tangent function */
//  template <class T>
//  Matrix<T> atanh (Matrix<T>);
  
  /* atan2 - returns angle whose tangent is y/x in the full angular
   * range [-pit,+pi].  Domain error if both x and y zero
   * The two matrices must have equal dimensions or one of the two
   * matrices must be scalar
   */
  template <class T>
  Matrix<T> atan2 (const Matrix<T> &, const Matrix<T> &);
  
  /* cbrt - cube root */
  template <class T>
  Matrix<T> cbrt (Matrix<T>);
  
  /* ceil - ceiling of a floating point number */
  template <class T>
  Matrix<T> ceil (Matrix<T>);
  
  /* copysign - return values with absval of 1st arg but sign of 2nd
   * The two matrices must have equal dimensions or one of the two
   * matrices must be scalar
   */
  template <class T>
  Matrix<T> copysign (const Matrix<T> &, const Matrix<T> &);

  /* cos - cosine function */
  template <class T>
  Matrix<T> cos (Matrix<T>);

  /* cosh - hyperbolic cosine function */
  template <class T>
  Matrix<T> cosh (Matrix<T>);

  /* erf - error function */
  template <class T>
  Matrix<T> erf (Matrix<T>);

  /* erfc - complementary error function */
  template <class T>
  Matrix<T> erfc (Matrix<T>);

  /* exp - Calculate the value of e^x for each  individual */
  template <class T>
  Matrix<T> exp (Matrix<T>);

//  /* expm1 - exponent minus 1 */
//  template <class T>
//  Matrix<T> expm1 (Matrix<T>);
  
  /* fabs - Calculate the absolute value of each Matrix element */
  template <class T>
  Matrix<T> fabs (Matrix<T>);

  /* floor - floor of floating point number */
  template <class T>
  Matrix<T> floor (Matrix<T>);

  /* fmod - return the remainder */
  template <class T>
  Matrix<T> fmod (const Matrix<T> &, const Matrix<T> &);

  /* frexp - returns fractional value of input, and fills int matrix
   * with exponent ex.  Frac on interval [1/2,1) x == frac * 2^ex
   */
  template <class T>
  Matrix<T> frexp(Matrix<T>, Matrix<int> &);

  /* hypot - euclidean distance function */
  template <class T>
  Matrix<T> hypot (const Matrix<T> &, const Matrix<T> &);
  
  /* ilogb - returns int verison of logb */
  template <class T>
  Matrix<int> ilogb (const Matrix<T> &);

  /* j0, j1, jn - bessel functions of the first kind of order
   * (May only support doubles, consult standard)
   */
  template <class T>
  Matrix<T> j0 (Matrix<T>);

  template <class T>
  Matrix<T> j1 (Matrix<T>);

  template <class T>
  Matrix<T> jn (const int &n, Matrix<T>);

  /* ldexp - returns x * 2^ex */
  template <class T>
  Matrix<T> ldexp(Matrix<T>, const int &);

  /* lgamma - returns natural log of the absval of the gamma function */
  template <class T>
  Matrix<T> lgamma (Matrix<T>);

  /* Log - Calculate the natural log of each Matrix element */
  template <class T>
  Matrix<T> log(Matrix<T>);
  
  /* Log10 - Calculate the Base 10 Log of each Matrix element */
  template <class T> 
  Matrix<T> log10(Matrix<T>);

  /* log1p - returns natrual log of 1 + x, domain error if x < -1 */
  template <class T>
  Matrix<T> log1p (Matrix<T>);

  /* logb - returns ex s.t. x == frac * ex^FLT_RADIX where frac is on
   * the interval [1,FLT_RADIX].  Domain error if x is 0.
   */
  template <class T>
  Matrix<T> logb (Matrix<T>);

  /* modf - x == frac + i where |frac| on [0,1) and both frac and i
   * have the same sign as x.  I is stored in the second matrix.
   */
  template <class T>
  Matrix<T> modf (Matrix<T>, Matrix<double> &);
  
  /* Pow - Raise each Matrix element to the power of the 2nd arg
   */
  template <class T, class S> 
  Matrix<T> pow(Matrix<T>, const S &);

  /* return the remainder of dividing */
  template <class T>
  Matrix<T> remainder (const Matrix<T> &, const Matrix<T> &);

  /* rint - returns x round to the nearest int using the current
   * rounding mode.  May rais and inexact floating-point exception if
   * the return value does not equal x???
   */
  template <class T>
  Matrix<T> rint (Matrix<T>);

  /* scalbn - returns x * FLT_RADIX^ex (ex is 2nd arg) */
  template <class T>
  Matrix<T> scalbn (Matrix<T>, const int &);

  /* sin - return the sine of x */
  template <class T>
  Matrix<T> sin (Matrix<T>);

  /* sinh - return the hyperbolic sine of x */
  template <class T>
  Matrix<T> sinh (Matrix<T>);

  /* Sqrt - Calculate the sqrt of each element of a Matrix */
  template <class T> 
  Matrix<T> sqrt (Matrix<T>);

  /* tan - return the tangent of x */
  template <class T>
  Matrix<T> tan (Matrix<T>);
  
  /* tanh - return the hyperbolic tangent of x */
  template <class T>
  Matrix<T> tanh (Matrix<T>);

  /* y0, y1, yn - bessel functions of the second kind of order
   * (May only support doubles, consult standard)
   */
  template <class T>
  Matrix<T> y0 (Matrix<T>);

  template <class T>
  Matrix<T> y1 (Matrix<T>);

  template <class T>
  Matrix<T> yn (const int &, Matrix<T>);


} // end namespace SCYTHE

#if defined (SCYTHE_COMPILE_DIRECT) && \
	  (defined (__GNUG__) || defined (__MWERKS__) || \
		 defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION))
#include "smath.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_MATH_H */
