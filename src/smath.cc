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
 * scythestat/math.cc
 *
 * Provides implementations of the template wrapper functions
 * that allow common math.h operation to be performed on
 * Scythe matrices.
 *
 */

#ifndef SCYTHE_MATH_CC
#define SCYTHE_MATH_CC

#include <cmath>

#ifdef SCYTHE_COMPILE_DIRECT
#include "smath.h"
#include "error.h"
#include "util.h"
#else
#include "scythestat/smath.h"
#include "scythestat/error.h"
#include "scythestat/util.h"
#endif

namespace SCYTHE {

  /* calc the inverse cosine of each element of a Matrix */
  template <class T>
  Matrix<T>
  acos (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::acos(A[i]);

    return A;
  }
  
//  /* calc the inverse hyperbolic cosine of each element of a Matrix */
//  template <class T>
//  Matrix<T>
//  acosh (Matrix<T> A)
//  {
//    for (int i = 0; i < A.size(); ++i)
//      A[i] = ::acosh(A[i]);
//
//    return A;
//  }
  
  /* calc the inverse sine of each element of a Matrix */
  template <class T>
  Matrix<T>
  asin (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::asin(A[i]);

    return A;
  }
  
//  /* calc the inverse hyperbolic sine of each element of a Matrix */
//  template <class T>
//  Matrix<T>
//  asinh (Matrix<T> A)
//  {
//    for (int i = 0; i < A.size(); ++i)
//      A[i] = ::asinh(A[i]);
//
//    return A;
//  }
  
  /* calc the inverse tangent of each element of a Matrix */
  template <class T>
  Matrix<T>
  atan (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::atan(A[i]);

    return A;
  }
  
//  /* calc the inverse hyperbolic tangent of each element of a Matrix */
//  template <class T>
//  Matrix<T>
// atanh (Matrix<T> A)
//  {
//    for (int i = 0; i < A.size(); ++i)
//      A[i] = ::atanh(A[i]);
//
//    return A;
//  }
  
  /* calc the angle whose tangent is y/x  */
  template <class T>
  Matrix<T>
  atan2 (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<T> temp;
    
    if (A.isScalar()) {
      temp = B;
      for (int i = 0; i < B.size(); ++i)
        temp[i] = ::atan2(A[0], B[i]);
    } else if (B.isScalar()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::atan2(A[i], B[0]);
    } else if (A.size() == B.size()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::atan2(A[i], B[i]);
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "A.size() != B.size() and neither A nor B is scalar");
    }

    return temp;
  }
  
  /* calc the cube root of each element of a Matrix */
  template <class T>
  Matrix<T>
  cbrt (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::cbrt(A[i]);

    return A;
  }
  
  /* calc the ceil of each element of a Matrix */
  template <class T>
  Matrix<T>
  ceil (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::ceil(A[i]);

    return A;
  }
  
  /* create a matrix containing the absval of the first input and the
   * sign of the second
   */
  template <class T>
  Matrix<T>
  copysign (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<T> temp;
    
    if (A.isScalar()) {
      temp = B;
      for (int i = 0; i < B.size(); ++i)
        temp[i] = ::copysign(A[0], B[i]);
    } else if (B.isScalar()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::copysign(A[i], B[0]);
    } else if (A.size() == B.size()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::copysign(A[i], B[i]);
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "A.size() != B.size() and neither A nor B is scalar");
    }

    return temp;
  }
  
  /* calc the cosine of each element of a Matrix */
  template <class T>
  Matrix<T>
  cos (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::cos(A[i]);

    return A;
  }
  
  /* calc the hyperbolic cosine of each element of a Matrix */
  template <class T>
  Matrix<T>
  cosh (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::cosh(A[i]);

    return A;
  }
  
  /* calc the error function of each element of a Matrix */
  template <class T>
  Matrix<T>
  erf (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::erf(A[i]);

    return A;
  }
  
  /* calc the complementary error function of each element of a Matrix */
  template <class T>
  Matrix<T>
  erfc (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::erfc(A[i]);

    return A;
  }
  
  /* calc the vaue e^x of each element of a Matrix */
  template <class T>
  Matrix<T>
  exp (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::exp(A[i]);

    return A;
  }
  
//  /* calc the exponent - 1 of each element of a Matrix */
//  template <class T>
//  Matrix<T>
//  expm1 (Matrix<T> A)
//  {
//    for (int i = 0; i < A.size(); ++i)
//     A[i] = ::expm1(A[i]);
//
//    return A;
//  }
  
  /* calc the absval of each element of a Matrix */
  template <class T>
  Matrix<T>
  fabs (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::fabs (A[i]);

    return A;
  }

  /* calc the floor of each element of a Matrix */
  template <class T>
  Matrix<T>
  floor (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::floor(A[i]);

    return A;
  }
  
  /* calc the remainder of the division of each matrix element */
  template <class T>
  Matrix<T>
  fmod (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<T> temp;
    
    if (A.isScalar()) {
      temp = B;
      for (int i = 0; i < B.size(); ++i)
        temp[i] = ::fmod(A[0], B[i]);
    } else if (B.isScalar()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::fmod(A[i], B[0]);
    } else if (A.size() == B.size()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::fmod(A[i], B[i]);
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "A.size() != B.size() and neither A nor B is scalar");
    }

    return temp;
  }
  
  /* calc the fractional val of input and return exponents in int
   * matrix reference
   */
  template <class T>
  Matrix<T>
  frexp (Matrix<T> A, Matrix<int> &ex)
  {
    if (A.size() != ex.size())
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "The input matrix sizes do not match");

    for (int i = 0; i < A.size(); ++i)
      A[i] = ::frexp(A[i], &(ex[i]));

    return A;
  }

  /* calc the euclidean distance between the two inputs */
  template <class T>
  Matrix<T>
  hypot (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<T> temp;
    
    if (A.isScalar()) {
      temp = B;
      for (int i = 0; i < B.size(); ++i)
        temp[i] = ::hypot(A[0], B[i]);
    } else if (B.isScalar()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::hypot(A[i], B[0]);
    } else if (A.size() == B.size()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::hypot(A[i], B[i]);
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "A.size() != B.size() and neither A nor B is scalar");
    }

    return temp;
  }

  /*  return (int) logb */
  template <class T>
  Matrix<int>
  ilogb (const Matrix<T> &A)
  {
    Matrix<int> temp(A.rows(), A.cols(), false);
    
    for (int i = 0; i < A.size(); ++i)
      temp[i] = ::ilogb(A[i]);

    return temp;
  }
  
  /* compute the bessel func of the first kind of the order 0 */
  template <class T>
  Matrix<T>
  j0 (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::j0(A[i]);

    return A;
  }
  
  /* compute the bessel func of the first kind of the order 1 */
  template <class T>
  Matrix<T>
  j1 (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::j1(A[i]);

    return A;
  }
  
  /* compute the bessel func of the first kind of the order n */
  template <class T>
  Matrix<T>
  jn (const int &n, Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::jn(n, A[i]);

    return A;
  }

  /* calc x * 2 ^ex */
  template <class T>
  Matrix<T>
  ldexp (Matrix<T> A, const int &ex)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::ldexp(A[i], ex);

    return A;
  }
  
  /*  compute the natural log of the absval of gamma function */
  template <class T>
  Matrix<T>
  lgamma (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::lgamma(A[i]);

    return A;
  }
  
  /* calc the natural log of each element of a Matrix */
  template <class T>
  Matrix<T>
  log (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::log(A[i]);

    return A;
  }
  
  /* calc the base-10 log of each element of a Matrix */
  template <class T>
  Matrix<T>
  log10 (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::log10(A[i]);

    return A;
  }
  
  /* calc the natural log of 1 + each element of a Matrix */
  template <class T>
  Matrix<T>
  log1p (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::log1p(A[i]);

    return A;
  }
  
  /* calc the logb of each element of a Matrix */
  template <class T>
  Matrix<T>
  logb (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::logb(A[i]);

    return A;
  }
  
  /* x = frac + i, return matrix of frac and place i in 2nd matrix
   */
  template <class T>
  Matrix<T>
  modf (Matrix<T> A, Matrix<double> &iret)
  {
    if (A.size() != iret.size())
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "The input matrix sizes do not match");

    for (int i = 0; i < A.size(); ++i)
      A[i] = ::modf(A[i], &(iret[i]));

    return A;
  }

  /* calc x^ex of each element of a Matrix */
  template <class T, class S>
  Matrix<T>
  pow (Matrix<T> A, const S &ex)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::pow(A[i], ex);

    return A;
  }

  /* calc rem == x - n * y */
  template <class T>
  Matrix<T>
  remainder (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<T> temp;
    
    if (A.isScalar()) {
      temp = B;
      for (int i = 0; i < B.size(); ++i)
        temp[i] = ::remainder(A[0], B[i]);
    } else if (B.isScalar()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::remainder(A[i], B[0]);
    } else if (A.size() == B.size()) {
      temp = A;
      for (int i = 0; i < A.size(); ++i)
        temp[i] = ::remainder(A[i], B[i]);
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "A.size() != B.size() and neither A nor B is scalar");
    }

    return temp;
  }

  /* return x rounded to nearest int */
  template <class T>
  Matrix<T>
  rint (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::rint(A[i]);

    return A;
  }

  /* returns x * FLT_RADIX^ex */
  template <class T>
  Matrix<T>
  scalbn (Matrix<T> A, const int &ex)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::scalbn(A[i], ex);

    return A;
  }

  /*  calc the sine of x */
  template <class T>
  Matrix<T>
  sin (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::sin(A[i]);

    return A;
  }

  /* calc the hyperbolic sine of x */
  template <class T>
  Matrix<T>
  sinh (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::sinh(A[i]);

    return A;
  }
  
  /* calc the sqrt of x */
  template <class T>
  Matrix<T>
  sqrt (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::sqrt(A[i]);

    return A;
  }


  /* calc the tangent of x */
  template <class T>
  Matrix<T>
  tan (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::tan(A[i]);

    return A;
  }

  /* calc the hyperbolic tangent of x */
  template <class T>
  Matrix<T>
  tanh (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::tanh(A[i]);

    return A;
  }

  /* bessel function of the second kind of order 0*/
  template <class T>
  Matrix<T>
  y0 (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::y0(A[i]);

    return A;
  }

  /* bessel function of the second kind of order 1*/
  template <class T>
  Matrix<T>
  y1 (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::y1(A[i]);

    return A;
  }

  /* bessel function of the second kind of order n*/
  template <class T>
  Matrix<T>
  yn (const int &n, Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = ::yn(n, A[i]);

    return A;
  }
  
} // end namespace SCYTHE

#ifndef SCYTHE_COMPILE_DIRECT
#include "scythestat/eti/smath.t"
#endif

#endif /* SCYTHE_MATH_CC */
