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
 * scythestat/ide.h
 *
 * Provides definitions for inversion and decomposition
 * template functions that operate on Scythe's Matrix class.
 *
 */

#ifndef SCYTHE_IDE_H
#define SCYTHE_IDE_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#else
#include "scythestat/matrix.h"
#endif

namespace SCYTHE {

  /* Cholesky decomposition of a sym pos-def matrix */
  template <class T>
  Matrix<T>
  cholesky (const Matrix<T> &);

  /* Solves Ax=b for x via backsubstitution using Cholesky
   * Decomposition  (NOTE: function is overloaded) A must be symmetric
   * and positive definite
   */
  template <class T>
  Matrix<T>
  chol_solve (const Matrix<T> &, const Matrix<T> &);

  /* Solves Ax=b for x via backsubstitution using Cholesky
   * Decomposition. This function takes in the lower triangular L as
   * input and does not depend upon cholesky() A must be symmetric and
   * positive definite
   */
  template <class T>
  Matrix<T> chol_solve (const Matrix<T> &, const Matrix<T> &,
                        const Matrix<T> &);

  /* Calculates the inverse of a Sym. Pos. Def. Matrix (NOTE: function
   * is overloaded)
   */
  template <class T>
  Matrix<T> invpd (const Matrix<T> &);

  /* Calculates the inverse of a Sym. Pos. Def. Matrix (NOTE: function
   * is overloaded)
   */
  template <class T>
  Matrix<T> invpd (const Matrix<T> &, const Matrix<T> &);

  /* Calculates the LU Decomposition of a square Matrix */
  template <class T>
  void lu_decomp (Matrix<T>, Matrix<T> &, Matrix<T> &,
                  Matrix<int> &);  

  /* Solve Ax=b for x via forward and backsubstitution using the LU
   * Decomp of Matrix A (NOTE: This function is overloaded)
   */
  template <class T> 
  Matrix<T> lu_solve(Matrix<T>, const Matrix<T> &);

  /* Solve Ax=b for x via forward and backsubstitution using the LU
   * Decomp of Matrix A (NOTE: This function is overloaded)
   */
  template <class T>
  Matrix<T> lu_solve (Matrix<T>, const Matrix<T> &, const Matrix<T> &,
                      const Matrix<T> &, const Matrix<int> &);

  /* Interchanges the rows of A with those in vector p and returns the
   * modified Matrix.
   */
  template <class T>
  Matrix<T> row_interchange(Matrix<T>, const Matrix<int> &);

  /* Calculate the Inverse of a square Matrix A via LU decomposition 
   *
   * DEPRECATED:  see operator^= in Scythe_Matrix
   */
  template <class T>
  inline Matrix<T> inv(Matrix<T> A) {
		return (A ^= -1);
	}

  /* Calculates the determinant of Matrix A via LU Decomposition */
  template <class T>
  inline T det(const Matrix<T> &A)
	{
		return ~A;
	}

}  // end namespace SCYTHE

#if defined (SCYTHE_COMPILE_DIRECT) && \
	  (defined (__GNUG__) || defined (__MWERKS__) || \
		 defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION))
#include "ide.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_IDE_H */
