/* Scythe_IDE.h
 *
 * This header provides definitions for inversion and decomposition
 * template functions for the Scythe Statistical Library.
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

#ifndef SCYTHE_IDE_H
#define SCYTHE_IDE_H

#include "Scythe_Matrix.h"

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
  Matrix<T> inv(Matrix<T>);

  /* Calculates the determinant of Matrix A via LU Decomposition */
  template <class T>
  T det(const Matrix<T> &);

}  // end namespace SCYTHE
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || \
    defined (EXPLICIT_TEMPLATE_INSTANTIATION)
  // Necessary for template instantiation with some compilers.
# include "Scythe_IDE.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_IDE_H */
