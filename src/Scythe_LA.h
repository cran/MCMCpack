/* Scythe_LA.h
 *
 * This header provides definitions of functions which deal with
 * linear algebra in the Scythe Statistical Library.
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

#ifndef SCYTHE_LA_H
#define SCYTHE_LA_H

#include "Scythe_Matrix.h"

namespace SCYTHE {

  /* Transpose  - computes the transpose of a Matrix */
  template <class T>
   Matrix<T> t (const Matrix<T> &);

  /* Ones - creates a Matrix of ones */
  template <class T>
   Matrix<T> ones (const int &, const int &);

  /* Eye - creates an Identity Matrix of size k x k */
  template <class T>
  Matrix<T> eye (const int &);

  /* Seqa - creates a vector additive sequence Matrix (size x 1) */
  template <class T>
  Matrix<T> seqa (T, const T &, const int &);

  /* sort - sorts all elements of a Matrix in row_major order.  This
   * function is DEPRECATED.  It is simply a wrapper to the STL sort
   * algorithm which you should use instead.
   * To sort with the STL in row_major simply do
   * sort(M.begin(), M.end());
   * To sort with the STL in col_major simply do
   * sort(M.beginc(), M.endc());
   *
   * Version two takes a compare function object while the first uses
   * primitive comparison operators.
   */
  template <class T> 
  Matrix<T> sort (Matrix<T>);

  /* sortc - sorts all columns of a Matrix using the STL sort
   * algorithm.  The second version takes a compare function object
   * while the first uses primitive comparison operators
   *
   * Unlike sort, this function is not deprecated
   */
  template <class T>
  Matrix<T> sortc (Matrix<T> A);

  /* Cbind - Column bind 2 matrices */
  template <class T>
  Matrix<T> cbind (const Matrix<T> &, const Matrix<T> &);

  /* FUNCTION: Rbind - Row bind 2 matrices */
  template <class T>
  Matrix<T> rbind (const Matrix<T> &, const Matrix<T> &);

  /* Order - Calculates the order of each element in a Matrix */
  // XXX - ask Quinn about this one
  template <class T>
  Matrix<int> order(const Matrix<T> &);

  /* Selif - Selects all the rows of Matrix A for which the col vector
   * has an element equal to 1 */
  template <class T>
  Matrix<T> selif(const Matrix<T> &, const Matrix<bool> &);

  /* Unique - Finds unique elements in a Matrix */
  template <class T>
  Matrix<T> unique(const Matrix<T> &);

  /* Vecr - Turn Matrix into Column vector by stacking rows */
  template <class T>
  inline Matrix<T> vecr(const Matrix<T> &A)
  {
    return (Matrix<T> (A.size(), 1, A.getArray()));
  }

  /* Vecc - Turn Matrix into Column vector by stacking columns */
  template <class T>
  inline Matrix<T> vecc(const Matrix<T> &A)
  {
    Matrix<T> temp(A.size(), 1, false);

    // Note we can use a row_major_iterator to write because we are
    // writing to a vector.  RMIs are a bit faster than CMIs.
    copy(A.beginc(), A.endc(), temp.begin());

    return temp;
  }

  /* Reshape - Reshapes a row major order Matrix or Vector */
  template <class T>
  Matrix<T> reshape(const Matrix<T> &, const int &, const int &);

  /* Vech - Make vector out of unique elements of a symmetric */
  template <class T>
  Matrix<T> vech(const Matrix<T> &);

  /* Xpnd - Get symmetric Matrix B back from A = vech(B) */
  template <class T>
  Matrix<T> xpnd(const Matrix<T> &);
  
  /* Diag - get the diagonal of a Matrix */
  template <class T>
  Matrix<T> diag(const Matrix<T> &);
  
  /* Gaxpy - Fast calculation of A*B + C */
  template <class T>
  Matrix<T> gaxpy(const Matrix<T> &, const Matrix<T> &,
                    const Matrix<T> &);
  
  /* Crossprod - Fast calculation of A'A */
   template <class T>
   Matrix<T> crossprod(const Matrix<T> &);

} // end namespace SCYTHE
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || \
    defined (EXPLICIT_TEMPLATE_INSTANTIATION)
  // Necessary for template instantiation with some compilers.
# include "Scythe_LA.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_LA_H */
