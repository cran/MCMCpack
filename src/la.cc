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
 * scythestat/la.cc
 *
 * Provides implementations of functions that perform common
 * linear algebra manipulation on Scythe matrices.
 *
 */

#ifndef SCYTHE_LA_CC
#define SCYTHE_LA_CC

#include <cmath>
#include <algorithm>
#include <numeric>
#include <set>

#ifdef SCYTHE_COMPILE_DIRECT
#include "error.h"
#include "util.h"
#include "la.h"
#include "stat.h"
#else
#include "scythestat/error.h"
#include "scythestat/util.h"
#include "scythestat/la.h"
#include "scythestat/stat.h"
#endif

namespace SCYTHE {
  /* Compute the transpose of a matrix. Kept for back-compatibility*/
  template <class T>
  Matrix<T>
  t (const Matrix<T> &old_matrix)
  {
    return (! old_matrix);
  }
  
  /* Create a matrix of ones from the given dimensions 
   * Note:  call is of from ones<double>(4,3) or ones<int>(5,8)
   */
  template <class T>
  Matrix<T>
  ones (const int& rows, const int& cols)
  {
    if (rows < 1 || cols < 1) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, std::string("Improper row (") & rows
          & ") or column (" & cols & ") dimension");
    }
    
    return Matrix<T> (rows, cols, true, (T) 1);
  }
  
  /* Create a k x k identity matrix
   * Note:  class is of form eye<double>(4) or eye<int>(7
   */
  template <class T>
  Matrix<T>
  eye (const int &k)
  {
    Matrix<T> temp(k, k, false);
    for (int i = 0; i < temp.rows(); ++i) {
      for (int j = 0; j < temp.cols(); ++j) {
        if (i == j)
          temp(i,j) = (T) 1.0;
        else
          temp(i,j) = (T) 0.0;
      }
    }
    
    return temp;
  }
  
  /* Create a k x 1 vector-additive sequence matrix */
  template <class T>
  Matrix<T>
  seqa (T start, const T& incr, const int& size)
  {
    Matrix<T> temp (size, 1, false);
    for (int i = 0; i < size; ++i) {
      temp[i] = start;
      start += incr;
    }
    
    return temp;
  }
  
  /* Uses the STL sort to sort a Matrix in ascending row-major order */
  template <class T>
  Matrix<T>
  sort (Matrix<T> A) {
    sort(A.begin(), A.end());
    return A;
  }
    
  template <class T>
  Matrix<T>
  sortc (Matrix<T> A)
  {
    for (typename Matrix<T>::col_major_iterator it = A.beginc();
        it < A.endc(); it.next_vec())
      sort(it, it + A.rows());

    return A;
  }
  
  /* Column bind two matrices */
  template <class T>
  Matrix<T>
  cbind (const Matrix<T> &A, const Matrix<T> &B)
  {
    if (A.rows() != B.rows()) {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrices have different number of rows");
    }
    
    Matrix<T> C(A.rows(), A.cols() + B.cols(), false);
    typename Matrix<T>::col_major_iterator write = C.beginc();
  
    for (typename Matrix<T>::const_col_major_iterator read = A.beginc();
        read < A.endc(); ++read)
      *(write++) = *read;
  
    for (typename Matrix<T>::const_col_major_iterator read = B.beginc();
        read < B.endc(); ++read)
      *(write++) = *read;
  
    return C;
  }

  
  /* Row bind two matrices: kept for backwards compatibility */
  template <class T>
  Matrix<T>
  rbind (const Matrix<T> &A, const Matrix<T> &B)
  {
    if (A.cols() != B.cols()) {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrices have different number of rows");
    }
    
    Matrix<T> C(A.rows() + B.rows(), A.cols(), false);
    typename Matrix<T>::row_major_iterator write = C.begin();
  
    for (typename Matrix<T>::const_row_major_iterator read = A.begin();
        read < A.end(); ++read)
      *(write++) = *read;
  
    for (typename Matrix<T>::const_row_major_iterator read = B.begin();
        read < B.end(); ++read)
      *(write++) = *read;
  
    return C;
  }
  
  /* Calculates the order of each element in a Matrix */
  template <class T>
  Matrix<int>
  order(const Matrix<T> &A){
    if (! A.isColVector()) {
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "A not a column vector");
    }
    Matrix<int> temp(A.rows(), 1, false);
    for (int i = 0;  i < A.rows(); ++i) {
      temp[i] = sumc(A << A[i])[0];
    }
    
    return temp;
  }
  
  /* Selects all the rows of Matrix A for which binary column vector e
   * has an element equal to 1
   */
  template <class T>
  Matrix<T>
  selif(const Matrix<T> &A, const Matrix<bool> &e)
  {
    if (A.rows() != e.rows()) {
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "A and e have different number of rows");
    }
  
    if (! e.isColVector()) {
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "e not a column vector");
    }
    
    // See how many rows are true
    int N = std::accumulate(e.begin(), e.end(), (int) 0);
    
    // declare and form output Matrix
    Matrix<T> temp(N, A.cols(), false);
    int cnt = 0;
    for (int i = 0; i < e.size(); ++i) {
      if (e[i])  {
        copy(A.vec(i), A.vec(i + 1), temp.vec(cnt++));
      }
    }
  
    return temp;
  }
  
  /* Find unique elements in a matrix and return a sorted row vector */
  template <class T>
  Matrix<T>
  unique(const Matrix<T> &A)
  {
    std::set<T> u(A.begin(), A.end());
    Matrix<T> temp(1, u.size(), false);
    
    copy(u.begin(), u.end(), temp.begin());

    return temp;
  }

  /* Reshape a matrix */
  template <class T>
  Matrix<T>
  reshape(const Matrix<T> &A, const int &r, const int &c) 
  {
    if (A.size() != r * c)
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
          std::string("Input dimensions (") & r & "," & c & ") not" &
          " consistent with size of input matrix (" & A.size() & ")");
  
    Matrix<T> temp(r, c, A.getArray());
    return temp;
  }
  
  /* Make vector out of unique elements of a symmetric Matrix.  
   * NOTE: DOES NOT CHECK FOR SYMMETRY!!!
   */
  template <class T>
  Matrix<T>
  vech(const Matrix<T> &A)
  {
    if (! A.isSquare()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrix not square");
    }
    Matrix<T> temp ((int) (0.5 * (A.size() - A.rows())) + A.rows(), 1,
        false);
    typename Matrix<T>::row_major_iterator iter = temp.begin();
    
    for (int i = 0; i < A.rows(); ++i)
      iter = copy(A.vecc(i) + i, A.vecc(i + 1), iter);
    
    
    return temp;
  }
  
  /* Expand xpnd(A) == B from A = vech(B) */
  template <class T>
  Matrix<T>
  xpnd(const Matrix<T> &A)
  {
    double newrowsize_d = -.5 + .5 * ::sqrt(1 + 8 * A.size());
    if (std::fmod(newrowsize_d, 1.0) != 0.0)
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
          "Can't turn input vector into a square matrix");
    
    int newrowsize = (int) newrowsize_d;
    Matrix<T> temp(newrowsize, newrowsize, false);
    int cnt = 0;
  
    for (int i = 0; i < newrowsize; ++i) {
      for (int j = i; j < newrowsize; ++j)
        temp(i, j) = temp(j, i) = A[cnt++];
    }
  
    return temp;
  }
  
  /* Get the diagonal of a Matrix. */
  template <class T>
  Matrix<T>
  diag(const Matrix<T> &A)
  {
    if (A.rows() != A.cols())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrix not square");
  
    Matrix<T> temp(A.rows(), 1, false);
    for (int i = 0; i < A.rows(); ++i)
      temp[i] = A(i, i);
    
    return temp;
  }

  template <class T>
  Matrix<T>
  gaxpy (const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C)
  {
    Matrix<T> temp;

    if (A.isScalar() && B.rows() == C.rows() && B.cols() == C.cols()) {
      // Case 1: 1 x 1  *  n x k  +  n x k
      temp = Matrix<T> (B.rows(), B.cols(), false);
      
      for (int i = 0; i  < B.size(); ++i)
        temp[i] = A[0] * B[i] + C[i];
      
    } else if (B.isScalar() && A.rows() == C.rows() && 
        A.cols() == C.cols()) {
      // Case 2: m x n  *  1 x 1  +  m x n
      temp = Matrix<T> (A.rows(), A.cols(), false);
      
      for (int i = 0; i  < A.size(); ++i)
        temp[i] = A[i] * B[0] + C[i];

    } else if (A.cols() == B.rows() && A.rows() == C.rows() &&
        B.cols() == C.cols()) {
      // Case 3: m x n  *  n x k  +  m x n
      temp = Matrix<T> (A.rows(), B.cols(), false);
      
      for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < B.cols(); ++j) {
          temp[i * B.cols() + j] = C[i * B.cols() +j];
          for (int k = 0; k < B.rows(); ++k)
            temp[i * B.cols() + j] += A[i * A.cols() + k] * 
              B[k * B.cols() + j];
        }
      }

    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, std::string("Expects (m x n  *  1 x 1  +  m x n)") &
            "or (1 x 1  *  n x k  +  n x k) or (m x n  *  n x k  +" &
            "  m x k");
    }

    return temp;
  }

  /* Fast calculation of A'A */
  template <class T>
  Matrix<T>
  crossprod (const Matrix<T> &A)
  {
		int rows = A.rows();
		int cols = A.cols();
    Matrix<T> result(cols, cols, true);
		T tmp;
		
		if (rows == 1) {
			for (int k = 0; k < rows; ++k) {
				for (int i = 0; i < cols; ++i) {
					tmp = A[k * cols + i];
					for (int j = i; j < cols; ++j) {
						result[j * cols +i] =
							result[i * cols + j] += tmp * A[k * cols + j];
					}
				}
			}
		} else {
			for (int k = 0; k < rows; ++k) {
				for (int i = 0; i < cols; ++i) {
					tmp = A[k * cols + i];
					for (int j = i; j < cols; ++j) {
							result[i * cols + j] += tmp * A[k * cols + j];
					}
				}
			}

			for (int i = 0; i < cols; ++i)
				for (int j = i + 1; j < cols; ++j)
					result[j * cols + i] = result[i * cols + j];
		}
  
    return result;
  }

} // end namespace SCYTHE

#ifndef SCYTHE_COMPILE_DIRECT
#include "scythestat/eti/la.t"
#endif

#endif /* SCYTHE_LA_CC */
