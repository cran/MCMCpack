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
 * scythestat/ide.cc
 *
 * Provides implementations for inversion and decomposition
 * template functions that operate on Scythe's Matrix class.
 *
 */

#ifndef SCYTHE_IDE_CC
#define SCYTHE_IDE_CC

#include <cmath>
#include <algorithm>

#ifdef SCYTHE_COMPILE_DIRECT
#include "ide.h"
#include "error.h"
#include "util.h"
#else
#include "scythestat/ide.h"
#include "scythestat/error.h"
#include "scythestat/util.h"
#endif

namespace SCYTHE {

  /* Cholesky Decomposition of a Symmetric Positive Definite Matrix */
  template <class T>
  Matrix<T>
  cholesky (const Matrix<T> &A){
  
    if (! A.isSquare()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrix not square");
    }
    
    Matrix<T> temp (A.rows(), A.cols(), false);
    register T h;
    
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = i; j < A.cols(); ++j) {
        h = A(i,j);
        for (int k = 0; k < i; ++k) {
          h -= temp(i, k) * temp(j, k);
        }
        if (i == j) {
          if (h <= (T) 0) {
            throw scythe_type_error(__FILE__, __PRETTY_FUNCTION__,
                __LINE__, "Matrix not positive definite");
          }
          temp(i,i) = std::sqrt(h);
        } else {
          temp(j,i) = (((T) 1) / temp(i,i)) * h;
          temp(i,j) = (T) 0;
        }
      }
    }
  
    return temp;
  }

  /* Solve Ax=b for x via backsubstitution using cholesky decomp */
  template <class T>
  Matrix<T>
  chol_solve (const Matrix<T> & A, const Matrix<T> & b)
  {
    /* NOTE: cholesky() call does check for square/posdef of A */
  
    if ((! b.isColVector()) || A.rows() != b.rows() || ! A.isSquare()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Inputs not proper dimension");
    }
  
    Matrix<T> M = cholesky (A);
    register T holder;
    register T *y = new T[A.rows()];
    register T *x = new T[A.rows()];
     
    // solve M*y = b
    for (int i = 0; i < A.rows(); ++i) {
      holder = (T) 0;
      for (int j = 0; j < i; ++j) {
        holder += M(i,j) * y[j];
      }
      y[i] = (((T) 1) / M(i,i)) * (b[i] - holder);
    }
  
    // solve M'*x = y
    for (int i = A.rows() - 1; i >= 0; --i) {
      holder = (T) 0;
      for (int j = i + 1; j < A.rows(); ++j) {
        holder += M(j,i) * x[j];
      }
      x[i] = (((T) 1) / M(i,i)) * (y[i] - holder);
    }    
      
    Matrix<T> temp (A.rows(), 1, x);
    delete[]y;
    delete[]x;
  
    return temp;
  }
  
  /* Solve Ax=b for x via backsub using cholesky decomp */
  template <class T>
  Matrix<T>
  chol_solve (const Matrix<T> &A, const Matrix<T> &b,
              const Matrix<T> &M)
  {
    if (b.cols() != 1 || A.rows() != b.rows() || A.rows() != M.rows()
      || ! A.isSquare() || ! M.isSquare()) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Inputs not proper dimension");
    }
    register T *y = new T[A.rows()];
    register T *x = new T[A.rows()];
    register T holder;
    
    // solve M*y = b
    for (int i = 0; i < A.rows(); ++i) {
      holder = 0.0;
      for (int j = 0; j < i; ++j) {
        holder += M(i,j) * y[j];
      }
      y[i] = (1.0 / M(i,i)) * (b[i] - holder);
    }
        
    // solve M'*x = y
    for (int i = A.rows() - 1; i >= 0; --i) {
      holder = 0.0;
      for (int j = i + 1; j < A.rows(); ++j) {
        holder += M(j,i) * x[j];
      }
      x[i] = (1.0 / M(i,i)) * (y[i] - holder);
    }
      
    Matrix<T> temp (A.rows(), 1, x);
    delete[]y;
    delete[]x;
    
    return temp;
  }


  /* Calculate the inverse of a symmetric positive definite matrix */
  template <class T>
  Matrix<T>
  invpd (const Matrix<T> &A)
  {  
    // SYMMETRY OF A *IS NOT* CHECKED
    if (! A.isSquare())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrix not square");
  
    // Cholesky decomp
    Matrix<T> M (A.rows(), A.cols(), false);
    register T h;
    
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = i; j < A.cols(); ++j) {
        h = A(i,j);
        for (int k = 0; k < i; ++k) {
          h -= M(i, k) * M(j, k);
        }
        if (i == j) {
          if (h <= (T) 0) {
            throw scythe_type_error(__FILE__, __PRETTY_FUNCTION__,
                __LINE__, "Matrix not positive definite");
          }
          M(i,i) = std::sqrt(h);
        } else {
          M(j,i) = (((T) 1) / M(i,i)) * h;
          M(i,j) = (T) 0;
        }
      }
    }
  
    // for chol_solve block
    register T *y = new T[A.rows()];
    register T *x = new T[A.rows()];
    Matrix<T> b(A.rows(), 1); // full of zeros
    
    // For final answer
    Matrix<T> Ainv(A.rows(), A.cols(), false);

    for (int k = 0; k < A.rows(); ++k) {
      b[k] = (T) 1;

      // begin chol_solve
      // solve M*y = b
      for (int i = 0; i < A.rows(); ++i) {
        h = (T) 0;
        for (int j = 0; j < i; ++j) {
          h += M(i,j) * y[j];
        }
        y[i] = (((T) 1) / M(i,i)) * (b[i] - h);
      }
    
      // solve M'*x = y
      for (int i = A.rows() - 1; i >= 0; --i) {
        h = (T) 0;
        for (int j = i + 1; j < A.rows(); ++j) {
          h += M(j,i) * x[j];
        }
        x[i] = (((T) 1) / M(i,i)) * (y[i] - h);
      }    
      // end chol_solve

      b[k] = (T) 0;
      for (int l = 0; l < A.rows(); ++l)
        Ainv(l,k) = x[l];
    }

    delete[] y;
    delete[] x;

    return Ainv;
  }


 /* Calculates the inverse of a Symmetric Positive Definite Matrix  */
  template <class T>
  Matrix<T>
  invpd (const Matrix<T> &A, const Matrix<T> &M)
  {
    if (A.rows() != M.cols() || A.cols() != M.rows())
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "A and M do not conform");
      
    register T h;

    // for chol_solve block
    register T *y = new T[A.rows()];
    register T *x = new T[A.rows()];
    Matrix<T> b(A.rows(), 1); // full of zeros
    
    // For final answer
    Matrix<T> Ainv(A.rows(), A.cols(), false);

    for (int k = 0; k < A.rows(); ++k) {
      b[k] = (T) 1;

      // begin chol_solve
      // solve M*y = b
      for (int i = 0; i < A.rows(); ++i) {
        h = (T) 0;
        for (int j = 0; j < i; ++j) {
          h += M(i,j) * y[j];
        }
        y[i] = (((T) 1) / M(i,i)) * (b[i] - h);
      }
    
      // solve M'*x = y
      for (int i = A.rows() - 1; i >= 0; --i) {
        h = (T) 0;
        for (int j = i + 1; j < A.rows(); ++j) {
          h += M(j,i) * x[j];
        }
        x[i] = (((T) 1) / M(i,i)) * (y[i] - h);
      }    
      // end chol_solve

      b[k] = (T) 0;
      for (int l = 0; l < A.rows(); ++l)
        Ainv(l,k) = x[l];
    }

    delete[] y;
    delete[] x;

    return Ainv;
  }

//  This code is based on  Algorithm 3.4.1 of Golub and Van Loan 
//  3rd edition, 1996. Major difference is in how the output is 
//  structured. 

   /* Calculates the LU Decomposition of a square Matrix */
  template <class T>
  void
  lu_decomp(Matrix<T> A, Matrix<T> &L, Matrix<T> &U,
            Matrix<int> &perm_vec)
  {
    if (! A.isSquare())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrix A not square");

    if (A.isRowVector()) {
      L = Matrix<T> (1, 1, true, 1); // all 1s
      U = A;
      perm_vec = Matrix<int>(1, 1);  // all 0s
      return;
    }
    
    L = U = Matrix<T>(A.rows(), A.cols(), false);
    perm_vec = Matrix<int> (A.rows() - 1, 1, false);

    int pivot;
    T temp;

    for (int k = 0; k < A.rows() - 1; ++k) {
      pivot = k;
      // find pivot
      for (int i = k; i < A.rows(); ++i) {
        if (std::fabs(A(pivot,k)) < std::fabs(A(i,k)))
          pivot = i;
      }
      
      if (A(pivot,k) == (T) 0)
        throw scythe_type_error(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Matrix is singular");

      // permute
      if (k != pivot) {
        for (int i = 0; i < A.rows(); ++i) {
          temp = A(pivot,i);
          A(pivot,i) = A(k,i);
          A(k,i) = temp;
        }
      }
      perm_vec[k] = pivot;

      for (int i = k + 1; i < A.rows(); ++i) {
        A(i,k) = A(i,k) / A(k,k);
        for (int j = k + 1; j < A.rows(); ++j)
          A(i,j) = A(i,j) - A(i,k) * A(k,j);
      }
    }

    L = A;

    for (int i = 0; i < A.rows(); ++i) {
      for (int j = i; j < A.rows(); ++j) {
        U(i,j) = A(i,j);
        L(i,j) = (T) 0;
        L(i,i) = (T) 1;
      }
    }
  }


  /* Solves A*x=b for x via lu_decomp */
  template <class T>
  Matrix<T>
  lu_solve(Matrix<T> A, const Matrix<T> &b)
  {
    if (! b.isColVector())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "b is not a column vector");
    
    if (! A.isSquare())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrix A not square");

    if (A.rows() != b.rows())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "A.rows() != b.rows()");
    
    // step 1 compute the LU factorization 
    Matrix<T> L, U;
    Matrix<int> perm_vec;
    
    if (A.isRowVector()) {
      L = Matrix<T> (1, 1, true, 1); // all 1s
      U = A;
      perm_vec = Matrix<int>(1, 1);  // all 0s
    } else {
    
      L = U = Matrix<T>(A.rows(), A.cols(), false);
      perm_vec = Matrix<int> (A.rows() - 1, 1, false);
  
      int pivot;
      T temp;
  
      for (int k = 0; k < A.rows() - 1; ++k) {
        pivot = k;
        // find pivot
        for (int i = k; i < A.rows(); ++i) {
          if (std::fabs(A(pivot,k)) < std::fabs(A(i,k)))
            pivot = i;
        }
        
        if (A(pivot,k) == (T) 0)
          throw scythe_type_error(__FILE__, __PRETTY_FUNCTION__,
              __LINE__, "Matrix is singular");
  
        // permute
        if (k != pivot) {
          for (int i = 0; i < A.rows(); ++i) {
            temp = A(pivot,i);
            A(pivot,i) = A(k,i);
            A(k,i) = temp;
          }
        }
        perm_vec[k] = pivot;
  
        for (int i = k + 1; i < A.rows(); ++i) {
          A(i,k) = A(i,k) / A(k,k);
          for (int j = k + 1; j < A.rows(); ++j)
            A(i,j) = A(i,j) - A(i,k) * A(k,j);
        }
      }
  
      L = A;
    
      for (int i = 0; i < A.rows(); ++i) {
        for (int j = i; j < A.rows(); ++j) {
          U(i,j) = A(i,j);
          L(i,j) = (T) 0;
          L(i,i) = (T) 1;
        }
      }
    }
    // step 2 solve L*y = Pb via forward substitution
    Matrix<T> bb = row_interchange(b, perm_vec);
    Matrix<T> y(A.rows(), 1, false);
    T sum;
  
    for (int i = 0; i < A.rows(); ++i) {
      sum = (T) 0;
      for (int j = 0; j < i; ++j)
        sum += L[i * A.cols() + j] *  y[j];
  
      y[i] = (bb[i] - sum) / L[i * A.cols() + i];
    }
  
    // step 3 solve U*x = y via backsubstitution
    Matrix<T> x(A.rows(), 1, false);
    for (int i = A.rows() - 1; i >= 0; --i) {
      sum = (T) 0;
      for (int j = i + 1; j < A.rows(); ++j)
        sum += U[i * A.cols() + j] *  x[j];
      
      x[i] = (y[i] - sum) / U[i * A.cols() + i];
    }
    
    return x;
  }

  /* lu_solve overrloaded: you need A, b + L, U, perm_vec from
   * lu_decomp
   */
  template <class T>
  Matrix<T>
  lu_solve (Matrix<T> A, const Matrix<T> &b, const Matrix<T> &L,
            const Matrix<T> &U, const Matrix<int> &perm_vec) 
  {
    if (! b.isColVector())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "b is not a column vector");

    if (! A.isSquare())
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "A is not square");

    if (A.rows() != b.rows())
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "A and b have different row sizes");

    if (A.rows() != L.rows() || A.rows() != U.rows() ||
        A.cols() != L.cols() || A.cols() != U.cols())
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "A, L, and U do not conform");

    if (perm_vec.rows() + 1 != A.rows())
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
        __LINE__, "perm_vec does not have exactly one less row than A");
  
  
    // step 1 solve L*y = Pb via forward substitution
    Matrix<T> bb = row_interchange(b, perm_vec);
    Matrix<T> y(A.rows(), 1, false);
    T sum;
  
    for (int i = 0; i < A.rows(); ++i) {
      sum = (T) 0;
      for (int j = 0; j < i; ++j)
        sum += L[i * A.cols() + j] *  y[j];
  
      y[i] = (bb[i] - sum) / L[i * A.cols() + i];
    }
  
    // step 2 solve U*x = y via backsubstitution
    Matrix<T> x(A.rows(), 1, false);
    for (int i = A.rows() - 1; i >= 0; --i) {
      sum = (T) 0;
      for (int j = i + 1; j < A.rows(); ++j)
        sum += U[i * A.cols() + j] *  x[j];
      
      x[i] = (y[i] - sum) / U[i * A.cols() + i];
    }
    
    return x;
  }

  /* Interchanges the rows of A with those in vector p */
	//XXX maybe I should inline this and get rid of .t lines
  template <class T>
  Matrix<T>
  row_interchange(Matrix<T> A, const Matrix<int> &p){
    if (! p.isColVector()) {
      throw scythe_dimension_error (__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "p not a column vector");
    }
    if (p.rows() + 1 != A.rows()) {
      throw scythe_conformation_error (__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "p must have one less row than A");
    }

    for (int i = 0; i < A.rows() - 1; ++i)
      swap_ranges(A.vec(i), A.vec(i + 1), A.vec(p[i]));
    
    return A;
  }
  
}  // end namespace SCYTHE

#ifndef SCYTHE_COMPILE_DIRECT
#include "scythestat/eti/ide.t"
#endif
  
#endif /* SCYTHE_IDE_CC */
