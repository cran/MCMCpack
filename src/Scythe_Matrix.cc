/* Scythe_Matrix.cc
 *
 * This file provides the class implementation of the Matrix class, part
 * of the SCYTHE project.
 *
 * Scythe C++ Library
 * Copyright (C) Kevin M. Quinn, Andrew D. Martin,
 * and Daniel B. Pemstein
 *
 * This code written by:
 *
 * Kevin Quinn
 * Assistant Professor
 * Dept. of Political Science and
 * Center for Statistics and Social Sciences
 * Box 354322
 * University of Washington
 * Seattle, WA 98195-4322
 * quinn@stat.washington.edu
 *
 * Andrew D. Martin
 * Assistant Professor
 * Dept. of Political Science
 * Campus Box 1063
 * Washington University
 * St. Louis, MO 63130
 * admartin@artsci.wustl.edu
 * 
 * Daniel B. Pemstein
 * dbpemste@artsci.wustl.edu
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */

#ifndef SCYTHE_MATRIX_CC
#define SCYTHE_MATRIX_CC

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <new>
#include <numeric>
#include "Scythe_Error.h"
#include "Scythe_Matrix.h"

namespace SCYTHE {

  /**** Constructors ****/
	
  /* Default Constructor */
  template <class T>
  Matrix<T>::Matrix ()
    :	rows_ (0),
			cols_ (0),
			size_ (0),
			alloc_ (0),
			data_ (0)
  {
    data_ = new (std::nothrow) T[alloc_];
    if (data_ == 0) {
      throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
			       __LINE__, "Failure allocating null Matrix");
    }
  }

  /* Parameterized Type Constructor */
  template <class T>
  Matrix<T>::Matrix (const T &e)
    :	rows_ (1),
			cols_ (1),
			size_ (1),
			alloc_ (1),
			data_ (0)
  {
    data_ = new (std::nothrow) T[alloc_];
    if (data_ == 0) {
      throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
			       __LINE__, "Failure allocating Matrix of size 1");
    } else
      data_[0] = e;
  }

  /* Standard Constructor */
	template <class T>
	Matrix<T>::Matrix ( const int &n, const int &m, const bool &fill,
											const T &fill_value)
		:	rows_ (n),
			cols_ (m),
			size_ (n * m),
			alloc_ (1),
			data_ (0)
	{
		while (alloc_ < size_)
			alloc_ <<= 1;
    data_ = new (std::nothrow) T[alloc_];
    if (data_ == 0) {
      throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
			       __LINE__, "Failure allocating Matrix of size 1");
		}

		if (fill) {
			for (register int i = 0; i < alloc_; ++i)
				data_[i] = fill_value;
		}
	}

  /* Array Constructor */
  template <class T>
  Matrix<T>::Matrix (const int &n, const int &m, const T *in,
		     IN_TYPE type, const int &a, const int &b,
		     const int &c, const int &d)
    :	rows_ (n),
			cols_ (m),
			size_ (n * m),
			alloc_ (1),
			data_ (0)
  {
    /* This constructor is an interface weakness.  There is no easy
     * way to ensure that the array is of the length required.
     * Incorrect arrays may cause a seg fault.  Worse yet, if the
     * addresses reside on a page in this program's memory space,
     * corrupted data could enter the matrix.  On the other hand, this
     * constructor has far higher utility than any other.  We should
     * consider switching to a safe array type in the future.
     */
		while (alloc_ < size_)
			alloc_ <<= 1;
    data_ = new (std::nothrow) T[alloc_];
    if (data_ == 0) {
      throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__, __LINE__, 
			       std::string("Failure allocating Matrix of size ") & 
						 	(n * m));
    } else if (type == NORMAL) {
			for (register int i = 0; i < size_; ++i)
					data_[i] = in[i];
    } else if (type == REPEAT) {
      if (a <= 0 || a > n * m) {
				throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
						__LINE__, "REPEAT requires a s.t. 0 < a <= n * m ");
      } else {
				int cnt = -1;
				for (register int i = 0; i < size_; ++i) {
					if (cnt == a - 1)
						cnt = -1;
					data_[i] = in[++cnt];
				}
			}
    } else if (type == DIAG) {
      int cnt = -1;
      for (register int i = 0; i < rows_; ++i) {
				for (register int j = 0; j < cols_; ++j) {
					if (i == j)
						data_[i * cols_ + j] = in[++cnt];
					else
						data_[i * cols_ + j] = 0;
				}
      }
    } else if (type == UTRIANG) {
      int cnt = -1;
      for (register int i = 0; i < rows_; ++i) {
				for (register int j = 0; j < cols_; ++j) {
					if (i <= j)
						data_[i * cols_ + j] = in[++cnt];
					else
						data_[i * cols_ + j] = 0;
				}
      }
    } else if (type == LTRIANG) {
      int cnt = -1;
      for (register int i = 0; i < rows_; ++i) {
				for (register int j = 0; j < cols_; ++j) {
					if (i >= j) 
						data_[i * cols_ + j] = in[++cnt];
					else
						data_[i * cols_ + j] = 0;
				}
      }
    } else if (type == BLOCK) {
      if (a < 0 || b < 0 || c < a || d < b || c >= n || d >= m) {
				throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
				 	__LINE__,"BLOCK requires (a, b, c, d) s.t. 0 <= a <= c < n; \
					0 <= b <= d < m");
      } else {
				int cnt = -1;
				for (int i = 0; i < rows_; ++i) {
					for (int j = 0; j < cols_; ++j) {
						if (i >= a && i <= c && j >= b && j <= d)
							data_[i * cols_ + j] = in[++cnt];
						else
							data_[i * cols_ + j] = 0;
					}
				}
      }		
    } else { // undefined IN_TYPE 
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
			       __LINE__, "Undefined IN_TYPE");
    }
  }
	
  /* Load from file constructor */
  template <class T>
  Matrix<T>::Matrix(const std::string &path)
    :	rows_ (0),
			cols_ (0),
			size_ (0),
			alloc_(1),
			data_ (0)
  {
    std::ifstream file(path.c_str());
    if (! file) {
      throw scythe_file_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      std::string("Could not open ") & path);
    } else {
      file >> rows_ >> cols_;
			size_ = rows_ * cols_;
      if (file.eof() || rows_ <= 0 || cols_ <= 0) {
				throw scythe_file_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
						"Bad file format");
      } else {
				while (alloc_ < size_)
					alloc_ <<= 1;
				data_ = new (std::nothrow) T[alloc_];
				if (data_ == 0) {
					throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
							__LINE__,
							std::string("Failure allocating Matrix of size ")
							& size_);
				} else {
					for (int i = 0; i < size_; ++i) {
						if (file.eof())
							throw scythe_file_error(__FILE__, __PRETTY_FUNCTION__,
									__LINE__, std::string("Reached end of file before ")
									& size_ & " values were read");

							file >> data_[i];
						}
				}
      }
      file.close();
    }
  }

  /* Copy constructor */
  template <class T>
  Matrix<T>::Matrix (const Matrix<T> &m, const bool &fill)
    :	rows_ (m.rows_),
			cols_ (m.cols_),
			size_ (m.size_),
			alloc_ (m.alloc_),
			data_ (0)
  {
    data_ = new (std::nothrow) T[alloc_];
    if (data_ == 0) {
      throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      std::string("Failure allocating Matrix of size ") & size_);
    } else if (fill) {
      for (int i = 0; i < size_; ++i) {
				data_[i] = m.data_[i];
      }
    }
  }

  template <class T>
  template <class S>
  Matrix<T>::Matrix (const Matrix<S> &m)
    :	rows_ (m.rows()),
			cols_ (m.cols()),
			size_ (m.size()),
			alloc_ (1),
			data_ (0)
  {
		while (alloc_ < size_)
			alloc_ <<= 1;
    data_ = new (std::nothrow) T[alloc_];
    if (data_ == 0) {
      throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			    std::string("Failure allocating Matrix of size ") & size_);
    } else {
      S *mdata = m.getArray();
      for (int i = 0; i < size_; ++i) {
				data_[i] = (T) mdata[i];
      }
    }
  }

  /* Destructor */
  template <class T>
  Matrix<T>::~Matrix ()
  {
    delete[] data_;
  }

	/* swap two matrices' internals */
	template <class T>
	void
	Matrix<T>::swap (Matrix<T> &M)
	{
		int trows = rows_;
		int tcols = cols_;
		int tsize = size_;
		int talloc = alloc_;
		T *tdata = data_;

		rows_ = M.rows_;
		cols_ = M.cols_;
		size_ = M.size_;
		alloc_ = M.alloc_;
		data_ = M.data_;

		M.rows_ = trows;
		M.cols_ = tcols;
		M.size_ = tsize;
		M.alloc_ = talloc;
		M.data_ = tdata;
	}
	
  /* Submatrix operator */
  template <class T>
  Matrix<T>
  Matrix<T>::operator() (const int &a, const int &b, const int &c,
			 const int &d) const
  {
#ifndef SCYTHE_NO_RANGE
    if (c < a || d < b || !inRange(a,b) || !inRange(c,d)) {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			       "() requires (a, b, c, d) s.t. 0 <= a <= c < rows_; \
0 <= b <= d < cols");
      return Matrix((c - a + 1), (d - b + 1));
    }
#endif

    int cnt = -1;

		Matrix<T> temp((c - a + 1), (d - b + 1), false);
		for (int i = a; i <= c; ++i)
			for (int j = b; j <= d; ++j)
				temp.data_[++cnt] = data_[i * cols_ + j];
    
		return temp;
  }

  /* all rows extraction with _ struct */
  template <class T>
  Matrix<T>
  Matrix<T>::operator() (const all_elements& a, const int& j) const
	{
#ifndef SCYTHE_NO_RANGE
    if (j >= cols_ || j < 0) {
      throw scythe_out_of_range_error (__FILE__,__PRETTY_FUNCTION__,
				       __LINE__, std::string("Index ") & j &
				       " out of range");
    }
#endif
		//XXX
		Matrix<T> temp(rows_, 1, false);
		int k = j;
    for (register int i=0; i<rows_; ++i){
      temp.data_[i] = data_[k];
			k += cols_;
    }
    
    return temp;
  }
  
  /* all cols extraction with _ struct */
  template <class T>
  Matrix<T>
  Matrix<T>::operator() (const int& i, const all_elements& a) const
	{
#ifndef SCYTHE_NO_RANGE
    if (i >= rows_ || i < 0) {
      throw scythe_out_of_range_error (__FILE__,__PRETTY_FUNCTION__,
				       __LINE__, std::string("Index ") & i &
				       " out of range");
    }
#endif
    //XXX
    Matrix<T> temp(1, cols_, false);
		int k = i * cols_ - 1;
		for (register int j=0; j<cols_; ++j){
      temp.data_[j] = data_[++k];
    }
    
    return temp;
  }

	
  /**** Assignment and Arithmetic Assignment Operators ****/
	
  /* Assignment operator */
  template <class T>
  Matrix<T>
  &Matrix<T>::operator= (const Matrix<T> &m)
  {
    resize2Match(m);
    for (register int i = 0; i < size_; ++i)
      data_[i] = m.data_[i];
    return *this;
  }

  template <class T>
  template <class S>
  Matrix<T>
  &Matrix<T>::operator= (const Matrix<S> &m)
  {
    resize(m.rows(), m.cols(), false);
    S *mdata = m.getArray();
    for (register int i = 0; i < size_; ++i)
      data_[i] = (T) mdata[i];
    return *this;
  }

  /* Addition/assignment operator */
  template <class T>
  Matrix<T>
  &Matrix<T>::operator+= (const Matrix<T> &m)
  {
    if (size_ == 1) {
      // Case 1: 1X1 += nXm
      T temp = data_[0];
      resize2Match(m);
      for (int i = 0; i < size_; ++i) 
				data_[i] = temp + m.data_[i];
    } else if (m.size_ == 1) {
      // Case 2: nXm += 1X1
      for (int i = 0; i < size_; ++i)
				data_[i] += m.data_[0];
    } else if (rows_ == m.rows_ && cols_ == m.cols_) {
      // Case 3: nXm += nXm
      for (int i = 0; i < size_; ++i)
				data_[i] += m.data_[i];
    } else { // error
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices are not addition conformable");
    }
    return *this;
  }
	
  /* Subtraction/assignment operator */
  template <class T>
  Matrix<T>
  &Matrix<T>::operator-= (const Matrix<T> &m)
  {
    if (size_ == 1) {
      // Case 1: 1X1 -= nXm
      T temp = data_[0];
      resize2Match(m);
      for (int i = 0; i < size_; ++i)
	data_[i] = temp - m.data_[i];
    } else if (m.size_ == 1) {
      // Case 2: nXm -= 1X1
      for (int i = 0; i < size_; ++i)
	data_[i] -= m.data_[0];
    } else if (rows_ == m.rows_ && cols_ == m.cols_) {
      // Case 3: nXm -= nXm
      for (int i = 0; i < size_; ++i)
	data_[i] -= m.data_[i];
    } else { // error
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices are not subtraction conformable");
    }
    return *this;
  }

  /* Multiplication/assignment operator */
  template <class T>
  Matrix<T>
  &Matrix<T>::operator*= (const Matrix<T> &m)
  {
		if (size_ == 1) {
			// Case 1: 1X1 *= nXm
			T temp = data_[0];
			resize2Match(m);
			for (int i = 0; i  < size_; ++i)
				data_[i] = temp * m.data_[i];
		} else if (m.size_ == 1) {
			// Case 2: nXm *= 1X1
			for (int i = 0; i < size_; ++i)
				data_[i] *= m.data_[0];
		} else if (cols_ == m.rows_) {
			// Case 4: nXm *= mXk
			alloc_ = 1;
			while (alloc_ < size_)
				alloc_ <<= 1;
			T* temp = new (std::nothrow) T[alloc_];
			if (temp == 0) {
				throw scythe_alloc_error(__FILE__,__PRETTY_FUNCTION__, __LINE__,
						"Failure allocating space for multiplication");
				return *this;
			}
			for (register int i = 0; i < rows_; ++i) {
				for (register int j = 0; j < m.cols_; ++j) {
					temp[i * m.cols_ + j] = (T) 0;
					for (register int k = 0; k < m.rows_; ++k) {
						temp[i * m.cols_ + j] += data_[i * cols_ + k] *
							m.data_[k * m.cols_ + j];
					}
				}
			}
      /*
      const_col_major_iterator cmi = m.beginc();
      for (int i = 0; i < rows_; ++i) {
	for (int j = 0; j < m.cols_; ++j) {
	  temp[i * m.cols_ + j] = inner_product(&data_[i * cols_], 
						&data_[i * cols_ + m.rows_], cmi + (m.rows_ * j), (T) 0);
	}
      }
			*/
			
      cols_ = m.cols_;
			size_ = rows_;
			size_ *= cols_;
      delete[] data_;
      data_ = temp;
    } else { // error
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices are not multiplication conformable");
    }
    return *this;
  }
	
  /* Kronecker Multiplication/assignment operator */
  template <class T>
  Matrix<T>
  &Matrix<T>::operator%= (const Matrix<T> &m)
  {
    grow(size_ * m.size_);
    int cnt = size_ * m.size_  - 1;
    for (int i = rows_ - 1; i >= 0; --i) {
      for (int j = m.rows_ - 1; j >= 0; --j) {
				for (int k = cols_ - 1; k >= 0; --k) {
					for (int n = m.cols_ - 1; n >= 0; --n)
						data_[cnt--] = data_[i * cols_ + k] *
						m.data_[j * m.cols_ + n];
				}
      }
    }
    rows_ *= m.rows_;
    cols_ *= m.cols_;
		size_ = rows_ * cols_;
    return *this;
  }
			
  /* Division/assignment operator */
  template <class T>
  Matrix<T>
  &Matrix<T>::operator/= (const Matrix<T> &m)
  {
    if (size_ == 1) {
      T temp = data_[0];
      resize2Match(m);
      for (int i = 0; i < size_; ++i)
	data_[i] = temp / m.data_[i];
    } else if (m.size_ == 1) {
      for (int i = 0; i < size_; ++i)
	data_[i] /= m.data_[0];
    } else if (rows_ == m.rows_ && cols_ == m.cols_) {
      for (int i = 0; i < size_; ++i)
	data_[i] /= m.data_[i];
    } else { // error
      throw scythe_conformation_error (__FILE__, __PRETTY_FUNCTION__,
				       __LINE__, "Matrices are not division conformable");
    }
    return *this;
  }

  /* Power/assignment operator */
  template <class T>
  Matrix<T>
  &Matrix<T>::operator^= (const int &e)
  {
    if (e > 0) {
      if (! isSquare()) {
	throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
				     __LINE__,
				     std::string("Matrix must be square to raise it to  the ")
				     & e & " power");
	return *this;
      }
      Matrix<T> temp = *this;
      for (int i = 1; i < e; ++i)
				*(this) *= temp;
    } else if (e == 0) {
      // Case 3: A^0 == identity matrix of this size
      for (int i = 0; i < rows_; ++i) {
				for (int j = 0; j < cols_; ++j) {
					if (i == j)
						data_[ijIndex(i, j)] = 1.0;
					else
						data_[ijIndex(i, j)] = 0.0;
				}
      }
    } else if (e == -1) { 
      // Case 3: A^-1 == inverse of this
      if (! isSquare()) {
	throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
				     __LINE__, "Matrix not square");
	return *this;
      }
      if (isNull()) {
	throw scythe_null_error(__FILE__, __PRETTY_FUNCTION__,
				__LINE__);
	return *this;
      }
      Matrix<T> b(rows_, 1);
      Matrix<T> A = *(this);
      Matrix<T> L, U, perm_vec;

      // step 1: compute the LU factorization
      if (A.rows_ == 1) {
	L = Matrix<T>(1);
	U = *(this);
	perm_vec = Matrix<T>(1,1);
      } else {
	int pivot;
	L = U = Matrix<T>(A.rows_, A.rows_);
	perm_vec = Matrix<T>(A.rows_ - 1, 1);

	for (int k = 0; k < A.rows_ - 1; ++k) {
	  pivot = k;
	  // find pivot
	  for (int i = k; i < A.rows_; ++i) {
	    if (::fabs(A(pivot,k)) < ::fabs(A(i, k)))
	      pivot = i;
	  }
	  // check for singularity
	  if (A(pivot, k) == 0) {
	    throw scythe_type_error(__FILE__, __PRETTY_FUNCTION__,
				    __LINE__, "Matrix is singular");
	    return *this;
	  }
	  // permute
	  if (k != pivot) {
	    for (int i = 0; i < A.rows_; ++i) {
	      T temp = A(pivot, i);
	      A(pivot,i) = A(k,i);
	      A(k,i) = temp;
	    }
	  }
	  perm_vec[k] = pivot;

	  for (int i = k + 1; i < A.rows_; ++i) {
	    A(i,k) = A(i,k) / A(k,k);
	    for (int j = k + 1; j < A.rows_; ++j)
	      A(i,j) = A(i,j) - A(i,k) * A(k,j);
	  }
	}

	L = A;
	for (int i = 0; i < A.rows_; ++i) {
	  for (int j = i; j < A.rows_; ++j) {
	    U(i,j) = A(i,j);
	    L(i,j) = 0;
	    L(i,i) = 1;
	  }
	}
      }

      // step 2: repeated solving of A*hold = b
      for (int j = 0; j < A.rows_; ++j) {
	b[j] = 1;
				// Matrix hold = lu_solve(A, b, L, U, p);
				
				// step 2.1: solve L*y = Pb via forward substitution
				// Do a row interchange
	Matrix<T> bb = b;
	for (int ci = 0; ci < b.rows_ - 1; ++ci) {
	  int swap_row = static_cast<int>(perm_vec[ci]);
	  for (int cj = 0; cj < b.cols_; ++cj) {
	    T temp = bb(ci,cj);
	    bb(ci,cj) = b(swap_row,cj);
	    b(swap_row,cj) = temp;
	  }
	}
	Matrix<T> y(A.rows_, 1);
	for (int i = 0; i < A.rows_; ++i) {
	  T sum = 0;
	  for (int j = 0; j < i; ++j) {
	    sum += L(i,j) * y[j];
	  }
	  y[i] = (bb[i] - sum) / L(i,i);
	}

				// step 2.2: solve U*x = y via backsubstitution
	Matrix<T> x(A.rows_,1);
	for (int i = A.rows_ - 1; i >= 0; --i) {
	  T sum = 0;
	  for (int j = i + 1; j < A.rows_; ++j) {
	    sum += U(i,j) * x[j];
	  }
	  x[i] = (y[i] - sum) / U(i,i);
	}

				// step 3: reset b and put the solution into this
	b[j] = 0;
	for (int k = 0; k < A.rows_; ++k)
	  (*this)(k,j) = x[k];
      }
    } else { // error A^=n not defined where n < -1
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			       "Invalid argument: -1");
    }
    return *this;
  }

  /**** Accessors ****/
  template <class T>
  bool
  Matrix<T>::isDiagonal() const
  {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
	if (i != j && data_[ijIndex(i,j)] != 0)
	  return false;
      }
    }
    return true;
  }

  template <class T>
  bool
  Matrix<T>::isIdentity() const
  {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
	if (i != j) {
	  if (data_[ijIndex(i,j)] != 0)
	    return false;
	} else if (data_[ijIndex(i,j)] != 1)
	  return false;
      }
    }
    return true;
  }
	
  template <class T>
  bool
  Matrix<T>::isLowerTriangular() const
  {
    for (int i = 0; i < rows_; ++i) {
      for (int j = i + 1; j < cols_; ++j) {
	if (data_[ijIndex(i,j)] != 0)
	  return false;
      }
    }
    return true;
  }
	
  template <class T>
  bool
  Matrix<T>::isUpperTriangular() const
  {
    for (int j = 0; j < cols_; ++j) {
      for (int i = j + 1; i < rows_; ++i) {
	if (data_[ijIndex(i,j)] != 0)
	  return false;
      }
    }
    return true;
  }

  template <class T>
  bool
  Matrix<T>::isSingular() const
  {
    if (! isSquare() || isNull())
      return false;
    if ((~(*(this))) == (T) 0)
      return true;
    return false;
  }

  template <class T>
  bool
  Matrix<T>::isSymmetric() const
  {
    if (! isSquare())
      return false;
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
	if (data_[ijIndex(i,j)] != data_[ijIndex(j,i)])
	  return false;
      }
    }
    return true;
  }

  template <class T>
  bool
  Matrix<T>::isSkewSymmetric() const
  {
    if (! isSquare())
      return false;
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
	if (data_[ijIndex(i,j)] != 0 - data_[ijIndex(j,i)])
	  return false;
      }
    }
    return true;
  }

  /**** Utilities ****/

  template <class T>
  std::string
  Matrix<T>::toString (const unsigned int &prec, 
		       const unsigned int &width,
		       const bool &dim, const bool &internal) const
  {
    std::ostringstream s;
    unsigned int mlen = width;

    /* 2 Passes so we can get things to line up nicely */
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
	s.str("");
	s << std::setw(width) << std::setprecision(prec)
	  << std::setiosflags(std::ios::fixed)
	  << data_[ijIndex(i,j)];
	if (s.str().length() > mlen)
	  mlen = s.str().length();
      }
    }
		
    s.str("");
    if (dim) {
      s << "Size: " << size_ << " (" << rows_ << " x " << cols_
	<< ")" << std::endl;
    }
    if (internal) {
      s << "Object: " << this << ", Data: " << data_
	<< ", Allocation: " << alloc_ << std::endl;
    }

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
	s << std::setw(mlen) << std::setprecision(prec)
	  << data_[ijIndex(i,j)];
	if (i < rows_ - 1 || j < cols_ - 1)
	  s << " ";
      }
      s << std::endl;
    }
    return s.str();
  }

  template <class T>
  void
  Matrix<T>::save(const std::string &path, const char &flag,
		  const bool &header, const int &prec,
		  const int &width) const
  {
    std::ofstream out;
    bool err = false;
    if (flag == 'n') {
      std::fstream temp(path.c_str(), std::ios::in);
      if (!temp) {
	out.open(path.c_str(), std::ios::out);
      } else {
	temp.close();
	err = true;
      }
    }
    else if (flag == 'o')
      out.open(path.c_str(), std::ios::out | std::ios::trunc);
    else if (flag == 'a')
      out.open(path.c_str(), std::ios::out | std::ios::app);
    else {
      throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			       std::string("Incorrect flag ") & flag);
      return;
    }
    if (! out || err) {
      throw scythe_file_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      std::string("Could not open file ") & path);
      return;
    }
    if (header) {
      out << rows_ << " " << cols_;
      for (int i = 0; i < rows_ * cols_; ++i) {
	out << std::setw(width) << std::setprecision(prec) <<  " "
	    << data_[i];
      }
      out << std::endl;
    } else {
      out << toString(prec,width);
    }
    out << std::endl;
    out.close();
  }

  /**** Private Helper Methods ****/

  template <class T>
  bool
  operator== (const Matrix<T> &A, const Matrix<T> &B)
  {
    if (&A == &B)
      return true;
    if (A.rows() != B.rows() || A.cols() != B.cols())
      return false;
    for (int i = 0; i < A.size(); ++i) {
      if (A[i] != B[i])
	return false;
    }
    return true;
  }

  template <class T>
  bool
  operator== (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A == Matrix<T>(b));
  }
	
  template <class T>
  bool
  operator== (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) == B);
  }

  template <class T>
  bool
  operator!= (const Matrix<T> &A, const Matrix<T> &B)
  {
    return !(A == B);
  }
	
  template <class T>
  bool
  operator!= (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return !(A == Matrix<T>(b));
  }
	
  template <class T>
  bool
  operator!= (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return !(Matrix<T>(a) == B);
  }

  template <class T>
  bool
  operator< (const Matrix<T> &A, const Matrix<T> &B)
  {
    if (A.size() < B.size())
      return true;
    return false;
  }
	
  template <class T>
  bool
  operator< (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A < Matrix<T>(b));
  }
	
  template <class T>
  bool
  operator< (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) < B);
  }

  template <class T>
  bool
  operator> (const Matrix<T> &A, const Matrix<T> &B)
  {
    if (A.size() > B.size())
      return true;
    return false;
  }
	
  template <class T>
  bool
  operator> (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A > Matrix<T>(b));
  }
	
  template <class T>
  bool
  operator> (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) > B);
  }

  template <class T>
  bool
  operator<= (const Matrix<T> &A, const Matrix<T> &B)
  {
    return (A.size() == B.size() || A < B );
  }
	
  template <class T>
  bool
  operator<= (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A <= Matrix<T>(b));
  }
	
  template <class T>
  bool
  operator<= (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) <= B);
  }

  template <class T>
  bool
  operator>= (const Matrix<T> &A, const Matrix<T> &B)
  {
    return (A.size() == B.size() || A > B );
  }

  template <class T>
  bool
  operator>= (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A >= Matrix<T>(b));
  }
	
  template <class T>
  bool
  operator>= (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) >= B);
  }

  /* Element-by-Element inequality operators.  We define all three 
   * ((Matrix, Matrix), (Matrix, T), (T, Matrix)) although, in most
   * cases, the (T, Matrix) option will throw an error.  It would be
   * nice to rule this option out entirely and catch the problem at
   * compile-time but such comparisons are valid when the Matrix
   * passed in is 1x1
   */
  template <class T>
  Matrix<bool>
  operator|= (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<bool> C(A.rows(), A.cols(), false);
    if (A.isNull() || B.isNull()) {
      throw scythe_null_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      "Invalid use of NULL Matrix");
    } else if (B.isScalar()) {
      // Case 1: Compare every element in A to B[0]
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) == B[0]);
      }
    } else if (A.rows() == B.rows() && A.cols() == B.cols()) {
      // Case 2: equal size matrices
      for (int i = 0; i < A.size(); ++i)
	C[i] = (A[i] == B[i]);
    } else if (A.rows() == B.rows() && B.cols() == 1) {
      // Case 3: Matrix == col vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) == B[i]);
      }
    } else if (A.cols() == B.cols() && B.rows() == 1) {
      // Case 4: Matrix == row vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) == B[j]);
      }
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices not conformable");
    }
    return C;
  }

  template <class T>
  Matrix<bool>
  operator|= (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A |= Matrix<T>(b));
  }

  template <class T>
  Matrix<bool>
  operator|= (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) |= B);
  }
				
  template <class T>
  Matrix<bool>
  operator| (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<bool> C(A.rows(), A.cols(), false);
    if (A.isNull() || B.isNull()) {
      throw scythe_null_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      "Invalid use of NULL Matrix");
    } else if (B.isScalar()) {
      // Case 1: Compare every element in A to B[0]
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) != B[0]);
      }
    } else if (A.rows() == B.rows() && A.cols() == B.cols()) {
      // Case 2: equal size matrices
      for (int i = 0; i < A.size(); ++i)
	C[i] = (A[i] != B[i]);
    } else if (A.rows() == B.rows() && B.cols() == 1) {
      // Case 3: Matrix == col vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) != B[i]);
      }
    } else if (A.cols() == B.cols() && B.rows() == 1) {
      // Case 4: Matrix == row vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) != B[j]);
      }
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices not conformable");
    }
    return C;
  }
	
  template <class T>
  Matrix<bool>
  operator| (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A | Matrix<T>(b));
  }
	
  template <class T>
  Matrix<bool>
  operator| (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) | B);
  }
				
  template <class T>
  Matrix<bool>
  operator<< (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<bool> C(A.rows(), A.cols(), false);
    if (A.isNull() || B.isNull()) {
      throw scythe_null_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      "Invalid use of NULL Matrix");
    } else if (B.isScalar()) {
      // Case 1: Compare every element in A to B[0]
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) < B[0]);
      }
    } else if (A.rows() == B.rows() && A.cols() == B.cols()) {
      // Case 2: equal size matrices
      for (int i = 0; i < A.size(); ++i)
	C[i] = (A[i] < B[i]);
    } else if (A.rows() == B.rows() && B.cols() == 1) {
      // Case 3: Matrix == col vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) < B[i]);
      }
    } else if (A.cols() == B.cols() && B.rows() == 1) {
      // Case 4: Matrix == row vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) < B[j]);
      }
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices not conformable");
    }
    return C;
  }
	
  template <class T>
  Matrix<bool>
  operator<< (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A << Matrix<T>(b));
  }
	
  template <class T>
  Matrix<bool>
  operator<< (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) << B);
  }
				
  template <class T>
  Matrix<bool>
  operator>> (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<bool> C(A.rows(), A.cols(), false);
    if (A.isNull() || B.isNull()) {
      throw scythe_null_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      "Invalid use of NULL Matrix");
    } else if (B.isScalar()) {
      // Case 1: Compare every element in A to B[0]
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) > B[0]);
      }
    } else if (A.rows() == B.rows() && A.cols() == B.cols()) {
      // Case 2: equal size matrices
      for (int i = 0; i < A.size(); ++i)
	C[i] = (A[i] > B[i]);
    } else if (A.rows() == B.rows() && B.cols() == 1) {
      // Case 3: Matrix == col vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) > B[i]);
      }
    } else if (A.cols() == B.cols() && B.rows() == 1) {
      // Case 4: Matrix == row vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) > B[j]);
      }
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices not conformable");
    }
    return C;
  }
	
  template <class T>
  Matrix<bool>
  operator>> (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A >> Matrix<T>(b));
  }
	
  template <class T>
  Matrix<bool>
  operator>> (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) >> B);
  }
				
  template <class T>
  Matrix<bool>
  operator<<= (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<bool> C(A.rows(), A.cols(), false);
    if (A.isNull() || B.isNull()) {
      throw scythe_null_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      "Invalid use of NULL Matrix");
    } else if (B.isScalar()) {
      // Case 1: Compare every element in A to B[0]
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) <= B[0]);
      }
    } else if (A.rows() == B.rows() && A.cols() == B.cols()) {
      // Case 2: equal size matrices
      for (int i = 0; i < A.size(); ++i)
	C[i] = (A[i] <= B[i]);
    } else if (A.rows() == B.rows() && B.cols() == 1) {
      // Case 3: Matrix == col vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) <= B[i]);
      }
    } else if (A.cols() == B.cols() && B.rows() == 1) {
      // Case 4: Matrix == row vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) <= B[j]);
      }
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices not conformable");
    }
    return C;
  }
	
  template <class T>
  Matrix<bool>
  operator<<= (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A <<= Matrix<T>(b));
  }
	
  template <class T>
  Matrix<bool>
  operator<<= (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) <<= B);
  }
				
  template <class T>
  Matrix<bool>
  operator>>= (const Matrix<T> &A, const Matrix<T> &B)
  {
    Matrix<bool> C(A.rows(), A.cols());
    if (A.isNull() || B.isNull()) {
      throw scythe_null_error(__FILE__, __PRETTY_FUNCTION__, __LINE__,
			      "Invalid use of NULL Matrix");
    } else if (B.isScalar()) {
      // Case 1: Compare every element in A to B[0]
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) >= B[0]);
      }
    } else if (A.rows() == B.rows() && A.cols() == B.cols()) {
      // Case 2: equal size matrices
      for (int i = 0; i < A.size(); ++i)
	C[i] = (A[i] >= B[i]);
    } else if (A.rows() == B.rows() && B.cols() == 1) {
      // Case 3: Matrix == col vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) >= B[i]);
      }
    } else if (A.cols() == B.cols() && B.rows() == 1) {
      // Case 4: Matrix == row vector
      for (int i = 0; i < A.rows(); ++i) {
	for (int j = 0; j < A.cols(); ++j)
	  C(i,j) = (A(i,j) >= B[j]);
      }
    } else {
      throw scythe_conformation_error(__FILE__, __PRETTY_FUNCTION__,
				      __LINE__, "Matrices not conformable");
    }
    return C;
  }
	
  template <class T>
  Matrix<bool>
  operator>>= (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return (A >>= Matrix<T>(b));
  }
	
  template <class T>
  Matrix<bool>
  operator>>= (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return (Matrix<T>(a) >>= B);
  }
	
  /**** Arithmetic Operators ****/
  /* When operations are commutative we pass both args in by reference
   * and check size before copying for (Matrix, Matrix) and always
   * copy the T where a T is involved.  For non-commutative, we pass
   * the first arg by value.  When possible, the user should try to
   * pass the smaller Matrix first although subsequent resizes may
   * make up for the time saved in copying in some cases.
   */
  template <class T>
  Matrix<T>
  operator+ (const Matrix<T> &A, const Matrix<T> &B)
  {
    // If A or B is 1x1 this can save some time
    if (A.size() < B.size())
      return Matrix<T>(A) += B;
    return Matrix<T>(B) += A;
  }

  template <class T>
  Matrix<T>
  operator+ (const Matrix<T> &A, const typename Matrix<T>::ttype &b)
  {
    return Matrix<T>(b) += A;
  }
	
  template <class T>
  Matrix<T>
  operator+ (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return Matrix<T>(a) += B;
  }
	
  template <class T>
  Matrix<T>
  operator- (Matrix<T> A, const Matrix<T> &B)
  {
    return A -= B;
  }
	
  template <class T>
  Matrix<T>
  operator- (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A -= Matrix<T>(b);
  }
	
  template <class T>
  Matrix<T>
  operator- (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return Matrix<T>(a) -= B;
  }
	
  // No overload here cause ambiguous with standard negation
  template <class T>
  Matrix<T>
  operator-(Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = -A[i];
    return A;
  }

  template <class T>
  Matrix<T>
  operator* (Matrix<T> A, const Matrix<T> &B)
  {
    return A *= B;
  }

  // 1X1 * nXm is commutative but we'd grow anyway so don't bother
  template <class T>
  Matrix<T>
  operator* (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A *= Matrix<T>(b);
  }

  template <class T>
  Matrix<T>
  operator* (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return Matrix<T>(a) *= B;
  }
	
  template <class T>
  Matrix<T>
  operator% (Matrix<T> A, const Matrix<T> &B)
  {
    return A %= B;
  }

  // commutative but we'd grow anyway
  template <class T>
  Matrix<T>
  operator% (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A %= Matrix<T>(b);
  }

  template <class T>
  Matrix<T>
  operator% (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return Matrix<T>(a) %= B;
  }
	
  template <class T>
  Matrix<T>
  operator/ (Matrix<T> A, const Matrix<T> &B)
  {
    return A /= B;
  }
	
  template <class T>
  Matrix<T>
  operator/ (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A /= Matrix<T>(b);
  }

  template <class T>
  Matrix<T>
  operator/ (const typename Matrix<T>::ttype &a, const Matrix<T> &B)
  {
    return Matrix<T>(a) /= B;
  }
	
  template <class T>
  Matrix<T>
  operator^ (Matrix<T> A, const int &i)
  {
    return A ^= i;
  }

  template <class T>
  Matrix<T>
  operator! (const Matrix <T> &m)
  {
		int rows = m.rows();
		int cols = m.cols();
    Matrix<T> ret(cols, rows, false);
    for (register int i = 0; i < rows; ++i) {
      for (register int j = 0; j < cols; ++j)
				//ret(j,i) = m(i,j);
				ret[j * rows + i] = m[i * cols + j];
    }
    return ret;
  }

  template <class T>
  T	
  operator~ (Matrix <T> A) // no reference because LU kills the matrix
  {
    if (! A.isSquare()) {
      throw scythe_dimension_error(__FILE__,__PRETTY_FUNCTION__,
				   __LINE__, "Matrix not square");
      return 0;
    }
    if (A.isNull()) {
      throw scythe_null_error(__FILE__,__PRETTY_FUNCTION__,
			      __LINE__, "Matrix is NULL");
      return 0;
    }
    Matrix<T> L(A.rows(), A.rows());
    Matrix<T> U = L;
    T sign = 1;
    int pivot;
    T temp;

    for (int k = 0; k < A.rows(); ++k) {
      pivot = k;
      // find pivot
      for (int i = k; i < A.rows(); ++i)
	if (A(pivot,k) < ::fabs(A(i,k)))
	  pivot = i;
			
      if (A(pivot,k) == 0)
	return 0;

      // permute
      if (k != pivot) {
	sign *= -1;
	for (int i = k; i < A.rows(); ++i) {
	  temp = A(pivot,i);
	  A(pivot,i) = A(k,i);
	  A(k,i) = temp;
	}
      }

      for (int i = k + 1; i < A.rows(); ++i) {
	A(i,k) = A(i,k) / A(k,k);
	for (int j = k + 1; j < A.rows(); ++j) {
	  A(i,j) = A(i,j) - A(i,k) * A(k,j);
	}
      }
    }

    T det = 1;
    for (int i = 0; i < A.rows(); ++i)
      det *= A(i,i);

    return sign * det;
  }

}	// end namespace SCYTHE

#endif /* SCYTHE_MATRIX_CC */
