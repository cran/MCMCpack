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
 * scythestat/matrix.h
 *
 * Provides the class definition for the Matrix class; this
 * data structure sits at the core of the library.  This class
 * behaves according to the Standard Template Library (STL)
 * standard for container classes and iterators for this class
 * are provided by include/Scythe_Matrix_Iterator.h
 *
 */

#ifndef SCYTHE_MATRIX_H
#define SCYTHE_MATRIX_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <new>
#include <numeric>
#include <string>
#include <climits>
#include <cmath>

#ifdef SCYTHE_COMPILE_DIRECT
#include "error.h"
#include "util.h"
#include "matrix_iterator.h"
#else
#include "scythestat/error.h"
#include "scythestat/util.h"
#include "scythestat/matrix_iterator.h"
#endif

namespace SCYTHE {
  
  struct all_elements{
  }  const _ = {};
 
  
  enum IN_TYPE {NORMAL, REPEAT, DIAG, UTRIANG, LTRIANG, BLOCK};

  template <class T>
  class Matrix
  {
  public:
    typedef T ttype;
		
		friend class matrix_iterator<ttype>;
    friend class const_matrix_iterator<ttype>;
    friend class row_major_iterator<ttype>;
    friend class const_row_major_iterator<ttype>;
    friend class col_major_iterator<ttype>;
    friend class const_col_major_iterator<ttype>;
    friend class reverse_row_major_iterator<ttype>;
    friend class const_reverse_row_major_iterator<ttype>;
    friend class reverse_col_major_iterator<ttype>;
    friend class const_reverse_col_major_iterator<ttype>;
    
		typedef matrix_iterator<ttype> iterator;
    typedef const_matrix_iterator<ttype> const_iterator;
    typedef row_major_iterator<ttype> row_major_iterator;
    typedef const_row_major_iterator<ttype> const_row_major_iterator;
    typedef col_major_iterator<ttype> col_major_iterator;
    typedef const_col_major_iterator<ttype> const_col_major_iterator;
    typedef reverse_row_major_iterator<ttype>
      reverse_row_major_iterator;
    typedef const_reverse_row_major_iterator<ttype>
      const_reverse_row_major_iterator;
    typedef reverse_col_major_iterator<ttype>
      reverse_col_major_iterator;
    typedef const_reverse_col_major_iterator<ttype>
      const_reverse_col_major_iterator;

    /**** Constructors ****/

    /* Default Constructor: Creates a Matrix of size 0.  This
     * Matrix cannot be used in operations but is useful when you
     * want to make arrays of matrices.
     */
    Matrix ()
			:  rows_ (0),
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

    /* Parameterized Type Constructor:  Create a 1x1 Matrix
     * containing one value.  While this is really just a scalar,
     * it has its uses:  Necessary for assignments such as Matrix A
     * = 3; or Matrix B = A[0];.  Also means we only have to define
     * operators for Matrix objects.
     */
    Matrix (const T& e) 
			:  rows_ (1),
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

    /* Standard Constructor:  Creates an n x m Matrix.  By default
     * it fills the Matrix with zeroes but you can turn this off or
     * change the fill value
     */
    explicit
    Matrix (const int &n, const int &m, const bool &fill = true,
            const T &fill_value = 0)
			:  rows_ (n),
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

    /* Array Constructor:  Creates an n x m Matrix from an array of
     * type T.  You may enter an IN_TYPE as enumerated above.  If
     * NORMAL (the default), DIAG, UTRIANG, LTRIANG, the
     * a,b,c,d should not be entered.  Doing so is not an error but
     * the values will be ignored.  The array should simply be the
     * correct length.  NORMAL requires n*m elements, DIAG n,
     * TRIDIAG 3n - 3, UTRIANG and LTRIANG require enough elements
     * to fill the upper and lower triangles of the matrix
     * respectively.  If you choose UTRIANG or LTRIANG for a 1 x 1
     * matrix you should have 1 element(although this isn't
     * technically a triangular matrix), 1 x m and n x 1 vectors
     * produce equally strange but logical results. The equations
     * for required array sizes for UTRIANG and LTRIANG are,
     * respectively: (n * m) - ((m * (m - 1)) / 2) - ((max(0, m - n)
     * * (m - n - 1)) / 2 ) and (n * m) - (n * (n - 1)) / 2) -
     * ((max(0, n - m) * (n - m - 1)) / 2). (Nick you could try to
     * simplify these equations for the documentation if possible).
     * REPEAT takes a small array and repeats it throughout the
     * matrix with a defined as the length of the array s.t. 0 < a
     * <= n * m.  We don't require that m * n is divisible by a but
     * one would expect this in most cases.  BLOCK places a block
     * array in a matrix of 0s s.t (a, b), (c, d) deliminate the
     * corners of the block;  0 <= a <= c < n; 0 <= b <= d < m.  The
     * array must be of size (c - a + 1) * (d - b + 1).
     */
    Matrix (const int &n, const int &m, const T* in,
            IN_TYPE type = NORMAL, const int &a = -1,
            const int &b = -1, const int &c = -1,
            const int &d = -1)
			:  rows_ (n),
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
						 __LINE__,
						 "BLOCK requires (a, b, c, d) s.t. 0 <= a <= c < n; \
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
    
    /* Create a matrix from a file: first two elements should be
     * # row and # col.  Followed by a list of elements.  Limited
     * range-checking is done: if the constructor reaches eof before
     * extracting the expected number of elements, an error will be
     * thrown.  Row/col <= 0 will be caught.  If you
     * forget the row/col it will consider the first 2 numbers as
     * row and col.  If your file is longer than expected, no error
     * will be thrown.
     */
    Matrix (const std::string &path)
			:  rows_ (0),
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
					throw scythe_file_error(__FILE__, __PRETTY_FUNCTION__,
							__LINE__, "Bad file format");
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

    /* Copy Constructor: Create a copy of an existing Matrix */
    Matrix (const Matrix<T> &m, const bool &fill=true)
			:  rows_ (m.rows_),
				cols_ (m.cols_),
				size_ (m.size_),
				alloc_ (m.alloc_),
				data_ (0)
		{
			data_ = new (std::nothrow) T[alloc_];
			if (data_ == 0) {
				throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
						__LINE__,
						std::string("Failure allocating Matrix of size ") & size_);
			} else if (fill) {
				for (int i = 0; i < size_; ++i) {
					data_[i] = m.data_[i];
				}
			}
		}

    template <class S>
    Matrix (const Matrix<S> &m)
			:  rows_ (m.rows()),
				cols_ (m.cols()),
				size_ (m.size()),
				alloc_ (1),
				data_ (0)
		{
			while (alloc_ < size_)
				alloc_ <<= 1;
			data_ = new (std::nothrow) T[alloc_];
			if (data_ == 0) {
				throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
						__LINE__,
						std::string("Failure allocating Matrix of size ") & size_);
			} else {
				S *mdata = m.getArray();
				for (int i = 0; i < size_; ++i) {
					data_[i] = (T) mdata[i];
				}
			}
		}

    /**** Destructor ****/
    ~Matrix ()
		{
			delete[] data_;
		}

    /**** STL container modifiers ****/
    /* Swap operator (sort of a dual copy constructor) */
    void swap (Matrix<T> &M)
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

    inline void clear ()
    {
      resize(0, 0);
    }

    /**** Indexing Operators ****/

    /* Retrieve the ith element in row major order */
    inline T &operator[] (const int &i)
    {
#ifndef SCYTHE_NO_RANGE
      if ( !inRange(i)) {
        throw scythe_out_of_range_error (__FILE__,__PRETTY_FUNCTION__,
             __LINE__, std::string("Index ") & i & 
             " out of range");
      }
#endif
      return data_[i];
    }

    /* Retrieve the (i,j)th element */
    inline T &operator() (const int &i, const int &j)
    {
#ifndef SCYTHE_NO_RANGE
      if (! inRange(i, j)) {
        throw scythe_out_of_range_error(__FILE__,__PRETTY_FUNCTION__,
            __LINE__, std::string("Index (") & i & "," &  j & 
            ") out of range");
      }
#endif
      return data_[i * cols_ + j];
    }

    /* Versions of the above two for const Matrix objects */
    inline T &operator[] (const int &i) const
    {
#ifndef SCYTHE_NO_RANGE
      if (! inRange(i)) {
        throw scythe_out_of_range_error (__FILE__,__PRETTY_FUNCTION__,
             __LINE__, std::string("Index ") & i &
             " out of range");
      }
#endif
      return data_[i];
    }
    
    inline T &operator() (const int &i, const int &j) const
    {
#ifndef SCYTHE_NO_RANGE
      if (! inRange(i, j)) {
        throw scythe_out_of_range_error(__FILE__,__PRETTY_FUNCTION__,
            __LINE__, std::string("Index (") & i & "," &  j & 
            ") out of range");
      }
#endif
      return data_[i * cols_ + j];
    }

    /* SubMatrix operator: returns a new Matrix with a,b,c,d
     * defining the bounds of the block s.t 0 <= a <= c < rows_;
     * 0 <= b <= d < cols_.
     */
    Matrix<T> operator() (const int &a, const int &b, const int &c,
                          const int &d) const
		{
#ifndef SCYTHE_NO_RANGE
			if (c < a || d < b || !inRange(a,b) || !inRange(c,d)) {
				throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
						__LINE__, "() requires (a, b, c, d) s.t. 0 <= a <= c < \
						rows_; 0 <= b <= d < cols");
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

    
    /* function for extracting all rows with the _ struct
     */
    Matrix<T> operator() (const all_elements& a, const int& j) const
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

    /* function for extracting all columns with the _ struct
     */
    Matrix<T> operator() (const int& i, const all_elements& a) const
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

    /**** Self-modifying arithmetic operators ****/

    /* Assignment operator */
    Matrix<T> &operator= (const Matrix<T> &m)
		{
			resize2Match(m);
			for (register int i = 0; i < size_; ++i)
				data_[i] = m.data_[i];
			return *this;
		}

    template <class S>
    Matrix<T> &operator= (const Matrix<S> &m)
		{
			resize(m.rows(), m.cols(), false);
			S *mdata = m.getArray();
			for (register int i = 0; i < size_; ++i)
				data_[i] = (T) mdata[i];
			return *this;
		}

    /* Matrix addition/assignment */
    Matrix<T> &operator+= (const Matrix<T> &m)
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

    /* Matrix subtraction/assignment */
    Matrix<T> &operator-= (const Matrix<T> &m)
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

    /* Matrix multiplication/assignment */
    Matrix<T> &operator*= (const Matrix<T> &m)
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
				while (alloc_ < rows_ * m.cols_)
					alloc_ <<= 1;
				T* temp = new (std::nothrow) T[alloc_];
				if (temp == 0) {
					throw scythe_alloc_error(__FILE__,__PRETTY_FUNCTION__,
							__LINE__, "Failure allocating space for multiplication");
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

    /* Kronecker Multiplication/assignment */
    Matrix<T> &operator%= (const Matrix<T> &m)
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

    /* Element-by-element division/assignment */
    Matrix<T> &operator/= (const Matrix<T> &m)
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

    /* Matrix power/assignment: ^ positive integer does matrix
     * power, ^0 returns an identity matrix of the input's size, ^-1
     * returns inverse via LU Decomposition.  Input must be
     * square for ^posint and ^-1.  You should use invpd in
     * Scythe_IDE instead of ^-1 if you know your matrix is positive
     * definite.
     */
    Matrix<T> &operator^= (const int &e)
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
				Matrix<T> L, U;
				Matrix<int> perm_vec;

				// step 1: compute the LU factorization
				if (A.rows_ == 1) {
					L = Matrix<T>(1);
					U = *(this);
					perm_vec = Matrix<int>(1,1);
				} else {
					int pivot;
					L = U = Matrix<T>(A.rows_, A.rows_);
					perm_vec = Matrix<int>(A.rows_ - 1, 1, false);

					for (int k = 0; k < A.rows_ - 1; ++k) {
						pivot = k;
						// find pivot
						for (int i = k; i < A.rows_; ++i) {
							if (::fabs(A(pivot,k)) < ::fabs(A(i, k)))
								pivot = i;
						}
						// check for singularity
						if (A(pivot, k) == (T) 0) {
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
							L(i,j) = (T) 0;
							L(i,i) = (T) 1;
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
					for (int ci = 0; ci < bb.rows_ - 1; ++ci) {
						int swap_row = static_cast<int>(perm_vec[ci]);
						for (int cj = 0; cj < bb.cols_; ++cj) {
							T temp = bb(ci,cj);
							bb(ci,cj) = bb(swap_row,cj);
							bb(swap_row,cj) = temp;
						}
					}
					/*
					Matrix<T> bb = b;
					for (int i = 0; i < bb.rows() - 1; ++i) {
						swap_ranges(bb.vec(i), bb.vec(i+1), bb.vec(perm_vec[i]));
					}*/

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
				throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
						__LINE__, "Invalid argument: -1");
			}
			return *this;
		}

    /**** Accessors ****/

    /* Return total length of the matrix */
    inline int size() const
    {
      return size_;
    }

    /* Return number of rows */
    inline int rows() const
    {
      return rows_;
    }

    /* Return number of columns */
    inline int cols() const
    {
      return cols_;
    }

    inline bool empty() const
    {
      // Equivalent to isNull but is the STL-compliant syntax
      return (rows_ == 0);
    }

    // Our max size is the biggest int or the limits of your mem
    inline int max_size () const {
      return INT_MAX;
    }

    /* Is this Matrix 0X0 (such a matrix is for array defs only) */
    inline bool isNull() const
    {
      // If rows_ == 0 so will cols_
      return (rows_ == 0);
    }

    /* Are all the elements in this Matrix == 0 */
    inline bool isZero() const
    {
      for (int i = 0; i < size_; ++i)
        if (data_[i] != 0)
          return false;
      return true;
    }
    
    /* 1X1 Matrix? */
    inline bool isScalar() const
    {
      return (rows_ == 1 && cols_ == 1);
    }

    /* 1Xm Matrix? */
    inline bool isRowVector() const
    {
      return (rows_ == 1);
    }

    /* nX1 Matrix? */
    inline bool isColVector() const
    {
      return (cols_ == 1);
    }
    
    /* nXn Matrix? Note that Null/Scalar Matrices are Square... */ 
    inline bool isSquare() const
    {
      return (rows_ == cols_);
    }
    
    /* M[i,j] == 0 when i != j */
    inline bool isDiagonal() const
		{
			for (int i = 0; i < rows_; ++i) {
				for (int j = 0; j < cols_; ++j) {
					if (i != j && data_[ijIndex(i,j)] != 0)
						return false;
				}
			}
			return true;
		}
    
    /* M[i,j] == 0 when i != j and 1 when i == j*/
    inline bool isIdentity() const
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

    /* M[i.j] == 0 when i < j */
    inline bool isLowerTriangular() const
		{
			for (int i = 0; i < rows_; ++i) {
				for (int j = i + 1; j < cols_; ++j) {
					if (data_[ijIndex(i,j)] != 0)
						return false;
				}
			}
			return true;
		}

    /*M[i,j] == 0 when i > j */
    inline bool isUpperTriangular() const
		{
			for (int j = 0; j < cols_; ++j) {
				for (int i = j + 1; i < rows_; ++i) {
					if (data_[ijIndex(i,j)] != 0)
						return false;
				}
			}
			return true;
		}

    /* This matrix has no inverse (iff its determinant == 0) */
    inline bool isSingular() const
		{
			if (! isSquare() || isNull())
				return false;
			if ((~(*(this))) == (T) 0)
				return true;
			return false;
		}

    /* This matrix is square and t(M) == M (inv(M) * t(M) == I) */
    inline bool isSymmetric() const
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

    /* This matrix is square and t(A) = -A */
    inline bool isSkewSymmetric() const
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

    /* Returns a pointer to the internal array.  User is responsible
     * for not messing things up.  This function should only be used
     * under special circumstances.
     */
    inline T *getArray() const
    {
      return data_;
    }

    /* Print matrix to a formatted string (very Java). Switched prec
     * and width from 0.1 cause prec is much more useful than width.
     * Also now can print out some internals.  Returns c++-style
     * string
     */
    std::string toString (const unsigned int &prec= 5,
                          const unsigned int &width = 0,
                          const bool &dim = false,
                          const bool &internal = false) const
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

    /* Save matrix to a file.  Flags can be "a":append,
     * "o":overwrite, "n":don't replace (all if file already
     * exists).  If header == true then the first elements written
     * are row and col size, and the matrix is saved as a flat
     * list.  Oherwise the matrix is written as a rectangular ascii
     * file. 
     * NOTE: Load is now a new type of constructor */
    void save (const std::string &path, const char &flag = 'n',
               const bool &header = 0, const int &prec = 5,
               const int &width = 0) const
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
				throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
						__LINE__,
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
    
    inline bool inRange (const int &i) const
    {
      return (i > -1 && i < size_ ? true : false);
    }
    
    inline bool inRange (const int &i, const int &j) const
    {
      int index = i * cols_ + j;
      return (index > -1 && index < size_ ? true : false);
    }

    /* Resizes the matrix by the ^2 1/2 alg.  */
    inline void resize (const int &rows, const int &cols,
                        const bool &fill=true)
    {
      if (rows < 0 || cols < 0)
        throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "Rows or cols < 0");

      resize(rows * cols, fill);
      rows_ = rows;
      cols_ = cols;
    }
    
    /**** Iterator Stuff ****/
    
    inline row_major_iterator begin()
    {
      return row_major_iterator(*this);
    }

    inline row_major_iterator end()
    {
      row_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline row_major_iterator vec(const int &n)
    {
      return row_major_iterator(*this).next_vec(n);
    }
    
    inline const_row_major_iterator begin() const
    {
      return const_row_major_iterator(*this);
    }

    inline const_row_major_iterator end() const
    {
      const_row_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline const_row_major_iterator vec(const int &n) const
    {
      return const_row_major_iterator(*this).next_vec(n);
    }

    inline col_major_iterator beginc()
    {
      return col_major_iterator(*this);
    }
    
    inline col_major_iterator endc()
    {
      col_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline col_major_iterator vecc(const int &n)
    {
      return col_major_iterator(*this).next_vec(n);
    }
    
    inline const_col_major_iterator beginc() const
    {
      return const_col_major_iterator(*this);
    }
    
    inline const_col_major_iterator endc() const
    {
      const_col_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline const_col_major_iterator vecc(const int &n) const
    {
      return const_col_major_iterator(*this).next_vec(n);
    }

    inline reverse_row_major_iterator rbegin ()
    {
      return reverse_row_major_iterator(*this);
    }
    
    inline reverse_row_major_iterator rend()
    {
      reverse_row_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline reverse_row_major_iterator rvec(const int &n)
    {
      return reverse_row_major_iterator(*this).next_vec(n);
    }
    
    inline const_reverse_row_major_iterator rbegin() const
    {
      return const_reverse_row_major_iterator(*this);
    }

    inline const_reverse_row_major_iterator rend() const
    {
      const_reverse_row_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline const_reverse_row_major_iterator rvec(const int &n) const
    {
      return const_reverse_row_major_iterator(*this).next_vec(n);
    }

    inline reverse_col_major_iterator rbeginc()
    {
      return reverse_col_major_iterator(*this);
    }
    
    inline reverse_col_major_iterator rendc()
    {
      reverse_col_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline reverse_col_major_iterator rvecc(const int &n)
    {
      return reverse_col_major_iterator(*this).next_vec(n);
    }
    
    inline const_reverse_col_major_iterator rbeginc() const
    {
      return const_reverse_col_major_iterator(*this);
    }
    
    inline const_reverse_col_major_iterator rendc() const
    {
      const_reverse_col_major_iterator temp(*this);
      return (temp + (size_));
    }
    
    inline const_reverse_col_major_iterator rvecc(const int &n) const
    {
      return const_reverse_col_major_iterator(*this).next_vec(n);
    }

  private:
    /**** Helper Functions ****/
    inline int ijIndex(const int &i, const int &j) const
    {
      return (i * cols_ + j);
    }
      
    inline void resize2Match (const Matrix<T> &m)
    {
      resize(m.size_, false);
      rows_ = m.rows_;
      cols_ = m.cols_;
    }

    inline void resize (const int &s, const bool &fill=true)
    {
      try {
        if (s > alloc_)
          grow(s, fill);
        else if (s < .25 * alloc_)
          shrink(fill);
      } catch (scythe_alloc_error &sae) {
        throw;
      }
      size_ = s;
    }

    inline void grow (const int &s, const bool &fill=true)
    {
      alloc_ = alloc_ ? alloc_ : 1;
      while (alloc_ < s)
        alloc_ <<= 1;

      T *temp = data_;
      data_ = new (std::nothrow) T[alloc_];
    
      if (data_ == 0) {
        throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "Failed to reallocate internal array");
      }
      
      if (fill) {
        for (int i =0; i < size_; ++i)
          data_[i] = temp[i];
      }
      
      delete[] temp;
    }
  
    inline void shrink (const bool &fill=true)
    {
      alloc_ >>= 1;
      //data_ = (T *) realloc(data_, sizeof(T) * alloc_);
      T *temp = data_;
      data_ = new (std::nothrow) T[alloc_];
      
      if (data_ == 0) {
        throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
           __LINE__, "Failed to reallocate internal array");
      }
      
      if (fill) {
        for (int i =0; i < alloc_; ++i)
          data_[i] = temp[i];
      }
      
      delete[] temp;
    }

    /**** Instance Variables ****/
    int rows_;        // # of rows
    int cols_;        // # of cols
    int size_;        // rows_ * cols_
    int alloc_;       // Total allocated size
    T *data_;         // The actual elements of the Matrix

  }; /* class Matrix */
      
  /***** (In)Equality operators ****/
  
  /* Matrix (in)equality (size and each element equal; <, > <=, >=
   * deal purely with size)
   */
  
  template <class T>
  bool operator== (const Matrix<T> &A, const Matrix<T> &B)
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

  template<class T>
  bool operator== (const Matrix<T> &A,
                   const typename Matrix<T>::ttype &b)
  {
    return (A == Matrix<T>(b));
  }

  template <class T>
  bool operator== (const typename Matrix<T>::ttype &a, 
                   const Matrix<T> &B)
  {
    return (Matrix<T>(a) == B);
  }

  template <class T>
  bool operator!= (const Matrix<T> &A, const Matrix<T> &B)
  {
    return !(A == B);
  }

  template<class T>
  bool operator!= (const Matrix<T> &A,
                   const typename Matrix<T>::ttype &b)
  {
    return !(A == Matrix<T>(b));
  }

  template <class T>
  bool operator!= (const typename Matrix<T>::ttype &a,
                   const Matrix<T> &B)
  {
    return !(Matrix<T>(a) == B);
  }
      
  template <class T>
  bool operator< (const Matrix<T> &A, const Matrix<T> &B)
  {
    if (A.size() < B.size())
      return true;
    return false;
  }
  
  template<class T>
  bool operator< (const Matrix<T> &A,
			            const typename Matrix<T>::ttype &b)
  {
    return (A < Matrix<T>(b));
  }

  template <class T>
  bool operator< (const typename Matrix<T>::ttype &a,
			            const Matrix<T> &B)
  {
    return (Matrix<T>(a) < B);
  }
      
  template <class T>
  bool operator> (const Matrix<T> &A, const Matrix<T> &B)
  {
    if (A.size() > B.size())
      return true;
    return false;
  }
  
  template<class T>
  bool operator> (const Matrix<T> &A,
			            const typename Matrix<T>::ttype &b)
  {
    return (A > Matrix<T>(b));
  }

  template <class T>
  bool operator> (const typename Matrix<T>::ttype &a,
			            const Matrix<T> &B)
  {
    return (Matrix<T>(a) > B);
  }
      
  template <class T>
  bool operator<= (const Matrix<T> &A, const Matrix<T> &B)
  {
    return (A.size() == B.size() || A < B );
  }
  
  template<class T>
  bool operator<= (const Matrix<T> &A,
			             const typename Matrix<T>::ttype &b)
  {
    return (A <= Matrix<T>(b));
  }

  template <class T>
  bool operator<= (const typename Matrix<T>::ttype &a,
			             const Matrix<T> &B)
  {
    return (Matrix<T>(a) <= B);
  }
      
  template <class T>
  bool operator>= (const Matrix<T> &A, const Matrix<T> &B)
  {
    return (A.size() == B.size() || A > B );
  }
  
  template<class T>
  bool operator>= (const Matrix<T> &A,
			             const typename Matrix<T>::ttype &b)
  {
    return (A >= Matrix<T>(b));
  }

  template <class T>
  bool operator>= (const typename Matrix<T>::ttype &a,
			             const Matrix<T> &B)
  {
    return (Matrix<T>(a) >= B);
  }
    
  /* Element-by-element (in)equality operators: return Matrices
   * filled with 0s and 1s.  Ex: If we have A |= B, the matrix
   * returned is always of the same dimensions as A and A conforms
   * with B in 4 cases: B.isScalar() A.rows() == B.rows() &&
   * A.cols() == B.cols() A.rows() == B.rows() && B.cols() == 1
   * A.cols() == B.cols() && B.rows() == 1
   *
	 * We define all three ((Matrix, Matrix), (Matrix, T), (T,
	 * Matrix)) although, in most cases, the (T, Matrix) option
	 * will throw an error.  It would be nice to rule this
	 * option out entirely and catch the problem at compile-time
	 * but such comparisons are valid when the Matrix
	 * passed in is 1x1
   */

  template <class T>
  Matrix<bool> operator|= (const Matrix<T> &A, const Matrix<T> &B)
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
  Matrix<bool> operator|= (const Matrix<T> &A,
			                     const typename Matrix<T>::ttype &b)
  {
    return (A |= Matrix<T>(b));
  }
  
  template <class T>
  Matrix<bool> operator|= (const typename Matrix<T>::ttype &a,
			                     const Matrix<T> &B)
  {
    return (Matrix<T>(a) |= B);
  }
  
  template <class T>
  Matrix<bool> operator| (const Matrix<T> &A, const Matrix<T> &B)
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
  Matrix<bool> operator| (const Matrix<T> &A,
			                    const typename Matrix<T>::ttype &b)
  {
    return (A | Matrix<T>(b));
  }
  
  template <class T>
  Matrix<bool> operator| (const typename Matrix<T>::ttype &a,
			                    const Matrix<T> &B)
  {
    return (Matrix<T>(a) | B);
  }
  
  template <class T>
  Matrix<bool> operator<< (const Matrix<T> &A, const Matrix<T> &B)
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
  Matrix<bool> operator<< (const Matrix<T> &A,
			                     const typename Matrix<T>::ttype &b)
  {
    return (A << Matrix<T>(b));
  }
  
  template <class T>
  Matrix<bool> operator<<(const typename Matrix<T>::ttype &a,
                          const Matrix<T> &B)
  {
    return (Matrix<T>(a) << B);
  }
  
  template <class T>
  Matrix<bool> operator>> (const Matrix<T> &A, const Matrix<T> &B)
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
  Matrix<bool> operator>>(const Matrix<T> &A,
                          const typename Matrix<T>::ttype &b)
  {
    return (A >> Matrix<T>(b));
  }
  
  template <class T>
  Matrix<bool> operator>>(const typename Matrix<T>::ttype &a,
                          const Matrix<T> &B)
  {
    return (Matrix<T>(a) >> B);
  }
  
  template <class T>
  Matrix<bool> operator<<= (const Matrix<T> &A, const Matrix<T> &B)
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
  Matrix<bool> operator<<= (const Matrix<T> &A,
                            const typename Matrix<T>::ttype &b)
  {
    return (A <<= Matrix<T>(b));
  }
  
  template <class T>
  Matrix<bool> operator<<= (const typename Matrix<T>::ttype &a,
                            const Matrix<T> &B)
  {
    return (Matrix<T>(a) <<= B);
  }
      
  template <class T>
  Matrix<bool> operator>>= (const Matrix<T> &A, const Matrix<T> &B)
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
  Matrix<bool> operator>>= (const Matrix<T> &A,
                            const typename Matrix<T>::ttype &b)
  {
    return (A >>= Matrix<T>(b));
  }
  
  template <class T>
  Matrix<bool> operator>>= (const typename Matrix<T>::ttype &a,
                            const Matrix<T> &B)
  {
    return (Matrix<T>(a) >>= B);
  }

  /**** Matrix arithmetic operators ****/
  /* When operations are commutative we pass both args in by reference
   * and check size before copying for (Matrix, Matrix) and always
   * copy the T where a T is involved.  For non-commutative, we pass
   * the first arg by value.  When possible, the user should try to
   * pass the smaller Matrix first although subsequent resizes may
   * make up for the time saved in copying in some cases.
   */

  /* Matrix addition */
  template <class T>
  Matrix<T> operator+ (const Matrix<T> &A, const Matrix<T> &B)
  {
    // If A or B is 1x1 this can save some time
    if (A.size() < B.size())
      return Matrix<T>(A) += B;
    return Matrix<T>(B) += A;
  }
  
  template <class T>
  Matrix<T> operator+(const Matrix<T> &A,
                      const typename Matrix<T>::ttype &b)
  {
    return Matrix<T>(b) += A;
  }
  
  template <class T>
  Matrix<T> operator+(const typename Matrix<T>::ttype &a,
                      const Matrix<T> &B)
  {
    return Matrix<T>(a) += B;
  }

  /* Matrix subtraction */
  template <class T>
  Matrix<T> operator- (Matrix<T> A, const Matrix<T> &B)
  {
    return A -= B;
	}
  
  template <class T>
  Matrix<T> operator- (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A -= Matrix<T>(b);
  }
  
  template <class T>
  Matrix<T> operator-(const typename Matrix<T>::ttype &a,
                      const Matrix<T> &B)
  {
    return Matrix<T>(a) -= B;
  }
    
  /* Negate all the elements in the matrix */
  template <class T>
  Matrix<T> operator- (Matrix<T> A)
  {
    for (int i = 0; i < A.size(); ++i)
      A[i] = -A[i];
    return A;
  }
    
  /* Matrix Multiplication */
  template <class T>
  Matrix<T> operator* (Matrix<T> A, const Matrix<T> &B)
  {
    return A *= B;
  }
  
  // 1X1 * nXm is commutative but we'd grow anyway so don't bother
  template <class T>
  Matrix<T> operator* (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A *= Matrix<T>(b);
  }
  
  template <class T>
  Matrix<T> operator*(const typename Matrix<T>::ttype &a,
                      const Matrix<T> &B)
  {
    return Matrix<T>(a) *= B;
  }
    
  /* Kronecker Multiplication */
  template <class T>
  Matrix<T> operator% (Matrix<T>A, const Matrix<T> &B)
  {
    return A %= B;
  }
  
  // commutative but we'd grow anyway
  template <class T>
  Matrix<T> operator% (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A %= Matrix<T>(b);
  }
  
  template <class T>
  Matrix<T> operator%(const typename Matrix<T>::ttype &a,
                      const Matrix<T> &B)
  {
    return Matrix<T>(a) %= B;
  }
    
  /* Element-by-element division */
  template <class T>
  Matrix<T> operator/ (Matrix<T> A, const Matrix<T> &B)
  {
    return A /= B;
  }
  
	template <class T>
  Matrix<T> operator/ (Matrix<T> A, const typename Matrix<T>::ttype &b)
  {
    return A /= Matrix<T>(b);
  }

  template <class T>
  Matrix<T> operator/ (const typename Matrix<T>::ttype &a,
										   const Matrix<T> &B)
  {
    return Matrix<T>(a) /= B;
  }
    
  /* Matrix power: ^0 returns identity matrix of this size, ^-1
   * returns inverse, otherwise must be a positive int
   */
  template <class T>
  Matrix<T> operator^ (Matrix<T> A, const int &i)
  {
    return A ^= i;
  }
    
  /* Return the transpose of this matrix */
  template <class T>
  Matrix<T> operator! (const Matrix<T> &M)
  {
    int rows = M.rows();
    int cols = M.cols();
    Matrix<T> ret(cols, rows, false);
    for (register int i = 0; i < rows; ++i) {
      for (register int j = 0; j < cols; ++j)
        //ret(j,i) = m(i,j);
        ret[j * rows + i] = M[i * cols + j];
    }
    return ret;
  }
    
  /* Return the determinant of a SQUARE matrix vi LU decomposition*/
  template <class T>
  T operator~ (Matrix <T> A) //no reference because LU kills the matrix
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



  template <class T>
    inline Matrix<T>
    r2scythe (const int &rows, const int &cols, const T *data)
    {
      Matrix<T> M(rows, cols, false);

      for (register int i = 0; i < cols; ++i) {
        for (register int j = 0; j < rows; ++j)
          M[j * cols + i] = data[i * rows + j];
      }

        return M;
    }

			
}  // end namespace SCYTHE

#endif /* SCYTHE_MATRIX_H */
