/* Scythe_Matrix.h
 *
 * This header provides the class definition of the Matrix class, part
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

#ifndef SCYTHE_MATRIX_H
#define SCYTHE_MATRIX_H

#include <string>
#include "Scythe_Error.h"
#include "Scythe_Util.h"
#include "Scythe_Matrix_Iterator.h"

namespace SCYTHE {
  
  inline  double sgn(const double& x){
    if (x>0)
      return 1.0;
    else if (x<0)
      return -1.0;
    else
      return 0.0;
  } 
  

  struct all_elements{
  }  const _ = {};
  
 
  
  enum IN_TYPE {	NORMAL, REPEAT, DIAG, UTRIANG, LTRIANG,
			BLOCK };

  template <class T>
    class Matrix
    {
    public:
      typedef T ttype;
      typedef matrix_iterator<ttype> iterator;
      typedef const_matrix_iterator<ttype> const_iterator;
      typedef row_major_iterator<ttype> row_major_iterator;
      typedef const_row_major_iterator<ttype> const_row_major_iterator;
      typedef col_major_iterator<ttype> col_major_iterator;
      typedef const_col_major_iterator<ttype> const_col_major_iterator;

      friend class iterator;
      friend class const_iterator;
      friend class row_major_iterator;
      friend class col_major_iterator;
      friend class const_row_major_iterator;
      friend class const_col_major_iterator;
			
      /**** Constructors ****/

      /* Default Constructor: Creates a Matrix of size 0.  This
       * Matrix cannot be used in operations but is useful when you
       * want to make arrays of matrices.
       */
      Matrix (); 

      /* Parameterized Type Constructor:  Create a 1x1 Matrix
       * containing one value.  While this is really just a scalar,
       * it has its uses:  Necessary for assignments such as Matrix A
       * = 3; or Matrix B = A[0];.  Also means we only have to define
       * operators for Matrix objects.
       */
      Matrix (const T&); 

      /* Standard Constructor:  Creates an n x m Matrix.  By default
       * it fills the Matrix with zeroes but you can turn this off or
       * change the fill value
       */
      explicit
	Matrix (const int &n, const int &m, const bool &fill = true,
		const T &fill_value = 0); 

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
      Matrix (const int &n, const int &m, const T*,
	      IN_TYPE type = NORMAL, const int &a = -1,
	      const int &b = -1, const int &c = -1,
	      const int &d = -1);
			
      /* Create a matrix from a file: first two elements should be
       * # row and # col.  Followed by a list of elements.  Limited
       * range-checking is done: if the constructor reaches eof before
       * extracting the expected number of elements, an error will be
       * thrown.  Row/col <= 0 will be caught.  If you
       * forget the row/col it will consider the first 2 numbers as
       * row and col.  If your file is longer than expected, no error
       * will be thrown.
       */
      Matrix (const std::string &path);

      /* Copy Constructor: Create a copy of an existing Matrix */
      Matrix (const Matrix<T> &);

      template <class S>
	Matrix (const Matrix<S> &);

      /**** Destructor ****/
      ~Matrix ();

      /**** Indexing Operators ****/

      /* Retrieve the ith element in row major order */
      inline T &operator[] (const int &i)
	{
	  if ( !inRange(i)) {
	    throw scythe_out_of_range_error (__FILE__,__PRETTY_FUNCTION__,
					     __LINE__, std::string("Index ") & i & 
					     " out of range");
	  }
	  return data_[i];
	}

      /* Retrieve the (i,j)th element */
      inline T &operator() (const int &i, const int &j)
	{
	  if (! inRange(i, j)) {
	    throw scythe_out_of_range_error(__FILE__,__PRETTY_FUNCTION__,
					    __LINE__, std::string("Index (") & i & "," &  j & 
					    ") out of range");
	  }
	  return data_[i * cols_ + j];
	}

      /* Versions of the above two for const Matrix objects */
      inline T &operator[] (const int &i) const
	{
	  if (! inRange(i)) {
	    throw scythe_out_of_range_error (__FILE__,__PRETTY_FUNCTION__,
					     __LINE__, std::string("Index ") & i &
					     " out of range");
	  }
	  return data_[i];
	}
			
      inline T &operator() (const int &i, const int &j) const
	{
	  if (! inRange(i, j)) {
	    throw scythe_out_of_range_error(__FILE__,__PRETTY_FUNCTION__,
					    __LINE__, std::string("Index (") & i & "," &  j & 
					    ") out of range");
	  }
	  return data_[i * cols_ + j];
	}

      /* SubMatrix operator: returns a new Matrix with a,b,c,d
       * defining the bounds of the block s.t 0 <= a <= c < rows_;
       * 0 <= b <= d < cols_.
       */
      Matrix<T> operator() (const int &a, const int &b, const int &c,
			    const int &d) const;

      
      /* functions for extracting all rows and all columns with the 
       * _ struct
       */
      Matrix<T> operator() (const all_elements& a, const int& j);
      Matrix<T> operator() (const int& i, const all_elements& a);

      /**** Self-modifying arithmetic operators ****/

      /* Assignment operator */
      Matrix<T> &operator= (const Matrix<T> &);

      template <class S>
	Matrix<T> &operator= (const Matrix<S> &);

      /* Matrix addition/assignment */
      Matrix<T> &operator+= (const Matrix<T> &);

      /* Matrix subtraction/assignment */
      Matrix<T> &operator-= (const Matrix<T> &);

      /* Matrix multiplication/assignment */
      Matrix<T> &operator*= (const Matrix<T> &);

      /* Kronecker Multiplication/assignment */
      Matrix<T> &operator%= (const Matrix<T> &);

      /* Element-by-element division/assignment */
      Matrix<T> &operator/= (const Matrix<T> &);

      /* Matrix power/assignment: ^ positive integer does matrix
       * power, ^0 returns an identity matrix of the input's size, ^-1
       * returns inverse via LU Decomposition.  Input must be
       * square for ^posint and ^-1.  You should use invpd in
       * Scythe_IDE instead of ^-1 if you know your matrix is positive
       * definite.
       */
      Matrix<T> &operator^= (const int &);
	
      /**** Accessors ****/

      /* Return total length of the matrix */
      inline int size() const
	{
	  return rows_ * cols_;
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

      /* Is this Matrxi 0X0 (such a matrix is for array defs only) */
      inline bool isNull() const
	{
	  // If rows_ == 0 so will cols_
	  return (rows_ == 0);
	}

      /* Are all the elements in this Matrix == 0 */
      inline bool isZero() const
	{
	  for (int i = 0; i < size(); ++i)
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
      bool isDiagonal() const;
			
      /* M[i,j] == 0 when i != j and 1 when i == j*/
      bool isIdentity() const;

      /* M[i.j] == 0 when i < j */
      bool isLowerTriangular() const;

      /*M[i,j] == 0 when i > j */
      bool isUpperTriangular() const;

      /* This matrix has no inverse (iff its determinant == 0) */
      bool isSingular() const;

      /* This matrix is square and t(M) == M (inv(M) * t(M) == I) */
      bool isSymmetric() const;

      /* This matrix is square and t(A) = -A */
      bool isSkewSymmetric() const;

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
			    const bool &internal = false) const;

      /* Save matrix to a file.  Flags can be "a":append,
       * "o":overwrite, "n":don't replace (all if file already
       * exists).  If header == true then the first elements written
       * are row and col size, and the matrix is saved as a flat
       * list.  Oherwise the matrix is written as a rectangular ascii
       * file. 
       * NOTE: Load is now a new type of constructor */
      void save(const std::string &path, const char &flag = 'n',
		const bool &header = 0, const int &prec = 5,
		const int &width = 0) const;
			
      inline bool inRange(const int &i) const
	{
	  return (i >= 0 && i < rows_ * cols_ ? true : false);
	}
			
      inline bool inRange (const int &i, const int &j) const
	{
	  return (i >= 0 && j >= 0 && i < rows_ && j < cols_
		  ? true : false);
	}

      /* Resize the matrix.  Note that end elements will be trimmed
       * off in row major order.  Use with caution.  Version leaves
       * junk at the end of a matrix that grows.  Version 2 fills it
       * with a default elements
       */
      inline void resize(const int &rows, const int &cols)
	{
	  if (rows < 0 || cols < 0)
	    throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
				     __LINE__, "Rows or cols < 0");

	  resize(rows * cols);
	  rows_ = rows;
	  cols_ = cols;
	}

      inline void resize(const int &rows, const int &cols, const T &def)
	{
	  int old_size = rows_ * cols_;
	  resize(rows, cols);
	  if (rows_ * cols_ > old_size) {
	    for (int i = old_size; i < rows_ * cols_; ++i)
	      data_[i] = def;
	  }
	}

      /**** Iterator Stuff ****/
			
      inline row_major_iterator begin()
	{
	  return row_major_iterator(*this);
	}

      inline row_major_iterator end()
	{
	  row_major_iterator temp(*this);
	  return (temp + (rows_ * cols_));
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
	  return (temp + (rows_ * cols_));
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
	  return (temp + (rows_ * cols_));
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
	  return (temp + (rows_ * cols_));
	}
			
      inline const_col_major_iterator vecc(const int &n) const
	{
	  return const_col_major_iterator(*this).next_vec(n);
	}

    private:
      /**** Helper Functions ****/
      inline int ijIndex(const int &i, const int &j) const
	{
	  return (i * cols_ + j);
	}
				
      inline int getAllocSize(const int &size) const
	{
	  if (size < 0) {
	    throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
				     __LINE__, "Can't allocate Matrix of size < 0");
	  } else if (size == 0) {
	    return 1;
	  } else if (size > alloc_) {
	    int x = 1;
	    while (size > x) {
	      x *= 2;
	    }
	    return x;
	  } else if (size < .25 * alloc_)
	    return (alloc_ / 2);
				
	  return alloc_;
	}
			
      inline void resize2Match (const Matrix<T> &m)
	{
	  resize(m.size());
	  rows_ = m.rows();
	  cols_ = m.cols();
	}
		
      inline void resize (const int &s)
	{
	  try {
	    if (s > size())
	      grow(s - size());
	    else if (s < size())
	      shrink(size() - s);
	  } catch (scythe_alloc_error &sae) {
	    throw;
	  }
	  // Ignore same size resizes
	}
		
      inline void grow (const int &extra)
	{
	  alloc_ = getAllocSize(size() + extra);
	  //data_ = (T *) realloc(data_, sizeof(T) * alloc_);
	  T *temp = data_;
	  data_ = new (std::nothrow) T[alloc_];
				
	  if (data_ == 0) {
	    throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
				     __LINE__, "Failed to reallocate internal array");
	  }
				
	  for (int i =0; i < rows_ * cols_; ++i)
	    data_[i] = temp[i];
				
	  delete[] temp;
	}
		
      inline void shrink (const int &difference)
	{
	  alloc_ = getAllocSize(size() - difference);
	  //data_ = (T *) realloc(data_, sizeof(T) * alloc_);
	  T *temp = data_;
	  data_ = new (std::nothrow) T[alloc_];
				
	  if (data_ == 0) {
	    throw scythe_alloc_error(__FILE__, __PRETTY_FUNCTION__,
				     __LINE__, "Failed to reallocate internal array");
	  }
				
	  for (int i =0; i < alloc_; ++i)
	    data_[i] = temp[i];
				
	  delete[] temp;
	}

      /**** Instance Variables ****/
      int rows_;				// # of rows
      int cols_;				// # of cols
      int alloc_;				// Total allocated size
      T *data_;						// The actual elements of the Matrix

    }; /* class Matrix */
			
  /***** (In)Equality operators ****/
	
  /* Matrix (in)equality (size and each element equal; <, > <=, >=
   * deal purely with size)
   */
	
  template <class T>
    bool operator== (const Matrix<T> &, const Matrix<T> &);

  template<class T>
    bool operator== (const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
    bool operator== (const typename Matrix<T>::ttype &, const Matrix<T> &);

  template <class T>
    bool operator!= (const Matrix<T> &, const Matrix<T> &);

  template<class T>
    bool operator!= (const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
    bool operator!= (const typename Matrix<T>::ttype &, const Matrix<T> &);
			
  template <class T>
    bool operator< (const Matrix<T> &, const Matrix<T> &);
	
  template<class T>
    bool operator< (const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
    bool operator< (const typename Matrix<T>::ttype &, const Matrix<T> &);
			
  template <class T>
    bool operator> (const Matrix<T> &, const Matrix<T> &);
	
  template<class T>
    bool operator> (const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
    bool operator> (const typename Matrix<T>::ttype &, const Matrix<T> &);
			
  template <class T>
    bool operator<= (const Matrix<T> &, const Matrix<T> &);
	
  template<class T>
    bool operator<= (const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
    bool operator<= (const typename Matrix<T>::ttype &, const Matrix<T> &);
			
  template <class T>
    bool operator>= (const Matrix<T> &, const Matrix<T> &);
	
  template<class T>
    bool operator>= (const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
    bool operator>= (const typename Matrix<T>::ttype &, const Matrix<T> &);
		
  /* Element-by-element (in)equality operators: return Matrices
   * filled with 0s and 1s.  Ex: If we have A |= B, the matrix
   * returned is always of the same dimensions as A and A conforms
   * with B in 4 cases: B.isScalar() A.rows() == B.rows() &&
   * A.cols() == B.cols() A.rows() == B.rows() && B.cols() == 1
   * A.cols() == B.cols() && B.rows() == 1
   */
  template <class T>
    Matrix<bool> operator|= (const Matrix<T> &, const Matrix<T> &);

  template <class T>
    Matrix<bool> operator|= (const Matrix<T> &, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<bool> operator|= (const typename Matrix<T>::ttype &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator| (const Matrix<T> &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator| (const Matrix<T> &, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<bool> operator| (const typename Matrix<T>::ttype &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator<< (const Matrix<T> &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator<< (const Matrix<T> &, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<bool> operator<< (const typename Matrix<T>::ttype &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator>> (const Matrix<T> &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator>> (const Matrix<T> &, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<bool> operator>> (const typename Matrix<T>::ttype &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator<<= (const Matrix<T> &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator<<=(const Matrix<T> &, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<bool> operator<<=(const typename Matrix<T>::ttype &, const Matrix<T> &);
			
  template <class T>
    Matrix<bool> operator>>= (const Matrix<T> &, const Matrix<T> &);
	
  template <class T>
    Matrix<bool> operator>>=(const Matrix<T> &, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<bool> operator>>=(const typename Matrix<T>::ttype &, const Matrix<T> &);

  /**** Matrix arithmetic operators ****/

  /* Matrix addition */
  template <class T>
    Matrix<T> operator+ (const Matrix<T> &, const Matrix<T> &);
	
  template <class T>
    Matrix<T> operator+ (const Matrix<T> &, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<T> operator+ (const typename Matrix<T>::ttype &, const Matrix<T> &);

  /* Matrix subtraction */
  template <class T>
    Matrix<T> operator- (Matrix<T>, const Matrix<T> &);
	
  template <class T>
    Matrix<T> operator- (Matrix<T>, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<T> operator- (const typename Matrix<T>::ttype &, const Matrix<T> &);
		
  /* Negate all the elements in the matrix */
  template <class T>
    Matrix<T> operator- (Matrix<T>);
		
  /* Matrix Multiplication */
  template <class T>
    Matrix<T> operator* (Matrix<T>, const Matrix<T> &);
	
  template <class T>
    Matrix<T> operator* (Matrix<T>, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<T> operator* (const typename Matrix<T>::ttype &, const Matrix<T> &);
		
  /* Kronecker Multiplication */
  template <class T>
    Matrix<T> operator% (Matrix<T>, const Matrix<T> &);
	
  template <class T>
    Matrix<T> operator% (Matrix<T>, const typename Matrix<T>::ttype &);
	
  template <class T>
    Matrix<T> operator% (const typename Matrix<T>::ttype &, const Matrix<T> &);
		
  /* Element-by-element division */
  template <class T>
    Matrix<T> operator/ (Matrix<T>, const Matrix<T> &);
		
  /* Matrix power: ^0 returns identity matrix of this size, ^-1
   * returns inverse, otherwise must be a positive int
   */
  template <class T>
    Matrix<T> operator^ (Matrix<T>, const int &);
		
  /* Return the transpose of this matrix */
  template <class T>
    Matrix<T> operator! (const Matrix<T> &);
		
  /* Return the determinant of a SQUARE matrix vi LU decomposition*/
  template <class T>
    T operator~ (Matrix<T>);
		

}	// end namespace SCYTHE
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || \
    defined (EXPLICIT_TEMPLATE_INSTANTIATION)
// Necessary for template instantiation with some compilers.
# include "Scythe_Matrix.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_MATRIX_H */
