/* Scythe_Matrix.h
 *
 * This header provides the class definition of the Matrix class, part
 * of the SCYTHE project.
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

#ifndef SCYTHE_MATRIX_H
#define SCYTHE_MATRIX_H

#include <string>
#include <climits>
#include <cmath>
#include "Scythe_Error.h"
#include "Scythe_Util.h"
#include "Scythe_Matrix_Iterator.h"

#define SCYTHE_10log2 .30102999

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
    Matrix (const Matrix<T> &, const bool &fill=true);

    template <class S>
    Matrix (const Matrix<S> &);

    /**** Destructor ****/
    ~Matrix ();

    /**** STL container modifiers ****/
    /* Swap operator (sort of a dual copy constructor) */
    void swap (Matrix<T> &);

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
                          const int &d) const;

    
    /* functions for extracting all rows and all columns with the 
     * _ struct
     */
    Matrix<T> operator() (const all_elements& a, const int& j) const;
    Matrix<T> operator() (const int& i, const all_elements& a) const;

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

    /* Is this Matrxi 0X0 (such a matrix is for array defs only) */
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
    void save (const std::string &path, const char &flag = 'n',
               const bool &header = 0, const int &prec = 5,
               const int &width = 0) const;
    
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
  bool operator== (const Matrix<T> &, const Matrix<T> &);

  template<class T>
  bool operator==(const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
  bool operator==(const typename Matrix<T>::ttype &, const Matrix<T> &);

  template <class T>
  bool operator!= (const Matrix<T> &, const Matrix<T> &);

  template<class T>
  bool operator!=(const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
  bool operator!=(const typename Matrix<T>::ttype &, const Matrix<T> &);
      
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
  bool operator<=(const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
  bool operator<=(const typename Matrix<T>::ttype &, const Matrix<T> &);
      
  template <class T>
  bool operator>= (const Matrix<T> &, const Matrix<T> &);
  
  template<class T>
  bool operator>=(const Matrix<T> &, const typename Matrix<T>::ttype &);

  template <class T>
  bool operator>=(const typename Matrix<T>::ttype &, const Matrix<T> &);
    
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
  Matrix<bool> operator|=(const Matrix<T> &,
                          const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<bool> operator|=(const typename Matrix<T>::ttype &,
                          const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator| (const Matrix<T> &, const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator| (const Matrix<T> &,
                          const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<bool> operator| (const typename Matrix<T>::ttype &,
                          const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator<< (const Matrix<T> &, const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator<<(const Matrix<T> &,
                          const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<bool> operator<<(const typename Matrix<T>::ttype &,
                          const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator>> (const Matrix<T> &, const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator>>(const Matrix<T> &,
                          const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<bool> operator>>(const typename Matrix<T>::ttype &,
                          const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator<<= (const Matrix<T> &, const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator<<= (const Matrix<T> &,
                            const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<bool> operator<<= (const typename Matrix<T>::ttype &,
                            const Matrix<T> &);
      
  template <class T>
  Matrix<bool> operator>>= (const Matrix<T> &, const Matrix<T> &);
  
  template <class T>
  Matrix<bool> operator>>= (const Matrix<T> &,
                            const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<bool> operator>>= (const typename Matrix<T>::ttype &,
                            const Matrix<T> &);

  /**** Matrix arithmetic operators ****/

  /* Matrix addition */
  template <class T>
  Matrix<T> operator+ (const Matrix<T> &, const Matrix<T> &);
  
  template <class T>
  Matrix<T> operator+(const Matrix<T> &,
                      const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<T> operator+(const typename Matrix<T>::ttype &,
                      const Matrix<T> &);

  /* Matrix subtraction */
  template <class T>
  Matrix<T> operator- (Matrix<T>, const Matrix<T> &);
  
  template <class T>
  Matrix<T> operator- (Matrix<T>, const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<T> operator-(const typename Matrix<T>::ttype &,
                      const Matrix<T> &);
    
  /* Negate all the elements in the matrix */
  template <class T>
  Matrix<T> operator- (Matrix<T>);
    
  /* Matrix Multiplication */
  template <class T>
  Matrix<T> operator* (Matrix<T>, const Matrix<T> &);
  
  template <class T>
  Matrix<T> operator* (Matrix<T>, const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<T> operator*(const typename Matrix<T>::ttype &,
                      const Matrix<T> &);
    
  /* Kronecker Multiplication */
  template <class T>
  Matrix<T> operator% (Matrix<T>, const Matrix<T> &);
  
  template <class T>
  Matrix<T> operator% (Matrix<T>, const typename Matrix<T>::ttype &);
  
  template <class T>
  Matrix<T> operator%(const typename Matrix<T>::ttype &,
                      const Matrix<T> &);
    
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

    
}  // end namespace SCYTHE
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || \
    defined (EXPLICIT_TEMPLATE_INSTANTIATION)
// Necessary for template instantiation with some compilers.
# include "Scythe_Matrix.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_MATRIX_H */
