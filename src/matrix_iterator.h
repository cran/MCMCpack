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
 * scythestat/matrix_iterator.h
 *
 * Provides definitions and implementations of iterators for
 * the Matrix class.  These iterators conform to the
 * expectations of functions from the C++ Standard Template
 * Library (STL) which operate over ranges within
 * container-type objects.
 *
 */

#ifndef SCYTHE_MATRIX_ITERATOR_H
#define SCYTHE_MATRIX_ITERATOR_H

#include <iterator>
#include <cmath>

#ifdef SCYTHE_COMPILE_DIRECT
#include "error.h"
#include "util.h"
#include "matrix.h"
#else
#include "scythestat/error.h"
#include "scythestat/util.h"
#include "scythestat/matrix.h"
#endif

namespace SCYTHE {

  template <class T>
  class Matrix;

  template <class T>
  class const_matrix_iterator;

  template <class T>
  class matrix_iterator
    : public std::iterator<std::random_access_iterator_tag, T>
  {
  public:
    friend class const_matrix_iterator<T>;
  
    virtual ~matrix_iterator ()
    {
    }

    /**** Forward iterator Facilities ****/

    // You should define an equals operator of the form:
    // matrix_iterator<T> &operator=(const matrix_iterator &);
    // in extending classes
      
    // Provide read/write access to referent
    inline T &operator* () const
    {
      return matrix_->data_[current_];
    }

    // Provide read access to a member (if any) of referent
    inline T *operator-> () const
    {
      return &(matrix_->data_[current_]);
    }

    /* Pre and postfix operators: note that the postfix operators
     * return non-reference values and therefore aren't technically
     * part of the interface because only pointers or references may
     * be covariant return values.
     */

    virtual matrix_iterator<T> &operator++ () = 0;

    //virtual matrix_iterator<T> operator++ (int) = 0;
      
    /**** Bidirectional Iterator Facilities ****/

    virtual matrix_iterator<T> &operator-- () = 0;

    //virtual matrix_iterator<T> operator-- (int) = 0;
      
    /**** Random Access Iterator Facilities ****/

    virtual T &operator[] (const int &) const = 0;

    virtual matrix_iterator<T> &operator+= (const int &) = 0;

    virtual matrix_iterator<T> &operator-= (const int &) = 0;

    // Extending classes should provide a distance operator of the
    // form:
    // std::ptrdiff_t operator- (const matrix_iterator<T> &rmi)


    /**** Matrix Iterator Facilities ****/

    virtual matrix_iterator<T> &plus_vec () = 0;

    virtual matrix_iterator<T> &plus_vec (const int &) = 0;

    virtual matrix_iterator<T> &minus_vec () = 0;

    virtual matrix_iterator<T> &minus_vec (const int &) = 0;

    virtual matrix_iterator<T> &next_vec () = 0;

    virtual matrix_iterator<T> &next_vec (const int &) = 0;

    virtual matrix_iterator<T> &prev_vec () = 0;

    virtual matrix_iterator<T> &prev_vec (const int &) = 0;

    /**** Figure out our current index ****/
    int get_row() const
    {
      return (int) (current_ / matrix_->cols());
    }

    int get_col() const
    {
      int row = (int) (current_ / matrix_->cols());
      return (current_ - (row * matrix_->cols()));
    }

    int get_index () const
    {
      return current_;
    }

  protected:
    matrix_iterator()
      :  matrix_ (0),
        current_ (0)
    {
    }
      
    explicit matrix_iterator (Matrix<T> &m)
      :  matrix_ (&m),
        current_ (0)
    {
    }

    matrix_iterator (const matrix_iterator<T> &mi)
      :  matrix_ (mi.matrix_),
        current_ (mi.current_)
    {
    }

    Matrix<T> *matrix_;
    int current_;
  };

  template <class T>
  class const_matrix_iterator
    : public std::iterator<std::random_access_iterator_tag, T>
  {
  public:
    virtual ~const_matrix_iterator ()
    {
    }

    /**** Forward iterator Facilities ****/

    // You should define an equals operator of the form:
    // matrix_iterator<T> &operator=(const matrix_iterator &);
    // in extending classes
    
    // Provide read/write access to referent
    inline const T &operator* () const
    {
      return matrix_->data_[current_];
    }

    // Provide read access to a member (if any) of referent
    inline const T *operator-> () const
    {
      return &(matrix_->data_[current_]);
    }

    /* Pre and postfix operators: note that the postfix operators
     * return non-reference values and therefore aren't technically
     * part of the interface because only pointers or references may
     * be covariant return values.
     */

    virtual const_matrix_iterator<T> &operator++ () = 0;

    //virtual const_matrix_iterator<T> operator++ (int) = 0;
    
    /**** Bidirectional Iterator Facilities ****/

    virtual const_matrix_iterator<T> &operator-- () = 0;

    //virtual const_matrix_iterator<T> operator-- (int) = 0;
    
    /**** Random Access Iterator Facilities ****/

    virtual const T &operator[] (const int &) const = 0;

    virtual const_matrix_iterator<T> &operator+= (const int &) = 0;

    virtual const_matrix_iterator<T> &operator-= (const int &) = 0;

    // Extending classes should provide a distance operator of the
    // form:
    // std::ptrdiff_t operator- (const const_matrix_iterator<T> &rmi)


    /**** Matrix Iterator Facilities ****/

    virtual const_matrix_iterator<T> &plus_vec () = 0;

    virtual const_matrix_iterator<T> &plus_vec (const int &) = 0;

    virtual const_matrix_iterator<T> &minus_vec () = 0;

    virtual const_matrix_iterator<T> &minus_vec (const int &) = 0;

    virtual const_matrix_iterator<T> &next_vec () = 0;

    virtual const_matrix_iterator<T> &next_vec (const int &) = 0;

    virtual const_matrix_iterator<T> &prev_vec () = 0;

    virtual const_matrix_iterator<T> &prev_vec (const int &) = 0;

    /**** Figure out our current index ****/
    int get_row() const
    {
      return (int) (current_ / matrix_->cols());
    }

    int get_col() const
    {
      int row = (int) (current_ / matrix_->cols());
      return (current_ - (row * matrix_->cols()));
    }

    int get_index () const
    {
      return current_;
    }


  protected:
    const_matrix_iterator()
      :  matrix_ (0),
        current_ (0)
    {
    }
    
    explicit const_matrix_iterator (const Matrix<T> &m)
      :  matrix_ (&m),
        current_ (0)
    {
    }

    const_matrix_iterator (const const_matrix_iterator<T> &mi)
      :  matrix_ (mi.matrix_),
        current_ (mi.current_)
    {
    }

    const_matrix_iterator (const matrix_iterator<T> &mi)
      : matrix_ (mi.matrix_),
        current_ (mi.current_)
    {
    }

    const Matrix<T> *matrix_;
    int current_;
  };

  /**** An iterator that does ops in row-major order ****/
  template <class T>
  class row_major_iterator : public matrix_iterator<T>
  {  
  public:

    /**** Constructors ****/

    row_major_iterator ()
      : matrix_iterator<T> ()
    {
    }

    explicit row_major_iterator (Matrix<T> &m)
      : matrix_iterator<T> (m)
    {
    }

    row_major_iterator (const row_major_iterator &rmi)
      : matrix_iterator<T> (rmi)
    {
    }

    virtual ~row_major_iterator ()
    {
    }

    /**** Forward Iterator Facilities ****/

    // Assignment operator
    inline row_major_iterator<T> &operator= (const
        row_major_iterator &rmi)
    {
      this->matrix_ = rmi.matrix_;
      this->current_ = rmi.current_;

      return *this;
    }
    
    // Step forward, return new position
    inline row_major_iterator<T> &operator++ ()
    {
      if (this->current_ < this->matrix_->size())
        ++this->current_;
      
      return *this;
    }

    // Step forward, return old position
    inline row_major_iterator<T> operator++ (int)
    {
      row_major_iterator<T> temp = *this;
      ++(*this);
      
      return temp;
    }

    /**** BiDirectional Iterator Facilities ****/

    // Step back, return new position
    inline row_major_iterator<T> &operator-- ()
    {
      if (this->current_ > 0)
        --this->current_;
      
      return *this;
    }

    // Step back, return old position
    inline row_major_iterator<T> operator-- (int)
    {
      row_major_iterator temp = *this;
      --(*this);
      
      return temp;
    }

    /**** Random Access Iterator Facilities ****/
    
    // Provide access to the [nth] element XXX int?
    inline T &operator[] (const int &n) const
    {
      return this->matrix_->data_[n];
    }

    // Step n elements
    inline row_major_iterator<T> &operator+= (const int &n)
    {
      if (this->current_ + n > this->matrix_->size())
        this->current_ = this->matrix_->size();
      else if (this->current_ + n < 0)
        this->current_ = 0;
      else
        this->current_ += n;
      
      return *this;
    }

    inline row_major_iterator<T> &operator-= (const int &n)
    {
      return (*this += -n);
    }

    /* Difference operators (for distance) */

    inline std::ptrdiff_t operator-
      (const row_major_iterator<T> &rmi) const
    {
      return this->current_ - rmi.current_;
    }

    /**** Matrix Iterator Facilities ****/

    // Jump forward the length of a row
    inline row_major_iterator<T> &plus_vec ()
    {
      return (*this += this->matrix_->cols());
    }

    // Jump forward the length of a row n times
    inline row_major_iterator<T> &plus_vec (const int &n)
    {
      return (*this += (n * this->matrix_->cols()));
    }

    // Jump backward the length of a row
    inline row_major_iterator<T> &minus_vec ()
    {
      return (*this -= this->matrix_->cols());
    }

    // Jump backward the length of a row n times
    inline row_major_iterator<T> &minus_vec (const int &n)
    {
      return (*this -= (n * this->matrix_->cols()));
    }

    // Jump to the beginnin of the next vector
    inline row_major_iterator<T> &next_vec ()
    {
      int cur_vec = (int) (this->current_ / this->matrix_->cols());
      if (cur_vec + 1 < this->matrix_->rows())
        this->current_ = (cur_vec + 1) * this->matrix_->cols();
      else
        this->current_ = this->matrix_->size();

      return *this;
    }

    // Jump to the beginning of the nth next vector
    inline row_major_iterator<T> &next_vec (const int &n)
    {
      int cur_vec = (int) (this->current_ / this->matrix_->cols());
      if (cur_vec + n >= this->matrix_->rows())
        this->current_ = this->matrix_->size();
      else if (cur_vec + n <= 0)
        this->current_ = 0;
      else
        this->current_ = (cur_vec + n) * this->matrix_->cols();

      return *this;
    }
      
    // Jump to the beginnin of the previous vector
    inline row_major_iterator<T> &prev_vec ()
    {
      int cur_vec = (int) (this->current_ / this->matrix_->cols());
      if (cur_vec > 0)
        this->current_ = (cur_vec - 1) * this->matrix_->cols();
      else
        this->current_ = 0;

      return *this;
    }

    // Jump to the beginning of the nth previous vector
    inline row_major_iterator<T> &prev_vec (const int &n)
    {
      return (this->next_vec(-n));
    }
    
    friend bool operator== (const row_major_iterator<T> &a,
                            const row_major_iterator<T> &b)
    {
      if (a.current_ == b.current_ && a.matrix_ == b.matrix_)
        return true;

      return false;
    }

    friend bool operator<(const row_major_iterator &a,
                          const row_major_iterator &b)
    {
      if (a.matrix_ != b.matrix_)
        throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, 
            "< Comparison on iterators to different matrices");
      
      if (a.current_ < b.current_)
        return true;

      return false;
    }

  };
  
  template <class T>
  class const_row_major_iterator : public const_matrix_iterator<T>
  {  
  public:

    /**** Constructors ****/

    const_row_major_iterator ()
      : const_matrix_iterator<T> ()
    {
    }

    explicit const_row_major_iterator (const Matrix<T> &m)
      : const_matrix_iterator<T> (m)
    {
    }

    const_row_major_iterator (const const_row_major_iterator<T> &rmi)
      : const_matrix_iterator<T> (rmi)
    {
    }
    
    const_row_major_iterator (const row_major_iterator<T> &rmi)
      : const_matrix_iterator<T> (rmi)
    {
    }

    virtual ~const_row_major_iterator ()
    {
    }

    /**** Forward Iterator Facilities ****/

    // Assignment operator
    inline const_row_major_iterator<T> &operator= (const
        const_row_major_iterator &rmi)
    {
      this->matrix_ = rmi.matrix_;
      this->current_ = rmi.current_;

      return *this;
    }
    
    // Step forward, return new position
    inline const_row_major_iterator<T> &operator++ ()
    {
      if (this->current_ < this->matrix_->size())
        ++this->current_;
      
      return *this;
    }

    // Step forward, return old position
    inline const_row_major_iterator<T> operator++ (int)
    {
      row_major_iterator<T> temp = *this;
      ++(*this);
      
      return temp;
    }

    /**** BiDirectional Iterator Facilities ****/

    // Step back, return new position
    inline const_row_major_iterator<T> &operator-- ()
    {
      if (this->current_ > 0)
        --this->current_;
      
      return *this;
    }

    // Step back, return old position
    inline const_row_major_iterator<T> operator-- (int)
    {
      const_row_major_iterator temp = *this;
      --(*this);
      
      return temp;
    }

    /**** Random Access Iterator Facilities ****/
    
    // Provide access to the [nth] element XXX int?
    inline const T &operator[] (const int &n) const
    {
      return this->matrix_->data_[n];
    }

    // Step n elements
    inline const_row_major_iterator<T> &operator+= (const int &n)
    {
      if (this->current_ + n > this->matrix_->size())
        this->current_ = this->matrix_->size();
      else if (this->current_ + n < 0)
        this->current_ = 0;
      else
        this->current_ += n;
      
      return *this;
    }

    inline const_row_major_iterator<T> &operator-= (const int &n)
    {
      return (*this += -n);
    }

    /* Difference operators (for distance) */

    inline std::ptrdiff_t operator-
      (const const_row_major_iterator<T> &rmi) const
    {
      return this->current_ - rmi.current_;
    }

    /**** Matrix Iterator Facilities ****/

    // Jump forward the length of a row
    inline const_row_major_iterator<T> &plus_vec ()
    {
      return (*this += this->matrix_->cols());
    }

    // Jump forward the length of a row n times
    inline const_row_major_iterator<T> &plus_vec (const int &n)
    {
      return (*this += (n * this->matrix_->cols()));
    }

    // Jump backward the length of a row
    inline const_row_major_iterator<T> &minus_vec ()
    {
      return (*this -= this->matrix_->cols());
    }

    // Jump backward the length of a row n times
    inline const_row_major_iterator<T> &minus_vec (const int &n)
    {
      return (*this -= (n * this->matrix_->cols()));
    }

    // Jump to the beginnin of the next vector
    inline const_row_major_iterator<T> &next_vec ()
    {
      int cur_vec = (int) (this->current_ / this->matrix_->cols());
      if (cur_vec + 1 < this->matrix_->rows())
        this->current_ = (cur_vec + 1) * this->matrix_->cols();
      else
        this->current_ = this->matrix_->size();

      return *this;
    }

    // Jump to the beginning of the nth next vector
    inline const_row_major_iterator<T> &next_vec (const int &n)
    {
      int cur_vec = (int) (this->current_ / this->matrix_->cols());
      if (cur_vec + n >= this->matrix_->rows())
        this->current_ = this->matrix_->size();
      else if (cur_vec + n <= 0)
        this->current_ = 0;
      else
        this->current_ = (cur_vec + n) * this->matrix_->cols();

      return *this;
    }
      
    // Jump to the beginnin of the previous vector
    inline const_row_major_iterator<T> &prev_vec ()
    {
      int cur_vec = (int) (this->current_ / this->matrix_->cols());
      if (cur_vec > 0)
        this->current_ = (cur_vec - 1) * this->matrix_->cols();
      else
        this->current_ = 0;

      return *this;
    }

    // Jump to the beginning of the nth previous vector
    inline const_row_major_iterator<T> &prev_vec (const int &n)
    {
      return (this->next_vec(-n));
    }
    
    friend bool operator== (const const_row_major_iterator<T> &a,
                            const const_row_major_iterator<T> &b)
    {
      if (a.current_ == b.current_ && a.matrix_ == b.matrix_)
        return true;

      return false;
    }

    friend bool operator<(const const_row_major_iterator &a,
                          const const_row_major_iterator &b)
    {
      if (a.matrix_ != b.matrix_)
        throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, 
            "< Comparison on iterators to different matrices");
      
      if (a.current_ < b.current_)
        return true;

      return false;
    }

  };
  
  /**** An iterator that does ops in col-major order ****/
  template <class T>
  class col_major_iterator : public matrix_iterator<T>
  {  
  public:

    /**** Constructors ****/

    col_major_iterator ()
      : matrix_iterator<T> ()
    {
    }

    explicit col_major_iterator (Matrix<T> &m)
      : matrix_iterator<T> (m)
    {
    }

    col_major_iterator (const col_major_iterator &cmi)
      : matrix_iterator<T> (cmi)
    {
    }

    virtual ~col_major_iterator ()
    {
    }

    /**** Forward Iterator Facilities ****/

    // Assignment operator
    inline col_major_iterator<T> &operator= (const
        col_major_iterator &cmi)
    {
      this->matrix_ = cmi.matrix_;
      this->current_ = cmi.current_;

      return *this;
    }
    
    // Step forward, return new position
    inline col_major_iterator<T> &operator++ ()
    {
      if (this->current_ >= this->matrix_->cols() * (this->matrix_->rows() - 1)) {
        if (this->current_ >= this->matrix_->size() - 1)
          this->current_ = this->matrix_->size();
        else
          this->current_ = (this->current_ + 1) -
            (this->matrix_->rows() - 1) * this->matrix_->cols();
      } else
        this->current_ += this->matrix_->cols();

      return *this;
    }

    // Step forward, return old position
    inline col_major_iterator<T> operator++ (int)
    {
      col_major_iterator<T> temp = *this;
      ++(*this);
      return temp;
    }

    /**** BiDirectional Iterator Facilities ****/

    // Step back, return new position
    inline col_major_iterator<T> &operator-- ()
    {
      if (this->current_ > 0) {
        if (this->current_ == this->matrix_->size())
          --this->current_;
        else if (this->current_ < this->matrix_->cols()) {
          this->current_ = (this->current_ - 1) + 
            (this->matrix_->rows() - 1) * this->matrix_->cols();
        } else
          this->current_ -= this->matrix_->cols();
      }
      
      return *this;
    }

    // Step back, return old position
    inline col_major_iterator<T> operator-- (int)
    {
      col_major_iterator temp = *this;
      --(*this);
      return temp;
    }

    /**** Random Access Iterator Facilities ****/
    
    // Provide access to the [nth] element XXX int?
    inline T &operator[] (const int &n) const
    {
      int col = (int) (n / this->matrix_->rows());
      int row = n - (col * this->matrix_->rows());

      return this->matrix_->data_[row * this->matrix_->cols_ + col];
    }

    // Step n elements
    inline col_major_iterator<T> &operator+= (const int &n)
    {
      int cm;

      if (this->current_ == this->matrix_->size())
        cm = this->current_;
      else {
        int row = (int) (this->current_ / this->matrix_->cols());
        int col = this->current_ - (row * this->matrix_->cols());
        cm = col * this->matrix_->rows() + row;
      }

      cm += n;

      if (cm >= this->matrix_->size())
        this->current_ = this->matrix_->size();
      else if (cm <= 0)
        this->current_ = 0;
      else {
        int col = (int) (cm / this->matrix_->rows());
        int row = cm - (col * this->matrix_->rows());
        this->current_ = row * this->matrix_->cols() + col;
      }

      return *this;
    }

    inline col_major_iterator<T> &operator-= (const int &n)
    {
      return (*this += -n);
    }

    /* Difference operators (for distance) */

    inline std::ptrdiff_t operator-
      (const col_major_iterator<T> &cmi) const
    {
      int cm, bcm;
      if (this->current_ == this->matrix_->size())
        cm = this->current_;
      else {
        int row = (int) (this->current_ / this->matrix_->cols());
        int col = this->current_ - (row * this->matrix_->cols());
        cm = col * this->matrix_->rows() + row;
      }

      if (cmi.current_ == this->matrix_->size())
        bcm = cmi.current_;
      else {
        int brow = (int) (cmi.current_ / this->matrix_->cols());
        int bcol = cmi.current_ - (brow * this->matrix_->cols());
        bcm = bcol * this->matrix_->rows() + brow;
      }
      
      return cm - bcm;
    }

    /**** Matrix Iterator Facilities ****/

    // Jump forward the length of a row
    inline col_major_iterator<T> &plus_vec ()
    {
      return (*this += this->matrix_->rows());
    }

    // Jump forward the length of a row n times
    inline col_major_iterator<T> &plus_vec (const int &n)
    {
      return (*this += (n * this->matrix_->rows()));
    }

    // Jump backward the length of a row
    inline col_major_iterator<T> &minus_vec ()
    {
      return (*this -= this->matrix_->rows());
    }

    // Jump backward the length of a row n times
    inline col_major_iterator<T> &minus_vec (const int &n)
    {
      return (*this -= (n * this->matrix_->rows()));
    }

    // Jump to the beginnin of the next vector
    inline col_major_iterator<T> &next_vec ()
    {
      int col = (int) (this->current_ - 
        ((int) (this->current_ / this->matrix_->cols()) * this->matrix_->cols()));
      if (col + 1 < this->matrix_->cols())
          this->current_ = col + 1;
      else
        this->current_ = this->matrix_->size();

      return *this;
    }

    // Jump to the beginning of the nth next vector
    inline col_major_iterator<T> &next_vec (const int &n)
    {
      int col = (int) (this->current_ - 
        ((int) (this->current_ / this->matrix_->cols()) * this->matrix_->cols()));
      if (col + n >= this->matrix_->cols())
        this->current_ = this->matrix_->size();
      else if (col + n <= 0)
        this->current_ = 0;
      else
        this->current_ = col + n;

      return *this;
    }
      
    // Jump to the beginnin of the previous vector
    inline col_major_iterator<T> &prev_vec ()
    {
      int col = (int) (this->current_ - 
        ((int) (this->current_ / this->matrix_->cols()) * this->matrix_->cols()));
      if (col - 1 > 0)
          this->current_ = col - 1;
      else
        this->current_ = 0;

      return *this;
    }

    // Jump to the beginning of the nth previous vector
    inline col_major_iterator<T> &prev_vec (const int &n)
    {
      return (this->next_vec(-n));
    }
    
    friend bool operator== (const col_major_iterator<T> &a,
                            const col_major_iterator<T> &b)
    {
      if (a.current_ == b.current_ && a.matrix_ == b.matrix_)
        return true;

      return false;
    }

    friend bool operator<(const col_major_iterator &a,
                          const col_major_iterator &b)  
    {
      if (a.matrix_ != b.matrix_)
        throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, 
            "< Comparison on iterators to different matrices");

      
      int cm, bcm;
      if (a.current_ == a.matrix_->size())
        cm = a.current_;
      else {
        int row = (int) (a.current_ / a.matrix_->cols());
        int col = a.current_ - (row * a.matrix_->cols());
        cm = col * a.matrix_->rows() + row;
      }

      if (b.current_ == a.matrix_->size())
        bcm = b.current_;
      else {
        int brow = (int) (b.current_ / a.matrix_->cols());
        int bcol = b.current_ - (brow * a.matrix_->cols());
        bcm = bcol * a.matrix_->rows() + brow;
      }

      if (cm < bcm)
        return true;

      return false;
    }
      
  };

  template <class T>
  class const_col_major_iterator : public const_matrix_iterator<T>
  {  
  public:

    /**** Constructors ****/

    const_col_major_iterator ()
      : const_matrix_iterator<T> ()
    {
    }

    explicit const_col_major_iterator (const Matrix<T> &m)
      : const_matrix_iterator<T> (m)
    {
    }

    const_col_major_iterator (const const_col_major_iterator<T> &cmi)
      : const_matrix_iterator<T> (cmi)
    {
    }

    const_col_major_iterator (const col_major_iterator<T> &cmi)
      : const_matrix_iterator<T> (cmi)
    {
    }

    virtual ~const_col_major_iterator ()
    {
    }

    /**** Forward Iterator Facilities ****/

    // Assignment operator
    inline const_col_major_iterator<T> &operator= (const
        const_col_major_iterator &cmi)
    {
      this->matrix_ = cmi.matrix_;
      this->current_ = cmi.current_;

      return *this;
    }
    
    // Step forward, return new position
    inline const_col_major_iterator<T> &operator++ ()
    {
      if (this->current_ >= this->matrix_->cols() * (this->matrix_->rows() - 1)) {
        if (this->current_ >= this->matrix_->size() - 1)
          this->current_ = this->matrix_->size();
        else
          this->current_ = (this->current_ + 1) -
            (this->matrix_->rows() - 1) * this->matrix_->cols();
      } else
        this->current_ += this->matrix_->cols();

      return *this;
    }

    // Step forward, return old position
    inline const_col_major_iterator<T> operator++ (int)
    {
      col_major_iterator<T> temp = *this;
      ++(*this);
      return temp;
    }

    /**** BiDirectional Iterator Facilities ****/

    // Step back, return new position
    inline const_col_major_iterator<T> &operator-- ()
    {
      if (this->current_ > 0) {
        if (this->current_ == this->matrix_->size())
          --this->current_;
        else if (this->current_ < this->matrix_->cols()) {
          this->current_ = (this->current_ - 1) + 
            (this->matrix_->rows() - 1) * this->matrix_->cols();
        } else
          this->current_ -= this->matrix_->cols();
      }
      
      return *this;
    }

    // Step back, return old position
    inline const_col_major_iterator<T> operator-- (int)
    {
      const_col_major_iterator temp = *this;
      --(*this);
      return temp;
    }

    /**** Random Access Iterator Facilities ****/
    
    // Provide access to the [nth] element XXX int?
    inline const T &operator[] (const int &n) const
    {
      int col = (int) (n / this->matrix_->rows());
      int row = n - (col * this->matrix_->rows());

      return this->matrix_->data_[row * this->matrix_->cols_ + col];
    }

    // Step n elements
    inline const_col_major_iterator<T> &operator+= (const int &n)
    {
      int cm;

      if (this->current_ == this->matrix_->size())
        cm = this->current_;
      else {
        int row = (int) (this->current_ / this->matrix_->cols());
        int col = this->current_ - (row * this->matrix_->cols());
        cm = col * this->matrix_->rows() + row;
      }

      cm += n;

      if (cm >= this->matrix_->size())
        this->current_ = this->matrix_->size();
      else if (cm <= 0)
        this->current_ = 0;
      else {
        int col = (int) (cm / this->matrix_->rows());
        int row = cm - (col * this->matrix_->rows());
        this->current_ = row * this->matrix_->cols() + col;
      }

      return *this;
    }

    inline const_col_major_iterator<T> &operator-= (const int &n)
    {
      return (*this += -n);
    }

    /* Difference operators (for distance) */

    inline std::ptrdiff_t operator-
      (const const_col_major_iterator<T> &cmi) const
    {
      int cm, bcm;
      if (this->current_ == this->matrix_->size())
        cm = this->current_;
      else {
        int row = (int) (this->current_ / this->matrix_->cols());
        int col = this->current_ - (row * this->matrix_->cols());
        cm = col * this->matrix_->rows() + row;
      }

      if (cmi.current_ == this->matrix_->size())
        bcm = cmi.current_;
      else {
        int brow = (int) (cmi.current_ / this->matrix_->cols());
        int bcol = cmi.current_ - (brow * this->matrix_->cols());
        bcm = bcol * this->matrix_->rows() + brow;
      }
      
      return cm - bcm;
    }

    /**** Matrix Iterator Facilities ****/

    // Jump forward the length of a row
    inline const_col_major_iterator<T> &plus_vec ()
    {
      return (*this += this->matrix_->rows());
    }

    // Jump forward the length of a row n times
    inline const_col_major_iterator<T> &plus_vec (const int &n)
    {
      return (*this += (n * this->matrix_->rows()));
    }

    // Jump backward the length of a row
    inline const_col_major_iterator<T> &minus_vec ()
    {
      return (*this -= this->matrix_->rows());
    }

    // Jump backward the length of a row n times
    inline const_col_major_iterator<T> &minus_vec (const int &n)
    {
      return (*this -= (n * this->matrix_->rows()));
    }

    // Jump to the beginnin of the next vector
    inline const_col_major_iterator<T> &next_vec ()
    {
      int col = (int) (this->current_ - 
        ((int)(this->current_ / this->matrix_->cols()) * this->matrix_->cols()));
      if (col + 1 < this->matrix_->cols())
          this->current_ = col + 1;
      else
        this->current_ = this->matrix_->size();

      return *this;
    }

    // Jump to the beginning of the nth next vector
    inline const_col_major_iterator<T> &next_vec (const int &n)
    {
      int col = (int) (this->current_ - 
        ((int)(this->current_ / this->matrix_->cols()) * this->matrix_->cols()));
      if (col + n >= this->matrix_->cols())
        this->current_ = this->matrix_->size();
      else if (col + n <= 0)
        this->current_ = 0;
      else
        this->current_ = col + n;

      return *this;
    }
      
    // Jump to the beginnin of the previous vector
    inline const_col_major_iterator<T> &prev_vec ()
    {
      int col = (int) (this->current_ - 
        ((int)(this->current_ / this->matrix_->cols()) * this->matrix_->cols()));
      if (col - 1 > 0)
          this->current_ = col - 1;
      else
        this->current_ = 0;

      return *this;
    }

    // Jump to the beginning of the nth previous vector
    inline const_col_major_iterator<T> &prev_vec (const int &n)
    {
      return (this->next_vec(-n));
    }
    
    friend bool operator== (const const_col_major_iterator<T> &a,
                            const const_col_major_iterator<T> &b)
    {
      if (a.current_ == b.current_ && a.matrix_ == b.matrix_)
        return true;

      return false;
    }

    friend bool operator<(const const_col_major_iterator &a,
                          const const_col_major_iterator &b)  
    {
      if (a.matrix_ != b.matrix_)
        throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, 
            "< Comparison on iterators to different matrices");

      
      int cm, bcm;
      if (a.current_ == a.matrix_->size())
        cm = a.current_;
      else {
        int row = (int) (a.current_ / a.matrix_->cols());
        int col = a.current_ - (row * a.matrix_->cols());
        cm = col * a.matrix_->rows() + row;
      }

      if (b.current_ == a.matrix_->size())
        bcm = b.current_;
      else {
        int brow = (int) (b.current_ / a.matrix_->cols());
        int bcol = b.current_ - (brow * a.matrix_->cols());
        bcm = bcol * a.matrix_->rows() + brow;
      }

      if (cm < bcm)
        return true;

      return false;
    }
      
  };
  
  /* reverse iterator adapters */
  template <class T>
  class reverse_row_major_iterator : public row_major_iterator<T>
  {
  public:
    
    /**** Constructors ****/
    reverse_row_major_iterator ()
      : row_major_iterator<T> ()
    {
    }
    explicit reverse_row_major_iterator (Matrix<T> &m)
      : row_major_iterator<T> (m)
    {
    }

    reverse_row_major_iterator
      (const row_major_iterator<T> &rmi)
      :  row_major_iterator<T> (rmi)
    {
    }

    virtual ~reverse_row_major_iterator ()
    {
    }

    /* Let people get a row_major_iterator by the standard */
    inline row_major_iterator<T> base () const
    {
      row_major_iterator<T> temp = *this;
      if (this->current_ == this->matrix_->size())
        temp -= this->matrix_->size();
      else if (this->current_ == 0)
        temp += this->matrix_->size();
      else {
        temp += this->matrix_->size();
        temp -= this->current_;
      }

      return temp;
    }

    /* Override these to get correct r-iter behavior */
    inline T &operator* () const
    {
      if (this->current_ == this->matrix_->size())
        return this->matrix_->data_[this->matrix_->size()];
          
      return this->matrix_->data_[this->matrix_->size() - this->current_ - 1];
    }

    inline T *operator-> () const
    {
      if (this->current_ == this->matrix_->size())
        return &(this->matrix_->data_[this->matrix_->size()]);
      
      return &(this->matrix_->data_[this->matrix_->size() - this->current_ - 1]);
    }
    
    /* These need overriding cause of the whole off-by-one issue
     * with reverse iterators
     */
    int get_row() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->cols());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      
      return (int) (cur / this->matrix_->cols());
    }

    int get_col() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->rows());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      int row = (int) (cur / this->matrix_->cols());
      return (cur - (row * this->matrix_->cols()));
    }

    int get_index () const
    {
      if (this->current_ == 0)
        return this->matrix_->size() - 1;
      
      return this->matrix_->size() - this->current_ - 1;
    }
  };

  template <class T>
  class const_reverse_row_major_iterator :
    public const_row_major_iterator<T>
  {
  public:
    
    /**** Constructors ****/
    const_reverse_row_major_iterator ()
      : const_row_major_iterator<T> ()
    {
    }
    explicit const_reverse_row_major_iterator (const Matrix<T> &m)
      : const_row_major_iterator<T> (m)
    {
    }

    const_reverse_row_major_iterator
      (const const_row_major_iterator<T> &rmi)
      :  const_row_major_iterator<T> (rmi)
    {
    }

    virtual ~const_reverse_row_major_iterator ()
    {
    }

    /* Let people get a row_major_iterator by the standard */
    inline const_row_major_iterator<T> base () const
    {
      const_row_major_iterator<T> temp = *this;
      if (this->current_ == this->matrix_->size())
        temp -= this->matrix_->size();
      else if (this->current_ == 0)
        temp += this->matrix_->size();
      else {
        temp += this->matrix_->size();
        temp -= this->current_;
      }

      return temp;
    }

    /* Override these to get correct r-iter behavior */
    inline T &operator* () const
    {
      if (this->current_ == this->matrix_->size())
        return this->matrix_->data_[this->matrix_->size()];
          
      return this->matrix_->data_[this->matrix_->size() - this->current_ - 1];
    }

    inline T *operator-> () const
    {
      if (this->current_ == this->matrix_->size())
        return &(this->matrix_->data_[this->matrix_->size()]);
      
      return &(this->matrix_->data_[this->matrix_->size() - this->current_ - 1]);
    }
    
    /* These need overriding cause of the whole off-by-one issue
     * with reverse iterators
     */
    int get_row() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->cols());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      
      return (int) (cur / this->matrix_->cols());
    }

    int get_col() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->rows());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      int row = (int) (cur / this->matrix_->cols());
      return (cur - (row * this->matrix_->cols()));
    }

    int get_index () const
    {
      if (this->current_ == 0)
        return this->matrix_->size() - 1;
      
      return this->matrix_->size() - this->current_ - 1;
    }
  };
  
  template <class T>
  class reverse_col_major_iterator : public col_major_iterator<T>
  {
  public:
    
    /**** Constructors ****/
    reverse_col_major_iterator ()
      : col_major_iterator<T> ()
    {
    }
    explicit reverse_col_major_iterator (Matrix<T> &m)
      : col_major_iterator<T> (m)
    {
    }

    reverse_col_major_iterator
      (const col_major_iterator<T> &rmi)
      :  col_major_iterator<T> (rmi)
    {
    }

    virtual ~reverse_col_major_iterator ()
    {
    }

    /* Let people get a row_major_iterator by the standard */
    inline col_major_iterator<T> base () const
    {
      col_major_iterator<T> temp = *this;
      --temp;
      if (this->current_ == this->matrix_->size())
        temp -= this->matrix_->size();
      else if (this->current_ == 0)
        temp += this->matrix_->size();
      else {
        temp += this->matrix_->size();
        temp -= this->current_;
      }

      return temp;
    }

    /* Override these to get correct r-iter behavior */
    inline T &operator* () const
    { 
      if (this->current_ == this->matrix_->size())
        return this->matrix_->data_[this->matrix_->size()];
      
      return this->matrix_->data_[this->matrix_->size() - this->current_ - 1];
    }

    inline T *operator-> () const
    {
      if (this->current_ == this->matrix_->size())
        return &(this->matrix_->data_[this->matrix_->size()]);
      
      return &(this->matrix_->data_[this->matrix_->size() - this->current_ - 1]);
    }
    
    /* These need overriding cause of the whole off-by-one issue
     * with reverse iterators
     */
    int get_row() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->cols());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      
      return (int) (cur / this->matrix_->cols());
    }

    int get_col() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->rows());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      int row = (int) (cur / this->matrix_->cols());
      return (cur - (row * this->matrix_->cols()));
    }

    int get_index () const
    {
      if (this->current_ == 0)
        return this->matrix_->size() - 1;
      
      return this->matrix_->size() - this->current_ - 1;
    }
  };
  
  template <class T>
  class const_reverse_col_major_iterator :
    public const_col_major_iterator<T>
  {
  public:
    
    /**** Constructors ****/
    const_reverse_col_major_iterator ()
      : const_col_major_iterator<T> ()
    {
    }
    explicit const_reverse_col_major_iterator (const Matrix<T> &m)
      : const_col_major_iterator<T> (m)
    {
    }

    const_reverse_col_major_iterator
      (const const_col_major_iterator<T> &rmi)
      :  const_col_major_iterator<T> (rmi)
    {
    }

    virtual ~const_reverse_col_major_iterator ()
    {
    }

    /* Let people get a row_major_iterator by the standard */
    inline const_col_major_iterator<T> base () const
    {
      const_col_major_iterator<T> temp = *this;
      --temp;
      if (this->current_ == this->matrix_->size())
        temp -= this->matrix_->size();
      else if (this->current_ == 0)
        temp += this->matrix_->size();
      else {
        temp += this->matrix_->size();
        temp -= this->current_;
      }

      return temp;
    }

    inline T &operator* () const
    { 
      if (this->current_ == this->matrix_->size())
        return this->matrix_->data_[this->matrix_->size()];
      
      return this->matrix_->data_[this->matrix_->size() - this->current_ - 1];
    }

    inline T *operator-> () const
    {
      if (this->current_ == this->matrix_->size())
        return &(this->matrix_->data_[this->matrix_->size()]);
      
      return &(this->matrix_->data_[this->matrix_->size() - this->current_ - 1]);
    }
    
    /* These need overriding cause of the whole off-by-one issue
     * with reverse iterators
     */
    int get_row() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->cols());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      
      return (int) (cur / this->matrix_->cols());
    }

    int get_col() const
    {
      if (this->current_ == 0)
        return (int) ((this->matrix_->size() - 1) / this->matrix_->rows());
      else if (this->current_ == this->matrix_->size())
        return 0;

      int cur = this->matrix_->size() - this->current_ - 1;
      int row = (int) (cur / this->matrix_->cols());
      return (cur - (row * this->matrix_->cols()));
    }

    int get_index () const
    {
      if (this->current_ == 0)
        return this->matrix_->size() - 1;
      
      return this->matrix_->size() - this->current_ - 1;
    }
  };

  template <class T>
  inline bool operator!= (const const_row_major_iterator<T> &a,
                          const const_row_major_iterator<T> &b)
  {
    return ! (a == b);
  }

  template <class T>
  inline bool operator!= (const const_col_major_iterator<T> &a,
                          const const_col_major_iterator<T> &b)
  {
    return ! (a == b);
  }

  template <class T>
  inline bool operator!= (const row_major_iterator<T> &a,
                          const row_major_iterator<T> &b)
  {
    return ! (a == b);
  }

  template <class T>
  inline bool operator!= (const col_major_iterator<T> &a,
                          const col_major_iterator<T> &b)
  {
    return ! (a == b);
  }

  template <class T>
  inline bool operator>(const const_row_major_iterator<T> &a,
                        const const_row_major_iterator<T> &b)
  {
    return ! (a < b);
  }

  template <class T>
  inline bool operator>(const const_col_major_iterator<T> &a,
                        const const_col_major_iterator<T> &b)
  {
    return ! (a < b);
  }

  template <class T>
  inline bool operator>(const row_major_iterator<T> &a,
                        const row_major_iterator<T> &b)
  {
    return ! (a < b);
  }

  template <class T>
  inline bool operator>(const col_major_iterator<T> &a,
                        const col_major_iterator<T> &b)
  {
    return ! (a < b);
  }

  template <class T>
  inline bool operator<= (const const_row_major_iterator<T> &a,
                          const const_row_major_iterator<T> &b)
  {
    return (a < b || a == b);
  }

  template <class T>
  inline bool operator<= (const const_col_major_iterator<T> &a,
                          const const_col_major_iterator<T> &b)
  {
    return (a < b || a == b);
  }

  template <class T>
  inline bool operator<= (const row_major_iterator<T> &a,
                          const row_major_iterator<T> &b)
  {
    return (a < b || a == b);
  }

  template <class T>
  inline bool operator<= (const col_major_iterator<T> &a,
                          const col_major_iterator<T> &b)
  {
    return (a < b || a == b);
  }

  template <class T>
  inline bool operator>= (const const_row_major_iterator<T> &a,
                          const const_row_major_iterator<T> &b)
  {
    return (a > b || a == b);
  }

  template <class T>
  inline bool operator>= (const const_col_major_iterator<T> &a,
                          const const_col_major_iterator<T> &b)
  {
    return (a > b || a == b);
  }
  
  template <class T>
  inline bool operator>= (const row_major_iterator<T> &a,
                          const row_major_iterator<T> &b)
  {
    return (a > b || a == b);
  }

  template <class T>
  inline bool operator>= (const col_major_iterator<T> &a,
                          const col_major_iterator<T> &b)
  {
    return (a > b || a == b);
  }
  
  /* Non-member arithmetic operators for various iterators */
  template <class T>
  inline row_major_iterator<T> operator+ (row_major_iterator<T> rmi,
                                          const int &n)
  {
    rmi += n;
    return rmi;
  }
  
  template <class T>
  inline row_major_iterator<T> operator+ (const int &n,
                                          row_major_iterator<T> rmi)
  {
    rmi += n;
    return rmi;
  }
  
  template <class T>
  inline row_major_iterator<T> operator- (row_major_iterator<T> rmi,
                                          const int &n)
  {
    rmi -= n;
    return rmi;
  }
  
  template <class T>
  inline col_major_iterator<T> operator+ (col_major_iterator<T> cmi,
                                          const int &n)
  {
    cmi += n;
    return cmi;
  }
  
  template <class T>
  inline col_major_iterator<T> operator+ (const int &n,
                                          col_major_iterator<T> cmi)
  {
    cmi += n;
    return cmi;
  }
  
  template <class T>
  inline col_major_iterator<T> operator- (col_major_iterator<T> cmi,
                                          const int &n)
  {
    cmi -= n;
    return cmi;
  }
  
  template <class T>
  inline const_row_major_iterator<T> operator+
    (const_row_major_iterator<T> rmi, const int &n)
  {
    rmi += n;
    return rmi;
  }
  
  template <class T>
  inline const_row_major_iterator<T> operator+ (const int &n,
                                        const_row_major_iterator<T> rmi)
  {
    rmi += n;
    return rmi;
  }
  
  template <class T>
  inline const_row_major_iterator<T> operator-
    (const_row_major_iterator<T> rmi, const int &n)
  {
    rmi -= n;
    return rmi;
  }
  
  template <class T>
  inline const_col_major_iterator<T> operator+
    (const_col_major_iterator<T> cmi, const int &n)
  {
    cmi += n;
    return cmi;
  }
  
  template <class T>
  inline const col_major_iterator<T> operator+ (const int &n,
                                      const_col_major_iterator<T> cmi)
  {
    cmi += n;
    return cmi;
  }
  
  template <class T>
  inline const_col_major_iterator<T> operator- 
    (const_col_major_iterator<T> cmi, const int &n)
  {
    cmi -= n;
    return cmi;
  }

}   // end namespace SCYTHE

#endif
