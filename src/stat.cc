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
 * scythestat/stat.cc
 *
 * Provides implementations of descriptive statistical
 * functions.
 *
 */

#ifndef SCYTHE_STAT_CC
#define SCYTHE_STAT_CC

#include <numeric>
#include <algorithm>
#include <cmath>

#ifdef SCYTHE_COMPILE_DIRECT
#include "stat.h"
#include "error.h"
#include "util.h"
#else
#include "scythestat/stat.h"
#include "scythestat/error.h"
#include "scythestat/util.h"
#endif

namespace SCYTHE {

  /* Calculate the sum of a Matrix */
  template <class T>
  T
  sum (const Matrix<T> &A)
  {
    return (std::accumulate(A.begin(), A.end(), (T) 0));
  }

  /* Calculate the sum of each column in a Matrix */
  template <class T>
  Matrix<T>
  sumc (const Matrix<T> &A)
  {
    Matrix<T> temp(1, A.cols(), false);
    
    for (int j = 0; j  < A.cols(); ++j)
      temp[j] = std::accumulate(A.vecc(j), A.vecc(j + 1), (T) 0);
  
    return temp;
  }
  
  /* Calculate the product of a Matrix */
  template <class T>
  T
  prod (const Matrix<T> &A)
  {
    T temp = (T) 1;
    
    for (int i = 0; i < A.size(); ++i)
      temp *= A[i];

    return temp;
  }

  /* Calculate the product of each column of a matrix */
  template <class T>
  Matrix<T>
  prodc (const Matrix<T> &A)
  {
    Matrix<T> temp(1, A.cols(), false);
    
    for (int j = 0; j  < A.cols(); ++j) {
      temp[j] = (T) 1;
      for (int i = 0; i < A.rows(); ++i)
        temp[j] *= A(i,j);
    }
    
    return temp;
  }
  
  /* Calculate the mean of a Matrix */
  template <class T>
  T
  mean (const Matrix<T> &A)
  {
    return (std::accumulate(A.begin(), A.end(), (T) 0) / A.size());
  }

  /* Calculate the mean of each column of a Matrix */
  template <class T>
  Matrix<T>
  meanc (const Matrix<T> &A)
  {
    Matrix<T> temp(1, A.cols(), false);
    
    for (int j = 0; j  < A.cols(); ++j) 
      temp[j] = std::accumulate(A.vecc(j), A.vecc(j + 1), (T) 0) / A.rows();
    
    return temp;
  }
  
  /* Calculate the median of a matrix.  Uses a sort but I'll implement
   * the randomized alg when I figure out how to generalize it to
   * even-length lists
   */
  template <class T>
  T
  median (const Matrix<T> &A)
  {
    Matrix<T> temp(A);
    int n = temp.size();

    sort(temp.begin(), temp.end());
    if (n % 2 == 0)
      return ((temp[n / 2] + temp[n / 2 - 1]) / 2);
    else
      return temp[(int) floor(n / 2)];
  }

  /* Calculate the median of each column of a matrix */
  template <class T>
  Matrix<T>
  medianc (const Matrix<T> &A)
  {
    Matrix<T> temp;
    Matrix<T> result(1, A.cols(), false);

    for (int i = 0; i < A.cols(); ++i) {
      temp = A(_, i);
      int n = temp.size();
      sort(temp.begin(), temp.end());
      if (n % 2 == 0)
        result[i] = ((temp[n / 2] +
              temp[n / 2 - 1]) / 2);
      else
        result[i] = temp[(int) floor(n / 2)];
    }

    return result;
  }

  /* Calculate the mode of a matrix */
  template <class T>
  T
  mode (const Matrix<T> &A)
  {
    Matrix<T> temp(A);
    
    sort(temp.begin(), temp.end());

    T last = temp[0];
    int cnt = 1;
    T cur_max = temp[0];
    int max_cnt = 1;
    
    for (int i = 1; i < temp.size(); ++i) {
      if (last == temp[i]) {
        ++cnt;
      } else {
        last = temp[i];
        cnt = 1;
      }
      if (cnt > max_cnt) {
        max_cnt = cnt;
        cur_max = temp[i];
      }
    }

    return cur_max;
  }

  template <class T>
  Matrix<T>
  modec (const Matrix<T> & A)
  {
    Matrix<T> temp;
    Matrix<T> result(1, A.cols(), false);

    for (int j = 0; j < A.cols(); ++j) {
      temp = A(_, j);
      T last = temp[0];
      int cnt = 1;
      T cur_max = temp[0];
      int max_cnt = 1;
      
      for (int i = 1; i < temp.size(); ++i) {
        if (last == temp[i]) {
          ++cnt;
        } else {
          last = temp[i];
          cnt = 1;
        }
        if (cnt > max_cnt) {
          max_cnt = cnt;
          cur_max = temp[i];
        }
      }
      result[j] = cur_max;
    }

    return result;
  }

  /* Calculate the skew of a Matrix */
  template <class T>
  T
  skew (const Matrix<T> &A)
  {
    T sde = sd(A);
    T mu = mean(A);
    T temp = (T) 0;

    for (int i = 0; i < A.size(); ++i) {
      temp += ::pow(A[i] - mu, 3);
    }
    temp /= A.size() * ::pow(sde, 3);

    return temp;
  }

  /* Calculate the skew of each column of a Matrix. */
  template <class T>
  Matrix<T>
  skewc (const Matrix<T> &A)
  {
    Matrix<T> sd = stdc(A);
    Matrix<T> mu = meanc(A);
    Matrix<T> temp(1, A.cols(), false);

    for (int j = 0; j < A.cols(); ++j) {
      temp[j] = 0;
      for (int i = 0; i < A.rows(); ++i) {
        temp[j] += ::pow(A(i,j) - mu[j], 3);
      }
      temp[j] /= A.cols() * ::pow(sd[j], 3);
    }

    return temp;
  }
  
  /* Calculate the kurtosis of a Matrix */
  template <class T>
  T
  kurtosis (const Matrix<T> &A)
  {
    T sde = sd(A);
    T mu = mean(A);
    T temp = (T) 0;

    for (int i = 0; i < A.size(); ++i) {
      temp += ::pow(A[i] - mu, 4);
    }
    temp /= A.size() * ::pow(sde, 4);
    temp -= 3;

    return temp;
  }
  
  /* Calculate the kurtosis of each column of a Matrix. */
  template <class T>
  Matrix<T>
  kurtosisc (const Matrix<T> &A)
  {
    Matrix<T> sd = stdc(A);
    Matrix<T> mu = meanc(A);
    Matrix<T> temp(1, A.cols(), false);

    for (int j = 0; j < A.cols(); ++j) {
      temp[j] = 0;
      for (int i = 0; i < A.rows(); ++i) {
        temp[j] += ::pow(A(i,j) - mu[j], 4);
      }
      temp[j] /= A.cols() * ::pow(sd[j], 4);
      temp[j] -= 3;
    }

    return temp;
  }
  
  /* Calculate the variance of a Matrix */
  template <class T>
  T
  var (const Matrix<T> &A)
  {
    T mu = mean(A);
    T temp = (T) 0;
    
    for (int i =0; i < A.size(); ++i)
      temp += ::pow(mu - A[i], 2) / (A.size() - 1);

    return temp;
  }

  /* Calculate the variances of each column of a Matrix. */
  template <class T>
  Matrix<T>
  varc (const Matrix<T> &A)
  {
    Matrix<T> mu = meanc (A);
    Matrix<T> temp(1, A.cols(), false);
  
    for (int j = 0; j < A.cols(); ++j) {
      temp[j] = 0;
      for (int i = 0; i < A.rows(); ++i)
        temp[j] += ::pow (mu[j] - A(i,j), 2) / (A.rows() - 1);
    }
  
    return temp;
  }
  
  /* Calculate the mean of a Matrix (not std cause of namespace std:: */
  template <class T>
  T
  sd (const Matrix<T> &A)
  {
    return ::sqrt(var(A));
  }
  
  /* Calculate the standard deviation of each column of a Matrix */
  template <class T>
  Matrix<T>
  stdc (const Matrix<T> &A)
  {
    Matrix<T> temp = varc(A);
    
    for (int i = 0; i < A.cols(); ++i)
      temp[i] = ::sqrt(temp[i]);
  
    return temp;
  }

  /* Calculates the maximum element in a Matrix */
  template <class T>
  T
  max (const Matrix<T> &A)
  {
    return *(max_element(A.begin(), A.end()));
  }
  
  /* Calculates the minimum element in a Matrix */
  template <class T>
  T
  min (const Matrix<T> &A)
  {
    return *(min_element(A.begin(), A.end()));
  }

  /* Find the index of the max element */
  template <class T>
  int
  maxind (const Matrix<T> &A)
  {
    return (max_element(A.begin(), A.end())).get_index();
  }
  
  /* Find the index of the min element */
  template <class T>
  int
  minind (const Matrix<T> &A)
  {
    return (min_element(A.begin(), A.end())).get_index();
  }
  
  /* Calculates the maximum of each Matrix column */
  template <class T>
  Matrix<T>
  maxc (const Matrix<T> &A)
  {
    Matrix<T> temp(1, A.cols(), false);
  
    for (int j = 0; j  < A.cols(); ++j) 
      temp[j] = *(max_element(A.vecc(j), A.vecc(j + 1)));
  
    return temp;
  }
  
  /* Calculates the minimum of each Matrix column */
  template <class T>
  Matrix<T>
  minc (const Matrix<T> &A)
  {
    Matrix<T> temp(1, A.cols(), false);
  
    for (int j = 0; j  < A.cols(); ++j) 
      temp[j] = *(min_element(A.vecc(j), A.vecc(j + 1)));
  
    return temp;
  }
  
  /* Finds the index of the maximum of each Matrix column */
  template <class T>
  Matrix<int>
  maxindc(const Matrix<T> &A)
  {
    Matrix<int> temp(1, A.cols(), false);
    
    for (int j = 0; j  < A.cols(); ++j) 
      temp[j] = (max_element(A.vecc(j), A.vecc(j + 1))).get_row();

    return temp;
  }
  
  /* Finds the index of the minimum of each Matrix column */
  template <class T>
  Matrix<int>
  minindc(const Matrix<T> &A)
  {
    Matrix<int> temp(1, A.cols(), false);
    
    for (int j = 0; j  < A.cols(); ++j) 
      temp[j] = (min_element(A.vecc(j), A.vecc(j + 1))).get_row();

    return temp;
  }
  
} // end namespace SCYTHE

#ifndef SCYTHE_COMPILE_DIRECT
#include "scythestat/eti/stat.t"
#endif
  
#endif /* SCYTHE_STAT_CC */
