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
 * scythestat/stat.h
 *
 * Provides declarations for descriptive statistical
 * functions.
 *
 */

#ifndef SCYTHE_STAT_H
#define SCYTHE_STAT_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#else
#include "scythestat/matrix.h"
#endif

namespace SCYTHE {
  
  /* Sum - Calculate the sum of a Matrix */
  template <class T>
  T sum (const Matrix<T> &);

  /* Sumc - Calculate the sum of each column of a Matrix */
  template <class T>
  Matrix<T> sumc (const Matrix<T> &);
  
  /* Prod - Calculate the product of a Matrix */
  template <class T>
  T prod (const Matrix<T> &);

  /* Prodc - Calculate the product of each column of a Matrix */
  template <class T>
  Matrix<T> prodc (const Matrix<T> &);
  
  /* Mean - Calculate the mean of a Matrix */
  template <class T>
  T mean (const Matrix<T> &);

  /* Meanc - Calculate the mean of each column of a Matrix */
  template <class T>
  Matrix<T> meanc (const Matrix<T> &);

  /* Median - Calculate the median of a Matrix */
  template <class T>
  T median (const Matrix<T> &);

  /* Medianc - Calculate the median of each column of a Matrix */
  template <class T>
  Matrix<T> medianc (const Matrix<T> &);

  /* Mode - Calcualte the mode of the Matrix */
  template <class T>
  T mode (const Matrix<T> &);
  
  /* Modec - Calculate the mode of each column of a Matrix */
  template <class T>
  Matrix<T> modec (const Matrix<T> &);
  
  /* Skew - Calcualte the skew of the Matrix */
  template <class T>
  T skew (const Matrix<T> &);
  
  /* Skewc - Calculate the skew of each column of a Matrix */
  template <class T>
  Matrix<T> skewc (const Matrix<T> &);
  
  /* Kurtosis - Calcualte the kurtosis of the Matrix */
  template <class T>
  T kurtosis (const Matrix<T> &);
  
  /* Kurtosisc - Calculate the kurtosis of each column of a Matrix */
  template <class T>
  Matrix<T> kurtosisc (const Matrix<T> &);
  
  /* Var - Calculate the variance of a Matrix */
  template <class T>
  T var (const Matrix<T> &);

  /* Varc - Calculate the variance of each Matrix column */
  template <class T>
  Matrix<T> varc (const Matrix<T> &);

  /* Std - Calculate the std deviation of a Matrix */
  template <class T> 
  T sd (const Matrix<T> &);
  
  /* Stdc - Calculate the std deviation of each Matrix column */
  template <class T> 
  Matrix<T> stdc (const Matrix<T> &);

  /* Max - Calculates the maximum element in a Matrix */
  template <class T> 
  T max (const Matrix<T> &);

  /* Min - Calculates the minimum element in a Matrix */
  template <class T> 
  T min (const Matrix<T> &);

  /* Maxind - Finds the index of the max element */
  template <class T>
  int maxind(const Matrix<T> &);

  /* Minind - Finds the index of the min element */
  template <class T>
  int minind(const Matrix<T> &);

  /* Maxc - Calculates the maximum of each Matrix column */
  template <class T>
  Matrix<T> maxc (const Matrix<T> &);

  /* Minc - Calculates the minimum of each Matrix column */
  template <class T>
  Matrix<T> minc (const Matrix<T> &);

  /* Maxindc - Finds the index of the max of each Matrix column */
  template <class T>
  Matrix<int> maxindc(const Matrix<T> &);
  
  /* Minindc - Finds the index of the min of each Matrix column */
  template <class T>
  Matrix<int> minindc(const Matrix<T> &);

} // end namespace SCYTHE

#if defined (SCYTHE_COMPILE_DIRECT) && \
	  (defined (__GNUG__) || defined (__MWERKS__) || \
		 defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION))
#include "stat.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_STAT_H */
