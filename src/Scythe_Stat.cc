/* Scythe_Stat.cc
 *
 * This file provides implementations of descriptive statistical
 * functions for the Scythe Statistical Library.
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
#ifndef SCYTHE_STAT_CC
#define SCYTHE_STAT_CC

#include <numeric>
#include <cmath>
#include "Scythe_Stat.h"
#include "Scythe_Error.h"
#include "Scythe_Util.h"

namespace SCYTHE {



  /* Calculate the sum of a Matrix */
  template <class T>
  T
  sum (const Matrix<T> &A)
  {
    return (accumulate(A.begin(), A.end(), (T) 0));
  }

  /* Calculate the sum of each column in a Matrix */
  template <class T>
  Matrix<T>
  sumc (const Matrix<T> &A)
  {
    Matrix<T> temp(1, A.cols(), false);
		
    for (int j = 0; j  < A.cols(); ++j)
      temp[j] = accumulate(A.vecc(j), A.vecc(j + 1), (T) 0);
	
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
    return (accumulate(A.begin(), A.end(), (T) 0) / A.size());
  }

  /* Calculate the mean of each column of a Matrix */
  template <class T>
  Matrix<T>
  meanc (const Matrix<T> &A)
  {
    Matrix<T> temp(1, A.cols(), false);
		
    for (int j = 0; j  < A.cols(); ++j) 
      temp[j] = accumulate(A.vecc(j), A.vecc(j + 1), (T) 0) / A.rows();
		
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
	
#endif /* SCYTHE_STAT_CC */
