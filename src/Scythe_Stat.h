/* Scythe_Stat.h
 *
 * This header provides declarations of descriptive statistical
 * functions for the Scythe Statistical Library.
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

#ifndef SCYTHE_STAT_H
#define SCYTHE_STAT_H

#include "Scythe_Matrix.h"

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
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || \
    defined (EXPLICIT_TEMPLATE_INSTANTIATION)
  // Necessary for template instantiation with some compilers.
# include "Scythe_Stat.cc"
#endif /* EXPLICIT_TEMPLATE_INSTANTIATION */


#endif /* SCYTHE_STAT_H */
