/* Scythe_Util.h
 *
 * This header provides definitions and implementations for some basic
 * utilities used within the Scythe Statistical Library.
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

#ifndef SCYTHE_UTIL_H
#define SCYTHE_UTIL_H

#include <string>
#include <iterator>
#include <sstream>

namespace SCYTHE
{

  /**** A couple of useful functions that make life easier but really
   * don't have anything to do with Scythe, per se.
   ****/

  template <class T>
  inline std::string operator& (const std::string & s, const T & v)
  {
    std::ostringstream ss;
    ss << s << v;
    return ss.str ();
  }

  inline std::ostream & operator<< (std::ostream & os,
                                    const scythe_exception & se)
  {
    os << se.what ();
    return os;
  }

  template <class T>
  inline T min (const T & a, const T & b)
  {
    return b < a ? b : a;
  }

  template <class T>
  inline T max (const T & a, const T & b)
  {
    return a < b ? b : a;
  }

  template <class T>
  inline T sgn (const T & x)
  {
    if (x > 0)
      return 1;
    else if (x < 0)
      return -1;
    else
      return 0;
  }

}  // end namespace SCYTHE

#endif /* SCYTHE_ERROR_H */
