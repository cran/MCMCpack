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
 * scythestat/util.h
 *
 * Provides definitions and implementations for some basic
 * utilities used within Scythe.
 *
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
