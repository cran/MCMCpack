/* Scythe_Util.h
 *
 * This header provides definitions and implementations for some basic
 * utilities used within the Scythe Statistical Library.
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


#ifndef SCYTHE_UTIL_H
#define SCYTHE_UTIL_H

#include <string>
#include <iterator>
#include <sstream>

namespace SCYTHE {

	/**** A couple of useful functions that make life easier but really
	 * don't have anything to do with Scythe, per se.
	 ****/

	template<class T>
	inline std::string
	operator& (const std::string &s, const T &v)
	{
		std::ostringstream ss;
		ss << s << v;
		return ss.str();
	}

	inline
	std::ostream &
	operator<< (std::ostream &os, const scythe_exception &se)
	{
		os << se.what();
		return os;
	}

	template <class T>
	inline
	T
	min (const T &a, const T &b)
	{
		return b < a ? b : a;
	}
	
	template <class T>
	inline
	T
	max (const T &a, const T &b)
	{
		return a < b ? b : a;
	}

	template <class T>
	inline
	T
	sgn (const T &x) {
		if (x > 0)
			return 1;
		else if (x < 0)
			return -1;
		else
			return 0;
	}

} // end namespace SCYTHE

#endif /* SCYTHE_ERROR_H */
