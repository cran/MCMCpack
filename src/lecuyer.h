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
 * scythestat/rng/lecuyer.h
 *
 * Provides the class definition for the L'Ecuyer random number
 * generator, a rng capable of generating many independent substreams.
 * This class extends the abstract rng class by implementing runif().
 * Based on RngStream.cpp, originally released under the following
 * license:
 *
 * Copyright:      Pierre L'Ecuyer, University of Montreal
 * Notice:         This code can be used freely for personnal,
 *                 academic, or non-commercial purposes. For commercial
 *                 purposes, please contact P. L'Ecuyer at:
 *                 lecuyer@iro.umontreal.ca
 * Date:           14 August 2001
 *
 */

#ifndef SCYTHE_LECUYER_H
#define SCYTHE_LECUYER_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "rng.h"
#else
#include "scythestat/rng.h"
#endif

namespace SCYTHE {
	
	class lecuyer : public rng
	{
		public:

			lecuyer (const char *name = "");


			static void SetPackageSeed (const unsigned long seed[6]);


			void ResetStartStream ();


			void ResetStartSubstream ();


			void ResetNextSubstream ();


			void SetAntithetic (bool);


			void IncreasedPrecis (bool);


			void SetSeed (const unsigned long seed[6]);


			void AdvanceState (long, long);


			void GetState (unsigned long seed[6]) const;


			void WriteState () const;


			void WriteStateFull () const;


			double runif ();

			/* We have to override the overloaded form of runif because
			 * overloading the no-arg runif() hides the base class
			 * definition; C++ stops looking once it finds the above.
			 */
			inline Matrix<double> runif(const int &rows, const int &cols)
			{
				return rng::runif(rows, cols);
			}

			long RandInt (long, long);

		protected:

				double Cg[6], Bg[6], Ig[6];


				bool anti, incPrec;


				std::string name;

				 
				static double nextSeed[6];

				 
				double U01 ();

				 
				double U01d ();

		};

	}

#endif /* SCYTHE_LECUYER_H */
