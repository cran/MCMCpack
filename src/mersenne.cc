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
 * scythestat/rng/mersenne.cc
 *
 * Provides the implementation for the mersenne class.  This is the
 * default random number generator in scythe.  See mersenne.h for
 * additional copyright information.
 *
 */

#ifndef SCYTHE_MERSENNE_CC
#define SCYTHE_MERSENNE_CC

#ifdef SCYTHE_COMPILE_DIRECT
#include "mersenne.h"
#else
#include "mersenne.h"
#endif

namespace SCYTHE {

#ifdef __MINGW32__
	/* constant vector a */
	static const unsigned long MATRIX_A = 0x9908b0dfUL;
	
	/* most significant w-r bits */
	static const unsigned long UPPER_MASK = 0x80000000UL;
	
	/* least significant r bits */
	static const unsigned long LOWER_MASK = 0x7fffffffUL;
#else
	namespace {
		/* constant vector a */
		const unsigned long MATRIX_A = 0x9908b0dfUL;
		
		/* most significant w-r bits */
		const unsigned long UPPER_MASK = 0x80000000UL;
		
		/* least significant r bits */
		const unsigned long LOWER_MASK = 0x7fffffffUL;
	}
#endif

	mersenne::mersenne ()
		:	mti (N + 1)
	{
	}

	mersenne::mersenne (const mersenne &m)
		: rng (),
			mti (m.mti)
	{
	}

	mersenne::~mersenne ()
	{
	}

  void
	mersenne::initialize (const unsigned long &s)
  {
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti] &= 0xffffffffUL;
      /* for >32 bit machines */
    }
  }

  /* generates a random number on [0,0xffffffff]-interval */
  unsigned long 
	mersenne::genrand_int32()
  {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
      int kk;

      if (mti == N+1)   /* if init_genrand() has not been called, */
        this->initialize(5489UL); /* a default initial seed is used */

      for (kk=0;kk<N-M;kk++) {
        y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for (;kk<N-1;kk++) {
        y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
      mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

      mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
  }
}

#endif /* SCYTHE_MERSENNE_CC */
