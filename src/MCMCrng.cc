// MCMCrng.h contains a function used to handle random number
// generator streams.
//
// Andrew D. Martin
// Dept. of Political Science
// Washington University in St. Louis
// admartin@wustl.edu
//
// Kevin M. Quinn
// Dept. of Government
// Harvard University
// kevin_quinn@harvard.edu
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2004 Andrew D. Martin and Kevin M. Quinn
// 
// ADM 7/23/04 (Lilja's birthday!)

#ifndef MCMCRNG_H
#define MCMCRNG_H

#include "rng.h"
#include "mersenne.h"
#include "lecuyer.h"

//#include <R.h>           // needed to use Rprintf()

using namespace std;
//using namespace SCYTHE;

namespace SCYTHE {
  
  // function to initiate random number streams
  // needs to be moved someplace else
  rng *MCMCpack_get_rng(const int lec, const int seed_array [],
			const int lecuyer_stream) {
    unsigned long u_seed_array[6];
    for(int i=0; i<6; ++i) u_seed_array[i] =
			     static_cast<unsigned long>(seed_array[i]);    
    if (lec==1) {
      lecuyer::SetPackageSeed (u_seed_array);      
      for(int i=0; i<(lecuyer_stream-1); ++i) {
	rng *retval = new lecuyer();
	delete retval; 
      }
      return new lecuyer();
    } else {
      rng *retval = new mersenne();
      dynamic_cast<mersenne *>(retval)->initialize(u_seed_array[0]);
      return retval;
    }
  }
  
}// end namespace SCYTHE

#endif
