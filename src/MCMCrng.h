// MCMCrng.h is the header file for MCMCrng.cc. It contains 
// functions used to handle random number generator streams.
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
// ADM 7/22/04

#ifndef MCMCRNG_H
#define MCMCRNG_H

#include "rng.h"
#include "mersenne.h"
#include "lecuyer.h"

//#include <R.h>           // needed to use Rprintf()

using namespace std;
//using namespace SCYTHE;

namespace SCYTHE {

   // function that sets rng streams
   rng *MCMCpack_get_rng(const int, const int [], const int);

}// end namespace SCYTHE


#endif
