//////////////////////////////////////////////////////////////////////////
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
// ADM 7/22/04
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////

#ifndef MCMCRNG_H
#define MCMCRNG_H

#include "mersenne.h"
#include "lecuyer.h"

/* This allows model handles to efficiently pass the appropriate rng
 * to object a model's implementation.  The first arg is the name of
 * the model implementation function.  The remaining arguments are the
 * arguments to the model implementation function.  
 *
 * The macro assumes that the function it is called in contains an int
 * pointer named uselecuyer holding a boolean indication of whether or
 * not to use the lecuyer rng.  Secondly, it assumes the function
 * contains a pointer to an array of integers called seedarray that
 * contains six random number seeds (or just one if using mersenne).
 * Finally, it assumes it contains a pointer to an integer called
 * lecuyerstream that indicates which of the lecuyer generator's
 * streams to use.
 */

#define MCMCPACK_PASSRNG2MODEL(MODEL_IMPL, ...)                           \
{                                                                         \
  unsigned long u_seed_array[6];                                          \
  for (int i = 0; i < 6; ++i)                                             \
    u_seed_array[i] = static_cast<unsigned long>(seedarray[i]);           \
                                                                          \
  if (*uselecuyer == 0) {                                                 \
    mersenne the_rng;                                                     \
    the_rng.initialize(u_seed_array[0]);                                  \
    MODEL_IMPL(the_rng, __VA_ARGS__);                                     \
  } else {                                                                \
    lecuyer::SetPackageSeed(u_seed_array);                                \
    for (int i = 0; i < (*lecuyerstream - 1); ++i)                        \
      lecuyer skip_rng;                                                   \
    lecuyer the_rng;                                                      \
    MODEL_IMPL(the_rng, __VA_ARGS__);                                     \
  }                                                                       \
}

#endif
