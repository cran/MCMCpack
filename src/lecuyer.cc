#ifndef SCYTHE_LECUYER_CC
#define SCYTHE_LECUYER_CC

#include "lecuyer.h"

namespace scythe {
  /* Default seed definition */
  double lecuyer::nextSeed[6] = 
      {
         12345.0, 12345.0, 12345.0, 12345.0, 12345.0, 12345.0
      };
}

#endif
