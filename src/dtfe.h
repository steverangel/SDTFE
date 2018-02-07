#ifdef __cplusplus
extern "C" {
#endif

#ifndef _DTFE_h
#define _DTFE_h

#include <stdint.h>
#include "constants.h"
#include "bitmap.h"
#include "io.h"

void compute_density( double *, int, int *, int, int, double, double, double **, float, double, double, double, \
  const int, const double ); 
#endif

#ifdef __cplusplus
}
#endif
