#ifdef __cplusplus
extern "C" {
#endif

#ifndef _UTILS_h
#define _UTILS_h

//#include "io.h"

double *convert_to_qhdata(float *, size_t); 
void   rotate3d(double **, size_t, double[3], double[3]);
void   shuffle_buff(double **, int , int); 
#endif

#ifdef __cplusplus
}
#endif
