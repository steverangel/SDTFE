#include <float.h>

#ifndef _CONSTANTS_
    #define _CONSTANTS_

    #define ZERO 1e4*DBL_EPSILON 
    #define DIM 3
    #define PERTURB_MAX 1000

    typedef enum { BDC, ACD, ADB, ABC } face_t;
    typedef enum { CB, AC, BA, DC, BD, DA } side_t;
    typedef enum { A, B, C, D } vertex_t;

#endif
