#ifdef __cplusplus
extern "C" {
#endif

#ifndef _TRIANGULATE_h
#define _TRIANGULATE_h

#include "io.h"
#include "bitmap.h"

int*   triangulate( double *, int, size_t *, const char*);
void   pre_compute_vol( double **, int, int *, int, float ); 
double *gradients( int *, int, double *);
int    find_forward_hull( int *, int, double *, int, word_t *);
int*   copy_fh_triangles( int *, int, double *, int, int, word_t *);

#endif

#ifdef __cplusplus
}
#endif
