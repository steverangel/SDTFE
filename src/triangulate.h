#ifdef __cplusplus
extern "C" {
#endif

#ifndef _TRIANGULATE_h
#define _TRIANGULATE_h

#include "io.h"
#include "bitmap.h"

int triangulate( double *, int, int **, const char*);
void pre_compute_vol( double **, int, int *, int, float ); 
void pre_compute_gradients( int *, int, double *, double **);
int find_forward_hull( int *, int, double *, int, char *, word_t *);
void copy_triangles( int *, int, double *, int, int *, int, word_t *);

#endif

#ifdef __cplusplus
}
#endif
