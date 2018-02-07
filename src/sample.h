#ifdef __cplusplus
extern "C" {
#endif

#ifndef _SAMPLE_h
#define _SAMPLE_h

#include "io.h"

void tree_order( particle_data *, int **, int, int, int);
void rotate3d( double **, int, double, double, double);
void find_partition_boundaries( int, int, int **);
void sub_sample_data( particle_data *, int *, int *, int, int *, double ** );
int sub_sample_data_by_vol( particle_data *, int *, int *, int, int *, int, double **, float, float, float, double );
void shuffle_partitions( int **, int **, int, int**);
int unique_rows( particle_data *, int );

int count_data_by_vol( particle_data *, int, float, float, float, double, double); 
int dtfe_particles_in_vol( particle_data *, int, double **, float , float , float , double, double );
void shuffle_buff( double **, int , int ); 

static int partition( particle_data *, int **, int, int, int, int);
static void nth_element( particle_data *, int **, int, int, int, int);
static int longest_dim( particle_data *, int *, int, int);

#endif

#ifdef __cplusplus
}
#endif
