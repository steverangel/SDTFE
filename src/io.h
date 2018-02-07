#ifndef _IO_FUNC_h
#define _IO_FUNC_h

#include <stdlib.h>
#include <stdint.h>

typedef enum Read_Opt{ POS, VEL, POSVEL } read_opt;

typedef struct fof_halo_properties {
   int32_t fof_halo_count;
   int64_t fof_halo_tag;
   float fof_halo_mass;
   float fof_halo_center_x;
   float fof_halo_center_y;
   float fof_halo_center_z;
   float fof_halo_mean_x;
   float fof_halo_mean_y;
   float fof_halo_mean_z;
   float fof_halo_mean_vx;
   float fof_halo_mean_vy;
   float fof_halo_mean_vz;
   float fof_halo_vel_disp;
   float redshift;
} fof_halo_properties;

typedef struct particle_data {
    float *x;
    float *y;
    float *z;
    float *vx;
    float *vy;
    float *vz;
} particle_data;

//void write_image( char *, double *, int );
int write_rho( char*, double **, int );
void print_halo_properties( fof_halo_properties *); 
int read_file( char *, fof_halo_properties *, float **);
int read_file2( char *, fof_halo_properties *, particle_data *, read_opt, int );

#endif
