#ifdef __cplusplus
extern "C" {
#endif

#ifndef _DTFE_h
#define _DTFE_h

#include <stdint.h>
#include "constants.h"
#include "bitmap.h"
#include "io.h"

void compute_density(double *particle_data, int n_particles, int *tetra_data, int n_tetra, int grid_dim, \
  double box_len, double box_len_z, double *rho, float p_mass, double center_x, double center_y, double center_z, \
  const int n_mc_samp, const double delta_sample);

void compute_3d_density(double *particle_data, int n_particles, int *tetra_data, int n_tetra, int grid_dim, \
  double box_len, double *rho, float p_mass, double center_x, double center_y, double center_z, unsigned supersampling);
#endif

#ifdef __cplusplus
}
#endif
