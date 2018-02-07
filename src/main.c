// Author: Esteban Rangel
// Description: This is a very simple example driver for using the SDTFE library.

#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#if defined(_OPENMP) 
#include <omp.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <stdint.h>

#include "sample.h"
#include "triangulate.h"
#include "dtfe.h"

char *basename(char const *path)
{
    char *s = strrchr(path, '/');
    if (!s)
        return strdup(path);
    else
        return strdup(s + 1);
}

int main( int argc, char *argv[] ) {

    char filename[512];
    if (argc<2) {
      printf("Usage: dtfe [ path_to_file n_particles grid_dim center_x center_y center_z field_width field_depth particle_mass]\n");
      printf("Enter the file location):\n");
      scanf("%s",filename);
    }
    else
      strcpy(filename,argv[1]);

    int n_particles;
    if (argc<3) {
      printf("Enter the number of particles in the file:\n");
      scanf("%d",&n_particles);
    }
    else
      n_particles = atoi(argv[2]);

    int grid_dim;
    if (argc<4) {
      printf("Enter the grid cells per dimension (square assumed):\n");
      scanf("%d",&grid_dim);
    }
    else
      grid_dim = atoi(argv[3]);

    double field_vol_center[3];
    if (argc<7) {
      printf("Enter the center of the 3D field volume (use spaces to separate coordinates):\n");
      scanf("%lf%lf%lf", &field_vol_center[0], &field_vol_center[1], &field_vol_center[2]);
    }
    else {
      field_vol_center[0] = atof(argv[4]);
      field_vol_center[1] = atof(argv[5]);
      field_vol_center[2] = atof(argv[6]);
    }

    double box_len;
    if (argc<8) {
      printf("Enter the width of the 2D field (square assumed):\n");
      scanf("%lf", &box_len);
    }
    else
      box_len = atof(argv[7]);

    double box_len_z;
    if (argc<9) {
      printf("Enter the depth of the 3D field (set limts of integration):\n");
      scanf("%lf", &box_len_z);
    }
    else
      box_len_z = atof(argv[8]);

    float pM;
    if (argc<10) {
      printf("Enter the particle mass\n");
      scanf("%f", &pM);
    }
    else
      pM = atof(argv[9]);

    float mc_box_width;
    if (argc<11) {
      printf("Enter the width for the Monte Carlo sample volume of each pixel:\n");
      scanf("%f", &mc_box_width);
    }
    else
      mc_box_width = atof(argv[10]);

    int n_mc_samples;
    if (argc<12) {
      printf("Enter number of Monte Carlo samples for each pixel:\n");
      scanf("%d", &n_mc_samples);
    }
    else
      n_mc_samples = atoi(argv[11]);

    srand(time(NULL));
    /* 
     * In this example the input file is a binary (row-major) single precision array of positions
     * in Cartesian space. Qhull and this code compute in double precision and conversion is needed.
     * */

    // Note: The format for qhull, however the 4th location is reused for on-site density estimate.

    float *buff = (float*)malloc(n_particles*sizeof(float)*3);
    
    FILE *ptr_myfile;
    ptr_myfile=fopen(filename,"rb");
    int t=fread(buff,sizeof(float),n_particles*3,ptr_myfile);
    fclose(ptr_myfile);

    double *sample_data = (double*)malloc(n_particles*sizeof(double)*4);
    // Read into a column-major array.
    int i;
    for (i=0;i<n_particles;++i) {
      sample_data[i*4+0]=(double)buff[0*n_particles+i];
      sample_data[i*4+1]=(double)buff[1*n_particles+i];
      sample_data[i*4+2]=(double)buff[2*n_particles+i];
      sample_data[i*4+3]=0.0;
    }
    free( buff );

    // Note: allocated by the "triangulate" function as 
    // ( A, B, C, D, neighbor_ba, neighbor_ac, neighbor_dc, neighbor_bd )
    int *tetra_data;
    int num_tetra;

    printf("Creating Delaunay triangulation...\n");
    num_tetra = triangulate( sample_data, n_particles, &tetra_data, "qhull d Qz Qt Qbb" );
    printf("Done creating %d tetrahedra.\n", num_tetra);


    printf("Creating %dx%d surface density field...\n", grid_dim, grid_dim);
    printf("Center: %f, %f, %f\n", field_vol_center[0], field_vol_center[1], field_vol_center[2] );
    printf("box_len:%f box_len_z: %f\n", box_len, box_len_z);

    // density field is stored in the rho array, must be zero initialized.
    double *rho = (double*)malloc(grid_dim*grid_dim*sizeof(double));
    for (i=0;i<grid_dim*grid_dim;++i)
      rho[i]=0.0;

    compute_density( sample_data, n_particles, tetra_data, num_tetra, grid_dim, \
        box_len, box_len_z, &rho, pM, field_vol_center[0], field_vol_center[1], field_vol_center[2], \
        n_mc_samples, mc_box_width );

    free( tetra_data );
    free( sample_data );

    char output[256] = "./";
    strcat(output,basename(filename));
    strcat(output,".rho.bin");

    char image[256] = "./";
    strcat(image,basename(filename));
    strcat(image,".rho.tiff");

    write_image(image, rho, grid_dim); 
    write_rho(output, &rho, grid_dim*grid_dim); 
    exit( EXIT_SUCCESS );
}
