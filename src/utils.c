#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "utils.h"

double *convert_to_qhdata(const float *rm_buff, size_t n_data) {
  double *qhdata= (double*)malloc(n_data*sizeof(double)*4);
  int i,j;
  for (i=0;i<n_data;++i) {
    for (j=0;j<3;++j)
      qhdata[i*4+j]=(double)rm_buff[j*n_data+i];
  }
  return qhdata;
}

void shuffle_buff(double *data, int n_particles_in_buff, int n_shuffle) {
  int i, j, r;
  double temp[4];
  for (i=1; i<n_shuffle;++i) {
    r = rand()%(n_particles_in_buff-i)+i;
    for (j=0;j<4;++j)
      temp[j]=(data)[r*4+j];
    for (j=0;j<4;++j)
      (data)[r*4+j]=(data)[(i-1)*4+j];
    for (j=0;j<4; ++j)
      (data)[(i-1)*4+j] = temp[j];
  }
}

// perform a 3-d rotation of the points
void rotate3d( double *data, size_t n_particles, double theta[3], double center[3] ) {

    int i, j, m;
    double x, y, z;
    double Qx[  3 ][ 3 ];
    double Qy[  3 ][ 3 ];
    double Qz[  3 ][ 3 ];
    double Qt1[ 3 ][ 3 ];
    double Qt2[ 3 ][ 3 ];

    Qx[ 0 ][ 0 ] = 1.0f;            Qx[ 0 ][ 1 ] = 0.0f;            Qx[ 0 ][ 2 ] = 0.0f;
    Qx[ 1 ][ 0 ] = 0.0f;            Qx[ 1 ][ 1 ] = cos(theta[0]);   Qx[ 1 ][ 2 ] = -sin(theta[0]);
    Qx[ 2 ][ 0 ] = 0.0f;            Qx[ 2 ][ 1 ] = sin(theta[0]);   Qx[ 2 ][ 2 ] = cos(theta[0]);

    Qy[ 0 ][ 0 ] = cos(theta[1]);   Qy[ 0 ][ 1 ] = 0.0f;            Qy[ 0 ][ 2 ] = sin(theta[1]);
    Qy[ 1 ][ 0 ] = 0.0f;            Qy[ 1 ][ 1 ] = 1.0f;            Qy[ 1 ][ 2 ] = 0.0f;
    Qy[ 2 ][ 0 ] = -sin(theta[1]);  Qy[ 2 ][ 1 ] = 0.0f;            Qy[ 2 ][ 2 ] = cos(theta[1]);

    Qz[ 0 ][ 0 ] = cos(theta[2]);   Qz[ 0 ][ 1 ] = -sin(theta[2]);  Qz[ 0 ][ 2 ] = 0.0f;
    Qz[ 1 ][ 0 ] = sin(theta[2]);   Qz[ 1 ][ 1 ] = cos(theta[2]);   Qz[ 1 ][ 2 ] = 0.0f;
    Qz[ 2 ][ 0 ] = 0.0f;            Qz[ 2 ][ 1 ] = 0.0f;            Qz[ 2 ][ 2 ] = 1.0f;

    for ( i = 0; i < 3; ++i ) {
        for( j = 0; j < 3; ++j ) {
            Qt1[ i ][ j ] = 0;
            for( m = 0; m < 3; ++m )
                Qt1[ i ][ j ] += Qx[ i ][ m ] * Qy[ m ][ j ];
        }
    }

    for ( i = 0; i < 3; ++i ) {
        for( j = 0; j < 3; ++j ) {
            Qt2[ i ][ j ] = 0;
            for( m = 0; m < 3; ++m )
                Qt2[ i ][ j ] += Qt1[ i ][ m ] * Qz[ m ][ j ];
        }
    }

    for ( i = 0; i < n_particles; ++i ) {
        x = (data)[ i * 4 + 0 ] - center[0];
        y = (data)[ i * 4 + 1 ] - center[1];
        z = (data)[ i * 4 + 2 ] - center[2];
        (data)[ i * 4 + 0 ] = x * Qt2[ 0 ][ 0 ] + y * Qt2[ 0 ][ 1 ] + z * Qt2[ 0 ][ 2 ] + center[0];
        (data)[ i * 4 + 1 ] = x * Qt2[ 1 ][ 0 ] + y * Qt2[ 1 ][ 1 ] + z * Qt2[ 1 ][ 2 ] + center[1];
        (data)[ i * 4 + 2 ] = x * Qt2[ 2 ][ 0 ] + y * Qt2[ 2 ][ 1 ] + z * Qt2[ 2 ][ 2 ] + center[2];
    }
}
