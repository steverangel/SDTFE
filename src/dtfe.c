#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "constants.h"
#include "triangulate.h"
#include "dtfe.h"

int dbl_comp( const void * a, const void * b ) {
  double da = *( double* )a;
  double db = *( double* )b;
  return ( da > db ) - ( da < db );
}

static inline int Sign(const double x) {
  return ( x > ZERO ? 1 : ( x < -ZERO ? -1 : 0) );
}

// Permuted inner product
static inline double vect_plucker_perm_inner_prod2( double *u, double *v ) {
  return u[2] * v[5] + v[0] * u[3] + v[1] * u[4];
}

// { CB, AC, BA, DC, BD, DA }
int rayTetraPlucker_full( double plRay2[], double *vert[],
  int *enterFace, int *leaveFace, double enterPoint[], double leavePoint[], int *isSpecialCase ) {

  double uEnter1, uEnter2, uLeave1, uLeave2;
  int ready = 0;
  int sigma[ 6 ];
  double u[ 6 ];
  double invVol;

  double tmp_pl2[ 6 ];

  *enterFace = -1;
  *leaveFace = -1;
  *isSpecialCase = 0;

  tmp_pl2[ 0 ] = vert[ B ][ 0 ] - vert[ D ][ 0 ];
  tmp_pl2[ 1 ] = vert[ B ][ 1 ] - vert[ D ][ 1 ];
  tmp_pl2[ 5 ] = vert[ B ][ 0 ] * vert[ D ][ 1 ] - vert[ B ][ 1 ] * vert[ D ][ 0 ];
  u[ BD ] = vect_plucker_perm_inner_prod2( plRay2, tmp_pl2 );
  sigma[ BD ] = Sign( u[ BD ] );

  tmp_pl2[ 0 ] = vert[ C ][ 0 ] - vert[ D ][ 0 ];
  tmp_pl2[ 1 ] = vert[ C ][ 1 ] - vert[ D ][ 1 ];
  tmp_pl2[ 5 ] = vert[ C ][ 0 ] * vert[ D ][ 1 ] - vert[ C ][ 1 ] * vert[ D ][ 0 ];
  u[ DC ] = vect_plucker_perm_inner_prod2( plRay2, tmp_pl2 );
  sigma[ DC ] = Sign( u[ DC ] );

  tmp_pl2[ 0 ] = vert[ C ][ 0 ] - vert[ B ][ 0 ];
  tmp_pl2[ 1 ] = vert[ C ][ 1 ] - vert[ B ][ 1 ];
  tmp_pl2[ 5 ] = vert[ C ][ 0 ] * vert[ B ][ 1 ] - vert[ C ][ 1 ] * vert[ B ][ 0 ];
  u[ CB ] = vect_plucker_perm_inner_prod2( plRay2, tmp_pl2 );
  sigma[ CB ] = Sign( u[ CB ] );

  // { BDC, ACD, ADB, ABC }
  // { CB, AC, BA, DC, BD, DA }
  // Face BDC
  // edge vectors: CD, BC, DB
  if ( ( sigma[ CB ] != 0 ) || ( sigma[ DC ] != 0 ) || ( sigma[ BD ] != 0 ) ) {
    invVol = 1.0f / ( -u[ CB ] + u[ DC ] - u[ BD ] );
    if ( ( *enterFace == -1 ) && ( -sigma[ CB ] >= 0 ) && ( sigma[ DC ] >= 0 ) && ( -sigma[ BD ] >= 0 ) ) {
      *enterFace = BDC;
      uEnter1 = -u[ BD ] * invVol;
      uEnter2 =  u[ DC ] * invVol;
      enterPoint[ 0 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ D ][ 0 ] + uEnter1 * vert[ C ][ 0 ] + uEnter2 * vert[ B ][ 0 ];
      enterPoint[ 1 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ D ][ 1 ] + uEnter1 * vert[ C ][ 1 ] + uEnter2 * vert[ B ][ 1 ];
      enterPoint[ 2 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ D ][ 2 ] + uEnter1 * vert[ C ][ 2 ] + uEnter2 * vert[ B ][ 2 ];
      if ( ( sigma[ CB ] == 0 ) || ( sigma[ DC ] == 0 ) || ( sigma[ BD ] == 0 ) )
        *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
        return 1;
    }
    else if ( ( *leaveFace == -1 ) && ( -sigma[ CB ] <= 0 ) && (sigma[ DC ] <= 0 ) && ( -sigma[ BD ] <= 0 ) ) {
      *leaveFace = BDC;
      uLeave1 = -u[ BD ] * invVol;
      uLeave2 =  u[ DC ] * invVol;
      leavePoint[ 0 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ D ][ 0 ] + uLeave1 * vert[ C ][ 0 ] + uLeave2 * vert[ B ][ 0 ];
      leavePoint[ 1 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ D ][ 1 ] + uLeave1 * vert[ C ][ 1 ] + uLeave2 * vert[ B ][ 1 ];
      leavePoint[ 2 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ D ][ 2 ] + uLeave1 * vert[ C ][ 2 ] + uLeave2 * vert[ B ][ 2 ];
      if ( ( sigma[ CB ] == 0 ) || ( sigma[ DC ] == 0 ) || ( sigma[ BD ] == 0 ) )
          *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
        return 1;
    }
  }
  // face ACD
  // edge vectors: DC, AD, CA
  tmp_pl2[ 0 ] = vert[ C ][ 0 ] - vert[ A ][ 0 ];
  tmp_pl2[ 1 ] = vert[ C ][ 1 ] - vert[ A ][ 1 ];
  tmp_pl2[ 5 ] = vert[ C ][ 0 ] * vert[ A ][ 1 ] - vert[ C ][ 1 ] * vert[ A ][ 0 ];
  u[ AC ] = vect_plucker_perm_inner_prod2( plRay2, tmp_pl2 );
  sigma[ AC ] = Sign( u[ AC ] );

  tmp_pl2[ 0 ] = vert[ D ][ 0 ] - vert[ A ][ 0 ];
  tmp_pl2[ 1 ] = vert[ D ][ 1 ] - vert[ A ][ 1 ];
  tmp_pl2[ 5 ] = vert[ D ][ 0 ] * vert[ A ][ 1 ] - vert[ D ][ 1 ] * vert[ A ][ 0 ];
  u[ DA ] = vect_plucker_perm_inner_prod2( plRay2, tmp_pl2 );
  sigma[ DA ] = Sign( u[ DA ] );

  if ( ( sigma[ AC ] != 0 ) || ( sigma[ DC ] != 0 ) || ( sigma[ DA ] != 0 ) ) {
    invVol = 1.0f / ( u[ AC ] - u[ DC ] - u[ DA ] );
    if ( ( *enterFace == -1 ) && ( sigma[ AC ] >= 0 ) && ( -sigma[ DC ] >= 0 ) && ( -sigma[ DA ] >= 0 ) ) {
      *enterFace = ACD;
      uEnter1 =  u[ AC ] * invVol;
      uEnter2 = -u[ DC ] * invVol;
      enterPoint[ 0 ] = ( 1 - uEnter1 - uEnter2) * vert[ C ][ 0 ] + uEnter1 * vert[ D ][ 0 ] + uEnter2 * vert[ A ][ 0 ];
      enterPoint[ 1 ] = ( 1 - uEnter1 - uEnter2) * vert[ C ][ 1 ] + uEnter1 * vert[ D ][ 1 ] + uEnter2 * vert[ A ][ 1 ];
      enterPoint[ 2 ] = ( 1 - uEnter1 - uEnter2) * vert[ C ][ 2 ] + uEnter1 * vert[ D ][ 2 ] + uEnter2 * vert[ A ][ 2 ];
      if ( ( sigma[ AC ] == 0 ) || ( sigma[ DC ] == 0 ) || ( sigma[ DA ] == 0 ) )
        *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
        return 1;
    }
    else if ( ( *leaveFace == -1 ) && ( sigma[AC] <= 0 ) && ( -sigma[DC] <= 0 ) && ( -sigma[DA] <= 0 ) ) {
      *leaveFace = ACD;
      uLeave1 =  u[ AC ] * invVol;
      uLeave2 = -u[ DC ] * invVol;
      leavePoint[ 0 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ C ][ 0 ] + uLeave1 * vert[ D ][ 0 ] + uLeave2 * vert[ A ][ 0 ];
      leavePoint[ 1 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ C ][ 1 ] + uLeave1 * vert[ D ][ 1 ] + uLeave2 * vert[ A ][ 1 ];
      leavePoint[ 2 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ C ][ 2 ] + uLeave1 * vert[ D ][ 2 ] + uLeave2 * vert[ A ][ 2 ];
      if ( ( sigma[ AC ] == 0 ) || ( sigma[ DC ] == 0 ) || ( sigma[ DA ] == 0 ) )
        *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
        return 1;
    }
  }
  // face ADB
  // edges: BD, AB, DA
  tmp_pl2[ 0 ] = vert[ B ][ 0 ] - vert[ A ][ 0 ];
  tmp_pl2[ 1 ] = vert[ B ][ 1 ] - vert[ A ][ 1 ];
  tmp_pl2[ 5 ] = vert[ B ][ 0 ] * vert[ A ][ 1 ] - vert[ B ][ 1 ] * vert[ A ][ 0 ];
  u[ BA ] = vect_plucker_perm_inner_prod2( plRay2, tmp_pl2 );
  sigma[ BA ] = Sign( u[ BA ]);

  if ( ( sigma[ BA ] != 0 ) || ( sigma[ BD ] != 0 ) || ( sigma[ DA ] != 0 ) ) {
    invVol = 1.0f / ( u[ DA ] + u[ BD ] - u[ BA ] );
    if ( ( *enterFace == -1 ) && ( -sigma[ BA ] >= 0 ) && (sigma[ BD ] >= 0 ) && ( sigma[ DA ] >= 0 ) ) {
      *enterFace = ADB;
      uEnter1 =  u[ BD ] * invVol;
      uEnter2 = -u[ BA ] * invVol;
      enterPoint[ 0 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ B ][ 0 ] + uEnter1 * vert[ A ][ 0 ] + uEnter2 * vert[ D ][ 0 ];
      enterPoint[ 1 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ B ][ 1 ] + uEnter1 * vert[ A ][ 1 ] + uEnter2 * vert[ D ][ 1 ];
      enterPoint[ 2 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ B ][ 2 ] + uEnter1 * vert[ A ][ 2 ] + uEnter2 * vert[ D ][ 2 ];
      if ( ( sigma[ BA ] == 0 ) || ( sigma[ BD ] == 0 ) || ( sigma[ DA ] == 0 ) )
        *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
        return 1;
    }
    else if ( ( *leaveFace == -1 ) && ( -sigma[ BA ] <= 0 ) && (sigma[ BD ] <= 0 ) && ( sigma[ DA ] <= 0 ) ) {
      *leaveFace = ADB;
      uLeave1 =  u[ BD ] * invVol;
      uLeave2 = -u[ BA ] * invVol;
      leavePoint[ 0 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ B ][ 0 ] + uLeave1 * vert[ A ][ 0 ] + uLeave2 * vert[ D ][ 0 ];
      leavePoint[ 1 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ B ][ 1 ] + uLeave1 * vert[ A ][ 1 ] + uLeave2 * vert[ D ][ 1 ];
      leavePoint[ 2 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ B ][ 2 ] + uLeave1 * vert[ A ][ 2 ] + uLeave2 * vert[ D ][ 2 ];
      if ( ( sigma[ BA ] == 0 ) || ( sigma[ BD ] == 0 ) || ( sigma[ DA ] == 0 ) )
        *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
        return 1;
    }
  }
  // face: ABC
  // edge: CB, AC, BA
  if ( ( sigma[ CB ] != 0 ) || ( sigma[ AC ] != 0 ) || ( sigma[ BA ] != 0 ) ) {
    invVol = 1.0f / ( u[ BA ] + u[ CB ] - u[ AC ] );
    if ( ( *enterFace == -1 ) && ( sigma[ CB ] >= 0 ) && ( -sigma[ AC ] >= 0 ) && ( sigma[ BA ] >= 0 ) ) {
      *enterFace = ABC;
      uEnter1 = -u[ AC ] * invVol;
      uEnter2 =  u[ BA ] * invVol;
      enterPoint[ 0 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ A ][ 0 ] + uEnter1 * vert[ B ][ 0 ] + uEnter2 * vert[ C ][ 0 ];
      enterPoint[ 1 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ A ][ 1 ] + uEnter1 * vert[ B ][ 1 ] + uEnter2 * vert[ C ][ 1 ];
      enterPoint[ 2 ] = ( 1 - uEnter1 - uEnter2 ) * vert[ A ][ 2 ] + uEnter1 * vert[ B ][ 2 ] + uEnter2 * vert[ C ][ 2 ];
      if ( ( sigma[ CB ] == 0 ) || ( sigma[ AC ] == 0 ) || ( sigma[ BA ] == 0 ) )
        *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
          return 1;
    }
    else if ( ( *leaveFace == -1 ) && ( sigma[ CB ] <= 0 ) && ( -sigma[ AC ] <= 0 ) && ( sigma[ BA ] <= 0 ) ) {
      *leaveFace = ABC;
      uLeave1 = -u[ AC ] * invVol;
      uLeave2 =  u[ BA ] * invVol;
      leavePoint[ 0 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ A ][ 0 ] + uLeave1 * vert[ B ][ 0 ] + uLeave2 * vert[ C ][ 0 ];
      leavePoint[ 1 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ A ][ 1 ] + uLeave1 * vert[ B ][ 1 ] + uLeave2 * vert[ C ][ 1 ];
      leavePoint[ 2 ] = ( 1 - uLeave1 - uLeave2 ) * vert[ A ][ 2 ] + uLeave1 * vert[ B ][ 2 ] + uLeave2 * vert[ C ][ 2 ];
      if ( ( sigma[ CB ] == 0 ) || ( sigma[ AC ] == 0 ) || ( sigma[ BA ] == 0 ) )
        *isSpecialCase = 1;
      if ( ready == 0 )
        ready = 1;
      else
        return 1;
    }
  }
  return 0;
}

// nudge the origin and destination points that define a ray such that it continues to
// intersect the tetrahedron it originally intersected.
static void perturb_ray( double orig[], double dest[], double **vert, double max_dist) {

  int i;
  // randomly choose a corner of the tetrahedron
  int r = rand() % 4; // interval [0, 4)

  // compute the ( x, y ) projection of ( vert[r] - origin )
  double displacement[ 2 ];
  for ( i = 0; i < 2; ++i )
    displacement[ i ] = vert[ r ][ i ] - orig[ i ];

  // if the magnitude of the displacement is larger than the maximum distance
  // scale the displacement vector.
  double mag = sqrt( displacement[ 0 ] * displacement[ 0 ] + displacement[ 1 ] * displacement[ 1 ] );
  if ( mag > max_dist ) {
    for ( i = 0; i < 2; ++i )
      displacement[ i ] = displacement[ i ] / mag;
    for ( i = 0; i < 2; ++i )
      displacement[ i ] = displacement[ i ] * max_dist;
  }

  // perturb the ray points
  for ( i = 0; i < 2; ++i ) {
    orig[ i ] = orig[ i ] + displacement[ i ];
    dest[ i ] = dest[ i ] + displacement[ i ];
  }
}

double shoot_3d_ray( int th_idx, int *tetra_data, double *nabla, double *particle_data, double x, double y,
    double offset_z, double box_len_z, double delta_xy, int grid_dim, int idx, int *tetra_crossed ){

  double rho = 0;
  double z;
  double *vert[4];
  double P_enter[3];
  double P_leave[3];
  int enterFace;
  int leaveFace;
  int is_intersection;
  int is_special_case;
  int perturb_counter;
  double orig[3];
  double dest[3];
  int tmp;
  int i;
  double dir[3];
  double plRay[ 6 ];
  plRay[0] = 0.0f;
  plRay[1] = 0.0f;
  plRay[5] = 0.0f;
  int dummy[8];
  int *tetrahedron;
  int last_th_idx;

  const double max_dist = delta_xy / 2.0;
  const double dist = max_dist / PERTURB_MAX;

  (*tetra_crossed)=0;

  if ( th_idx != -1 ) {

    assert( th_idx >= -1 );

    orig[ 0 ] = x;
    orig[ 1 ] = y;
    orig[ 2 ] = offset_z;

    dest[ 0 ] = x;
    dest[ 1 ] = y;
    dest[ 2 ] = offset_z + box_len_z;

    leaveFace = 0;
    dummy[ 4 + leaveFace ] = th_idx;
    tetrahedron = dummy;
    last_th_idx = -1;

    // dir, dest - orig
    plRay[ 2 ] = box_len_z;
    // cross, dir x orig
    plRay[ 3 ] = -box_len_z * orig[ 1 ];
    plRay[ 4 ] =  box_len_z * orig[ 0 ];

    while( tetrahedron[ 4 + leaveFace ] != -1 ) {
      // detect a turn around, and, turn around
      if ( tetrahedron[ 4 + leaveFace ] == last_th_idx ) {
        tmp = enterFace;
        enterFace = leaveFace;
        leaveFace = tmp;
        if ( tetrahedron[ 4 + leaveFace ] == -1 )
          break;
      }

      last_th_idx = th_idx;
      th_idx = tetrahedron[ 4 + leaveFace ];
      assert( th_idx >= -1 );
      tetrahedron = &tetra_data[th_idx*8];

      vert[0] = &particle_data[tetrahedron[A]*4];
      vert[1] = &particle_data[tetrahedron[B]*4];
      vert[2] = &particle_data[tetrahedron[C]*4];
      vert[3] = &particle_data[tetrahedron[D]*4];

      if ( !rayTetraPlucker_full( plRay, vert, &enterFace, &leaveFace, P_enter, P_leave, &is_special_case ) ) {
        fprintf( stderr, "Warning: rayTetraPlucker_full no intersection -- returned density 0.0 at location [ %f, %f ]\n", x, y );
        (*tetra_crossed)=0;
        return 0.0f;
      }

      assert( enterFace != leaveFace );

      perturb_counter = 0;
      while( is_special_case ) {
        // x and y coords of orig and dest are modified
        perturb_ray( orig, dest, vert, dist );

        plRay[ 3 ] = -box_len_z * orig[ 1 ];
        plRay[ 4 ] =  box_len_z * orig[ 0 ];
        plRay[ 5 ] = 0.0f;

        if ( ++perturb_counter > PERTURB_MAX ) {
          fprintf( stderr, "Warning: rayTetraPlucker_full perturb_counter max -- returned density 0.0 at location [ %f, %f ]\n", x, y );
          (*tetra_crossed)=0;
          return 0.0f;
        }

        if ( !rayTetraPlucker_full( plRay, vert, &enterFace, &leaveFace, P_enter, P_leave, &is_special_case ) ) {
          fprintf( stderr, "Warning: rayTetraPlucker_full no intersection after perturb -- returned density 0.0 at location [ %f, %f ]\n", x, y );
          (*tetra_crossed)=0;
          return 0.0f;
        }
      }
      ++(*tetra_crossed);

      // interpolate the point on the ray within the tetrahedron
      // for surface density, take the center of the ( enter, leave ) interval
      // ******************************************************************
      // check is the box starts (or ends) within the triangulation
      double enter_z = P_enter[ 2 ];
      double leave_z = P_leave[ 2 ];

      // the actual interpolation of the point
      double dz = ( leave_z - enter_z );
      double z_interp = dz / 2 + enter_z;
      if ( z_interp > offset_z && z_interp < offset_z + box_len_z ) {
        double tetra_rho = particle_data[ tetrahedron[ A ] * 4 + 3 ];
        double nab_tmp = nabla[ th_idx * 3 + 2 ];
        double fixed_part_sum = tetra_rho + \
            ( nabla[ th_idx * 3  + 0 ] * ( orig[ 0 ] - vert[ 0 ][ 0 ] ) ) + ( nabla[ th_idx * 3 + 1 ] * ( orig[ 1 ] - vert[ 0 ][ 1 ] ) ) - ( nab_tmp * vert[ 0 ][ 2 ] );

        rho += ( fixed_part_sum + nab_tmp * z_interp ) * dz;
      }
      // ******************************************************************
    }
  }
  return rho;
}

void evaluate_3d_ray_on_grid( int th_idx, int *tetra_data, double *nabla, double *particle_data, double x, double y,
    double offset_z, double box_len_z, int nz, double delta_xy, int *tetra_crossed, double* rho_vector){
  double z;
  double *vert[4];
  double P_enter[3];
  double P_leave[3];
  int enterFace;
  int leaveFace;
  int is_intersection;
  int is_special_case;
  int perturb_counter;
  double orig[3];
  double dest[3];
  int tmp;
  int i;
  double dir[3];
  double plRay[ 6 ];
  plRay[0] = 0.0f;
  plRay[1] = 0.0f;
  plRay[5] = 0.0f;
  int dummy[8];
  int *tetrahedron;
  int last_th_idx;

  double dz = box_len_z / nz;
  const double max_dist = delta_xy / 2.0;
  const double dist = max_dist / PERTURB_MAX;

  (*tetra_crossed)=0;

  if ( th_idx != -1 ) {

    assert( th_idx >= -1 );

    orig[ 0 ] = x;
    orig[ 1 ] = y;
    orig[ 2 ] = offset_z;

    dest[ 0 ] = x;
    dest[ 1 ] = y;
    dest[ 2 ] = offset_z + box_len_z;

    leaveFace = 0;
    dummy[ 4 + leaveFace ] = th_idx;
    tetrahedron = dummy;
    last_th_idx = -1;

    // dir, dest - orig
    plRay[ 2 ] = box_len_z;
    // cross, dir x orig
    plRay[ 3 ] = -box_len_z * orig[ 1 ];
    plRay[ 4 ] =  box_len_z * orig[ 0 ];

    while( tetrahedron[ 4 + leaveFace ] != -1 ) {
      // detect a turn around, and, turn around
      if ( tetrahedron[ 4 + leaveFace ] == last_th_idx ) {
        tmp = enterFace;
        enterFace = leaveFace;
        leaveFace = tmp;
        if ( tetrahedron[ 4 + leaveFace ] == -1 )
          break;
      }

      last_th_idx = th_idx;
      th_idx = tetrahedron[ 4 + leaveFace ];
      assert( th_idx >= -1 );
      tetrahedron = &tetra_data[th_idx*8];

      vert[0] = &particle_data[tetrahedron[A]*4];
      vert[1] = &particle_data[tetrahedron[B]*4];
      vert[2] = &particle_data[tetrahedron[C]*4];
      vert[3] = &particle_data[tetrahedron[D]*4];

      if ( !rayTetraPlucker_full( plRay, vert, &enterFace, &leaveFace, P_enter, P_leave, &is_special_case ) ) {
        fprintf( stderr, "Warning: rayTetraPlucker_full no intersection -- returned density 0.0 at location [ %f, %f ]\n", x, y );
        (*tetra_crossed)=0;
        return;
      }

      assert( enterFace != leaveFace );

      perturb_counter = 0;
      while( is_special_case ) {
        // x and y coords of orig and dest are modified
        perturb_ray( orig, dest, vert, dist );

        plRay[ 3 ] = -box_len_z * orig[ 1 ];
        plRay[ 4 ] =  box_len_z * orig[ 0 ];
        plRay[ 5 ] = 0.0f;

        if ( ++perturb_counter > PERTURB_MAX ) {
          fprintf( stderr, "Warning: rayTetraPlucker_full perturb_counter max -- returned density 0.0 at location [ %f, %f ]\n", x, y );
          (*tetra_crossed)=0;
          return;
        }

        if ( !rayTetraPlucker_full( plRay, vert, &enterFace, &leaveFace, P_enter, P_leave, &is_special_case ) ) {
          fprintf( stderr, "Warning: rayTetraPlucker_full no intersection after perturb -- returned density 0.0 at location [ %f, %f ]\n", x, y );
          (*tetra_crossed)=0;
          return;
        }
      }
      ++(*tetra_crossed);

      // interpolate the point on the ray within the tetrahedron
      // for surface density, take the center of the ( enter, leave ) interval
      // ******************************************************************
      // check is the box starts (or ends) within the triangulation
      double enter_z = P_enter[ 2 ];
      double leave_z = P_leave[ 2 ];

      // start and past-the-end index of gridpoints along z within tetrahedron
      int enter_iz = ceil((enter_z - offset_z)/dz);
      enter_iz = enter_iz < 0 ? 0 : enter_iz;
      int leave_iz = ceil((leave_z - offset_z)/dz);
      leave_iz = leave_iz > nz ? nz : leave_iz;

      if(leave_iz > enter_iz) {
        double tetra_rho = particle_data[ tetrahedron[ A ] * 4 + 3 ];
        double nab_tmp = nabla[ th_idx * 3 + 2 ];
        double fixed_part_sum = tetra_rho + \
            ( nabla[ th_idx * 3  + 0 ] * ( orig[ 0 ] - vert[ 0 ][ 0 ] ) ) + ( nabla[ th_idx * 3 + 1 ] * ( orig[ 1 ] - vert[ 0 ][ 1 ] ) ) - ( nab_tmp * vert[ 0 ][ 2 ] );
        // printf("  tetra=%10i z_in=%10.4e z_out=%10.4e iz_in=%3i iz_out=%3i\n", *tetra_crossed, enter_z, leave_z, enter_iz, leave_iz);
        for(int iz=enter_iz; iz < leave_iz; ++iz) {
          // z-coordinate of grid-point
          double z_interp = offset_z + dz * iz;
          rho_vector[iz] = fixed_part_sum + nab_tmp * z_interp;
        }
      }
    }
  }
}

int pt_loc_2d(int *tri_data, double *particle_data, int *start_tri, double q_x, double q_y) {
  double *vert[3];
  int next_tri = *start_tri;
  int tetra;
  do{
    vert[A]=&particle_data[tri_data[next_tri*8+A]*4];
    vert[B]=&particle_data[tri_data[next_tri*8+B]*4];
    vert[C]=&particle_data[tri_data[next_tri*8+C]*4];
    register double d_ins;
    register double d_qry;
    d_ins=(vert[C][0]-vert[A][0])*(vert[B][1]-vert[A][1]) - (vert[C][1]-vert[A][1])*(vert[B][0]-vert[A][0]);
    d_qry=(q_x       -vert[A][0])*(vert[B][1]-vert[A][1]) - (q_y-       vert[A][1])*(vert[B][0]-vert[A][0]);
    if (d_ins*d_qry<0.0) {
      next_tri = tri_data[ next_tri * 8 + ( 3 + BA ) ];
      continue;
    }
    d_ins=(vert[A][0]-vert[B][0])*(vert[C][1]-vert[B][1]) - (vert[A][1]-vert[B][1])*(vert[C][0]-vert[B][0]);
    d_qry=(q_x       -vert[B][0])*(vert[C][1]-vert[B][1]) - (q_y-       vert[B][1])*(vert[C][0]-vert[B][0]);
    if (d_ins*d_qry<0.0) {
      next_tri = tri_data[ next_tri * 8 + ( 3 + CB ) ];
      continue;
    }
    d_ins=(vert[B][0]-vert[C][0])*(vert[A][1]-vert[C][1]) - (vert[B][1]-vert[C][1])*(vert[A][0]-vert[C][0]);
    d_qry=(q_x       -vert[C][0])*(vert[A][1]-vert[C][1]) - (q_y-       vert[C][1])*(vert[A][0]-vert[C][0]);
    if (d_ins*d_qry<0.0) {
      next_tri = tri_data[ next_tri * 8 + ( 3 + AC ) ];
      continue;
    }
    break;
  }
  while(next_tri!=-1);
  *start_tri = (next_tri == -1) ? 0  : next_tri;
  tetra      = (next_tri == -1) ? -1 : tri_data[next_tri*8+6];
  return tetra;
}

void compute_density(double *particle_data, int n_particles, int *tetra_data, int n_tetra, int grid_dim, \
  double box_len, double box_len_z, double *rho, float p_mass, double center_x, double center_y, double center_z, \
  const int n_mc_samp, const double delta_sample) {

  // precompute the on-site density values and gradients
  pre_compute_vol(&particle_data, n_particles, tetra_data, n_tetra, p_mass);
  double *nabla = gradients(tetra_data, n_tetra, particle_data);

  //*******************************************************************************
  // find the forward hull facets of the 3D Delaunay triangulation convex hull
  // a forward facet is one who's outward surface normal dotted with the z unit
  // vector is positive
  // the bitmap has value 1 for a tetrahedron with forward hull face
  // here 4 is the number of tetrahedron faces
  word_t *bitmap = create_bitmap(n_tetra*4);
  int num_forward_faces = find_forward_hull(tetra_data, n_tetra, particle_data, n_particles, bitmap);

  // tri row ( A, B, C, N1, N2, N3 )
  int *tri_data = copy_fh_triangles(tetra_data, n_tetra, particle_data, n_particles, num_forward_faces, bitmap);

  // delete the bitmap
  free(bitmap);

  //*******************************************************************************
  // the pixel width
  double delta_xy = box_len / grid_dim;

  // The coordinate offsets of the grid
  double offset_x = -(box_len/2.0) + center_x;
  double offset_y = -(box_len/2.0) + center_y;
  double offset_z = -(box_len_z/2.0) + center_z;

  // if the MC sample area is bigger than the pixel, keep it, otherwise use the pixel area
  const double sample_scale  = delta_sample > delta_xy ? delta_sample : delta_xy;
  // adjust the offset for MC sampling
  const double sample_offset = delta_sample > delta_xy ? (delta_sample-delta_xy)/2.0 : 0.0;

  int num_threads = 1;
  int my_thread_id = 1;

#if defined(_OPENMP)
  #pragma omp parallel \
  default( none ) \
  shared( tri_data, tetra_data, particle_data, rho, offset_z, box_len_z, delta_xy, grid_dim, \
  num_threads, i, offset_x, offset_y, my_thread_id, nabla, tetra_crossed, tetra_crossed_srt ) \
  private(sub_dim, num_sub_elem )
{
  num_threads = omp_get_num_threads();    // for OPENMP update the number of threads
  my_thread_id = omp_get_thread_num();
#endif
  int sub_dim = 4; // must be poower of 2
  int num_sub_elem = sub_dim * sub_dim;
  int i;
#if defined(_OPENMP)
  #pragma omp for schedule( dynamic, 1 )
#endif
  for (i=0;i<grid_dim*grid_dim;i+=num_sub_elem) {
    int index = i - ( ( i % ( grid_dim * sub_dim ) ) / num_sub_elem ) * ( num_sub_elem - sub_dim );
    int j;
    for (j=0;j<sub_dim;++j) {
      int k;
      for (k=j*grid_dim;k<j*grid_dim+sub_dim;++k) {
        double x = ((index+k)%grid_dim)*delta_xy+offset_x-sample_offset;
        double y = ((index+k)/grid_dim)*delta_xy+offset_y-sample_offset;
        double rho_tmp[n_mc_samp+1]; // keep a running sum total in the last position
        int r;
        for (r=0;r<n_mc_samp+1;++r)
          rho_tmp[r]=0.0;
        int n_avg=n_mc_samp;
        for (r=0;r<n_mc_samp;++r) {
          double x_t=x+((double)rand()/RAND_MAX)*sample_scale;
          double y_t=y+((double)rand()/RAND_MAX)*sample_scale;
          int ft=0;
          int start_tri=0;
          ft=pt_loc_2d(tri_data, particle_data, &start_tri, x_t, y_t);
          if (ft==-1) {
            rho_tmp[r]=0.0;
            --n_avg;
          }
          else {
            int tmp; // the number of tetra intersected
            rho_tmp[r]+=shoot_3d_ray(ft, tetra_data, nabla, particle_data, x_t, y_t, \
                offset_z, box_len_z, delta_xy, grid_dim, index+k, &tmp);
            rho_tmp[n_mc_samp]+=rho_tmp[r];
            if (tmp==0)
              --n_avg;
          }
        }
        if (n_avg>0) {
          // you can do something more interesting here to filter outliers
          //qsort(rho_tmp,n_mc_samp,sizeof(double),dbl_comp);
          (rho)[index+k]+=rho_tmp[n_mc_samp]/n_avg;
        }
        else
          (rho)[index+k]+=0.0;
      }
    }
  }
#if defined(_OPENMP)
}
#endif
free(nabla);
free(tri_data);
}

void compute_3d_density(double *particle_data, int n_particles, int *tetra_data, int n_tetra, int grid_dim, \
  double box_len, double *rho, float p_mass, double center_x, double center_y, double center_z) {

  // precompute the on-site density values and gradients
  pre_compute_vol(&particle_data, n_particles, tetra_data, n_tetra, p_mass);
  double *nabla = gradients(tetra_data, n_tetra, particle_data);

  //*******************************************************************************
  // find the forward hull facets of the 3D Delaunay triangulation convex hull
  // a forward facet is one who's outward surface normal dotted with the z unit
  // vector is positive
  // the bitmap has value 1 for a tetrahedron with forward hull face
  // here 4 is the number of tetrahedron faces
  word_t *bitmap = create_bitmap(n_tetra*4);
  int num_forward_faces = find_forward_hull(tetra_data, n_tetra, particle_data, n_particles, bitmap);

  // tri row ( A, B, C, N1, N2, N3 )
  int *tri_data = copy_fh_triangles(tetra_data, n_tetra, particle_data, n_particles, num_forward_faces, bitmap);

  // delete the bitmap
  free(bitmap);

  //*******************************************************************************
  // the pixel width
  double delta_xy = box_len / grid_dim;

  // The coordinate offsets of the grid
  double offset_x = -(box_len/2.0) + center_x;
  double offset_y = -(box_len/2.0) + center_y;
  double offset_z = -(box_len/2.0) + center_z;

  int num_threads = 1;
  int my_thread_id = 1;

#if defined(_OPENMP)
  #pragma omp parallel \
  default( none ) \
  shared( tri_data, tetra_data, particle_data, rho, offset_z, box_len_z, delta_xy, grid_dim, \
  num_threads, i, offset_x, offset_y, my_thread_id, nabla, tetra_crossed, tetra_crossed_srt ) \
  private(sub_dim, num_sub_elem )
{
  num_threads = omp_get_num_threads();    // for OPENMP update the number of threads
  my_thread_id = omp_get_thread_num();
#endif
#if defined(_OPENMP)
  #pragma omp for schedule( dynamic, 1 )
#endif

  for(int i=0; i<grid_dim; ++i) {
    for(int j=0; j<grid_dim; ++j) {
      double x = i*delta_xy + offset_x;
      double y = j*delta_xy + offset_y;
      size_t index0 = ((i*grid_dim) + j)*grid_dim;
      int ft=0, start_tri=0;
      ft=pt_loc_2d(tri_data, particle_data, &start_tri, x, y);
      int tmp; // the number of tetra intersected
      // printf("Loop i=%3i j=%3i\n", i, j);
      evaluate_3d_ray_on_grid(
        ft, tetra_data, nabla, particle_data,
        x, y, offset_z, box_len, grid_dim, delta_xy,
        &tmp, rho + index0);
    }
  }
  #if defined(_OPENMP)
  }
  #endif
  free(nabla);
  free(tri_data);
}