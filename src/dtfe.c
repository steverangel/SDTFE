#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#if defined(_OPENMP) 
#include <omp.h>
#endif
#include "triangulate.h"
#include "qhull_a.h"
#include "constants.h"
#include "dtfe.h"

int int_comp( const void * a, const void * b ) {

    int ia = *( int* )a;
    int ib = *( int* )b;

    return ( ia > ib ) - ( ia < ib );
}

int dbl_comp( const void * a, const void * b ) {

    double da = *( double* )a;
    double db = *( double* )b;

    return ( da > db ) - ( da < db );
}

inline int Sign(const double x) {
    return ( x > ZERO ? 1 : ( x < -ZERO ? -1 : 0) );
}

// Permuted inner product 
inline double vect_plucker_perm_inner_prod2( double *u, double *v ) {
    return u[ 2 ] * v[ 5 ] + v[ 0 ] * u[ 3 ] + v[ 1 ] * u[ 4 ];
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
// time of day difference - timer for code performance evaluation
static inline long long int toddiff( struct timeval *tod1, struct timeval *tod2 )
{
    long long t1, t2;
    t1 = tod1->tv_sec * 1000000 + tod1->tv_usec;
    t2 = tod2->tv_sec * 1000000 + tod2->tv_usec;
    return t1 - t2;
}

typedef union Dbl_t{
   double d;
   long int i;
} Dbl_t;

int almostEqual( double A, double B ) {
    double abs_diff = fabs( A - B );
    if ( abs_diff <= MAX_DIFF )
        return 1;

    Dbl_t uA, uB;
    uA.d = A;
    uB.d = B;

    // check the signs
    if ( ( uA.i >> 63 ) != ( uB.i >> 63 ) )
        return 0;

    long int ulpsDiff = labs( uA.i - uB.i );
    if ( ulpsDiff <= MAX_ULP ) 
        return 1;

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
    //#pragma simd
    for ( i = 0; i < 2; ++i )
        displacement[ i ] = vert[ r ][ i ] - orig[ i ];

    // if the magnitude of the displacement is larger than the maximum distance
    // scale the displacement vector. 
    double mag = sqrt( displacement[ 0 ] * displacement[ 0 ] + displacement[ 1 ] * displacement[ 1 ] ); 
    if ( mag > max_dist ) {
        //#pragma simd
        for ( i = 0; i < 2; ++i )
            displacement[ i ] = displacement[ i ] / mag;
        //#pragma simd
        for ( i = 0; i < 2; ++i ) 
            displacement[ i ] = displacement[ i ] * max_dist;
    }

    // perturb the ray points
    //#pragma simd
    for ( i = 0; i < 2; ++i ) {
        orig[ i ] = orig[ i ] + displacement[ i ];
        dest[ i ] = dest[ i ] + displacement[ i ];
    }
}

double shoot_3d_ray( int th_idx, int *tetra_data, double *nabla, double *sample_data, double x, double y, 
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

            vert[0] = &sample_data[tetrahedron[A]*4];
            vert[1] = &sample_data[tetrahedron[B]*4];
            vert[2] = &sample_data[tetrahedron[C]*4];
            vert[3] = &sample_data[tetrahedron[D]*4];

            if ( !rayTetraPlucker_full( plRay, vert, &enterFace, &leaveFace, P_enter, P_leave, &is_special_case ) ) {
                fprintf( stderr, "Warning: rayTetraPlucker_full no intersection -- returned density 0.0 at location [ %f, %f ]\n", x, y );
                (*tetra_crossed)=0;
                return 0.0f;
            }

            assert( enterFace != leaveFace );

            perturb_counter = 0;
            while( is_special_case ) {
                //fprintf( stdout, "perturb ray\n" );
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
              double tetra_rho = sample_data[ tetrahedron[ A ] * 4 + 3 ];
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

int pt_loc_2d(int *tri_data, double *sample_data, int *start_tri, double q_x, double q_y) {
  double *vert[3];
  int next_tri = *start_tri;
  int tetra;
  do{
    vert[A]=&sample_data[tri_data[next_tri*8+A]*4];
    vert[B]=&sample_data[tri_data[next_tri*8+B]*4];
    vert[C]=&sample_data[tri_data[next_tri*8+C]*4];
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

void compute_density( double *sample_data, int num_samples, int *tetra_data, int num_tetra, int grid_dim, \
    double box_len, double box_len_z, double **rho, float p_mass, double center_x, double center_y, double center_z, \
    const int n_mc_samp, const double delta_sample ) { 

    int i, num_threads, my_thread_id;
    double x, y, offset_x, offset_y, offset_z, delta_xy;

    // the grid spacing
    delta_xy = box_len / grid_dim;

    // The coordinate offsets of the grid
    offset_x = -( box_len / 2 ) + center_x;
    offset_y = -( box_len / 2 ) + center_y;
    offset_z = -( box_len_z / 2 ) + center_z;

    // actually estimates the on-site density and stores the value in the point sample 4th array slot
    pre_compute_vol( &sample_data, num_samples, tetra_data, num_tetra, p_mass );
    
    // we compute 3 gradients for eaach tetrahedron
    double *nabla = ( double* )malloc( num_tetra * sizeof( double ) * 3 );
    pre_compute_gradients( tetra_data, num_tetra, sample_data, &nabla );

    //*******************************************************************************
    // find the forward hull facets of the 3-d delaunay triangulation convex hull
    // a forward facet is one who's outward surface normal dotted with the z unit 
    // vector is positive
    int k, j;
    int num_forward_faces = 0;

    // the bitmap has value 1 for a tetrahedron with forward hull face
    // here 4 is the number of tetrahedron faces
    word_t *bitmap = create_bitmap( num_tetra * 4 );

    num_forward_faces = find_forward_hull( tetra_data, num_tetra, sample_data, num_samples, NULL, bitmap );

    // tri row ( A, B, C, N1, N2, N3 )
    int *tri_data = (int*)malloc( num_forward_faces * 8 * sizeof( int ) );

    copy_triangles( tetra_data, num_tetra, sample_data, num_samples, tri_data, num_forward_faces, bitmap );

    // delete the bitmap 
    delete_bitmap( bitmap );

    num_threads = 1;
    my_thread_id = 1;

    int sub_dim = 1;
    int num_sub_elem = 1;

    int *tetra_crossed     = (int*)calloc( grid_dim * grid_dim, sizeof(int) );
    int *tetra_crossed_srt = (int*)calloc( grid_dim * grid_dim, sizeof(int) );

    const double sample_scale  = delta_sample > delta_xy ? delta_sample : delta_xy;
    const double sample_offset = delta_sample > delta_xy ? (delta_sample-delta_xy)/2.0 : 0.0;
#if defined(_OPENMP)
    #pragma omp parallel \
    default( none ) \
    shared( tri_data, tetra_data, sample_data, rho, offset_z, box_len_z, delta_xy, grid_dim, \
    num_threads, i, offset_x, offset_y, my_thread_id, nabla, tetra_crossed, tetra_crossed_srt ) \
    private( x, y, j, k, sub_dim, num_sub_elem ) 
{
    num_threads = omp_get_num_threads();    // for OPENMP update the number of threads
    my_thread_id = omp_get_thread_num();
#endif
    sub_dim = 4; // must be poower of 2
    num_sub_elem = sub_dim * sub_dim;
#if defined(_OPENMP)
    #pragma omp for schedule( dynamic, 1 ) 
#endif
    for ( i = 0; i < grid_dim * grid_dim; i += num_sub_elem ) {
        int index = i - ( ( i % ( grid_dim * sub_dim ) ) / num_sub_elem ) * ( num_sub_elem - sub_dim );
        for ( j = 0; j < sub_dim; ++j ) {
            for ( k = j * grid_dim; k < j * grid_dim + sub_dim; ++k ) {
              x = ( ( index + k ) % grid_dim ) * delta_xy + offset_x - sample_offset;
              y = ( ( index + k ) / grid_dim ) * delta_xy + offset_y - sample_offset;
              double rho_tmp[n_mc_samp+1]; // keep a running sum total in the last position
              int r;
              for (r=0;r<n_mc_samp+1;++r)
                rho_tmp[r]=0.0;
              int n_avg = n_mc_samp;
              for (r=0;r<n_mc_samp;++r) {
                double x_t = x + ((double)rand()/RAND_MAX)*sample_scale;
                double y_t = y + ((double)rand()/RAND_MAX)*sample_scale;
                int ft = 0;
                int start_tri = 0;
                ft =  pt_loc_2d(tri_data, sample_data, &start_tri, x_t, y_t); 
                if ( ft == -1 ) {
                  rho_tmp[r] = 0.0;
                  --n_avg;
                }
                else {
                  int tmp; // the number of tetra intersected
                  rho_tmp[r] += shoot_3d_ray( ft, tetra_data, nabla, sample_data, x_t, y_t, \
                      offset_z, box_len_z, delta_xy, grid_dim, ( index + k ), &tmp );
                  rho_tmp[n_mc_samp]+=rho_tmp[r];
                  if (tmp==0)
                    --n_avg;
                }
              }
              if (n_avg>0) {
                // you can do something more interesting here to filter outliers
                //qsort(rho_tmp,n_mc_samp,sizeof(double),dbl_comp);
                ( *rho )[ ( index + k ) ] += rho_tmp[n_mc_samp]/n_avg;
              }
              else
                ( *rho )[ ( index + k ) ] += 0.0;
            }
        }
    }
#if defined(_OPENMP)
}
#endif
    free(tetra_crossed);
    free(tetra_crossed_srt);

}
