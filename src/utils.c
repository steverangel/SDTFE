#include <assert.h>
#include <stdio.h>
#include <float.h>
#include "io.h"
#include "sample.h"
#include "bitmap.h"
#include <math.h>

#define SWAP( a, b, type ) { type tmp; tmp = ( a ); ( a ) = ( b ); ( b ) = tmp; }
#define INNER_BOX_OVERLOAD 0.2f
#define FLOAT_MAX_ULP 2
#define SIMD_WIDTH 4

float *base_arr;

typedef enum { false, true } bool;    

typedef union Flt_t{
   float f;
   int i;
} Flt_t;

bool f_almost_eq( float A, float B ) {
    double abs_diff = fabs( A - B );
    if ( abs_diff <= 2*FLT_EPSILON )
        return true;

    Flt_t uA, uB;
    uA.f = A;
    uB.f = B;

    // check the signs
    if ( ( uA.i >> 31 ) != ( uB.i >> 31 ) )
        return false;

    int ulpsDiff = abs( uA.i - uB.i );
    if ( ulpsDiff <= FLOAT_MAX_ULP ) 
        return true;

    return false;
}

int coord_comptor( const void * a, const void * b ) {

    int *ptr_a = *( int** )a;
    int *ptr_b = *( int** )b;

    float fa = base_arr[ *ptr_a ];
    float fb = base_arr[ *ptr_b ];

    if ( f_almost_eq( fa, fb ) ) {
        return (int)( ptr_a - ptr_b );
    }

    return ( fa > fb ) - ( fa < fb );
}

// stable sort of the particles in parallel arrays {x, y, z} structure-of-arrays
void sort_rows( particle_data *particles, int num_particles, int **idx ) {
    
    int i;

    int **idx_ptr = (int**)malloc( num_particles * sizeof( int* ) );
    int *idx_tmp = (int*)malloc( num_particles * sizeof( int ) );

    // Z dim
    for ( i = 0; i < num_particles; ++i )
        idx_ptr[ i ] = &(*idx)[ i ];
    base_arr = particles->z;
    qsort( idx_ptr, num_particles, sizeof( int* ), coord_comptor );
    for ( i = 0; i < num_particles; ++i )
        idx_tmp[ i ] = *(idx_ptr[ i ]);
    for ( i = 0; i < num_particles; ++i )
        (*idx)[ i ] = idx_tmp[ i ];

    // Y dim
    for ( i = 0; i < num_particles; ++i )
        idx_ptr[ i ] = &(*idx)[ i ];
    base_arr = particles->y;
    qsort( idx_ptr, num_particles, sizeof( int* ), coord_comptor );
    for ( i = 0; i < num_particles; ++i )
        idx_tmp[ i ] = *(idx_ptr[ i ]);
    for ( i = 0; i < num_particles; ++i )
        (*idx)[ i ] = idx_tmp[ i ];

    // X dim
    for ( i = 0; i < num_particles; ++i )
        idx_ptr[ i ] = &(*idx)[ i ];
    base_arr = particles->x;
    qsort( idx_ptr, num_particles, sizeof( int* ), coord_comptor );
    for ( i = 0; i < num_particles; ++i )
        idx_tmp[ i ] = *(idx_ptr[ i ]);
    for ( i = 0; i < num_particles; ++i )
        (*idx)[ i ] = idx_tmp[ i ];

    free( idx_ptr );
    free( idx_tmp );
}

int unique_rows( particle_data *particles, int num_particles ) {
    int num_unique;
    int i;

    int *coord_idx = (int*)malloc( num_particles * sizeof( int ) );

    for ( i = 0; i < num_particles; ++i ) {
        coord_idx[ i ] = i;
    }

    sort_rows( particles, num_particles, &coord_idx );

    bool eq[ 3 ];

    word_t *delete_map = create_bitmap( num_particles );
    for ( i = 1; i < num_particles; ++i ) {

        eq[ 0 ] = f_almost_eq( particles->x[ coord_idx[ i ] ], particles->x[ coord_idx[ i - 1 ] ] );
        eq[ 1 ] = f_almost_eq( particles->y[ coord_idx[ i ] ], particles->y[ coord_idx[ i - 1 ] ] );
        eq[ 2 ] = f_almost_eq( particles->z[ coord_idx[ i ] ], particles->z[ coord_idx[ i - 1 ] ] );

        if ( eq[ 0 ] && eq[ 1 ] && eq[ 2 ] ) {
            set_bit( delete_map, i - 1 );
        }
    }
    
    num_unique = num_particles;
    for ( i = 0; i < num_particles; ++i ) {
        if ( get_bit( delete_map, i ) ) {
            /*while ( get_bit( delete_map, num_unique - 1 ) )
                --num_unique;
            particles->x[ i ] = particles->x[ num_unique - 1 ];
            particles->y[ i ] = particles->y[ num_unique - 1 ];
            particles->z[ i ] = particles->z[ num_unique - 1 ];
            --num_unique;*/
            particles->x[ i ] += rand() % 2 ? -1e-6 : 1e-6;
            particles->y[ i ] += rand() % 2 ? -1e-6 : 1e-6;
            particles->z[ i ] += rand() % 2 ? -1e-6 : 1e-6;
        }
    }

    free( coord_idx );
    free( delete_map );

    return num_unique;
}

int comptor( const void * a,const void * b ) {
    return ( *( int* )a - *( int* )b );
}

void find_partition_boundaries( int num_points, int depth, int **partition_boundaries ) {

    int i, j, mid, num_partitions;

    (*partition_boundaries)[ 0 ] = 0;
    (*partition_boundaries)[ 1 ] = num_points - 1;
    
    for ( i = 0; i < depth; ++i ) {
        num_partitions = ( int )pow( 2.0, i );

        for ( j = 0; j < num_partitions; ++j ) {
            mid = (*partition_boundaries)[ j * 2 ] + 
                ( ( (*partition_boundaries)[ j * 2 + 1 ] - (*partition_boundaries)[ j * 2 ] ) >> 1 );

            (*partition_boundaries)[ ( 2 * num_partitions ) + ( 2 * j ) + 0 ] = mid;
            (*partition_boundaries)[ ( 2 * num_partitions ) + ( 2 * j ) + 1 ] = mid + 1;
        }
        qsort( *partition_boundaries, 2 * ( 2 * num_partitions ), sizeof( int ), comptor );
    }
}

int partition( particle_data *particles, int **tree_index, int order_dim, int left, int right, int pivot_index ) {
    
    int i;
    float pivot_value;
    float *data;

    if ( order_dim == 0 )
        data = particles->x;
    else if ( order_dim == 1 )
        data = particles->y;
    else
        data = particles->z;

    pivot_value = data[ (*tree_index)[ pivot_index ] ];
    int store_index = left;

    SWAP( (*tree_index)[ pivot_index ], (*tree_index)[ right ], int );

    for ( i = left; i < right; ++i ) {
        if ( data[ (*tree_index)[ i ] ] < pivot_value ) {
            SWAP( (*tree_index)[ store_index ], (*tree_index)[ i ], int );
            ++store_index;
        }
    }

    SWAP( (*tree_index)[ right ], (*tree_index)[ store_index ], int );

    return store_index;
}

void nth_element( particle_data *particles, int **tree_index, int order_dim, int left, int right, int n ) {

    int pivot_index;
    int rnd_pivot_index;
    if ( left == right )
        return;

    while ( 1 ) {
        assert( right >= left );
        rnd_pivot_index = left + ( rand() % ( right - left + 1 ) );
        pivot_index = partition( particles, tree_index, order_dim, left, right, rnd_pivot_index );

        if ( n == pivot_index )
            return;
        else if ( n < pivot_index )
            right = pivot_index - 1;
        else
            left = pivot_index + 1;
    }
}

int longest_dim( particle_data *particles, int *tree_index, int left, int right ) {

    float min[3], max[3], dist[3];
    int i, j;


    min[0] = particles->x[ tree_index[ left ] ];
    min[1] = particles->y[ tree_index[ left ] ];
    min[2] = particles->z[ tree_index[ left ] ];

    max[0] = min[0];
    max[1] = min[1];
    max[2] = min[2];

    float x, y, z;
    for ( j = left + 1; j <= right; ++j ) {
        x = particles->x[ tree_index[ j ] ];
        y = particles->y[ tree_index[ j ] ];
        z = particles->z[ tree_index[ j ] ];

        max[0] = ( max[0] < x ) ? x : max[0];
        max[1] = ( max[1] < y ) ? y : max[1];
        max[2] = ( max[2] < z ) ? z : max[2];

        min[0] = ( min[0] > x ) ? x : min[0];
        min[1] = ( min[1] > y ) ? y : min[1];
        min[2] = ( min[2] > z ) ? z : min[2];
    }

    dist[0] = max[0] - min[0];
    dist[1] = max[1] - min[1];
    dist[2] = max[2] - min[2];

    float max_dist = dist[0];
    int max_dim = 0;
    if ( max_dist < dist[1] ) {
        max_dist = dist[1];
        max_dim = 1;
    }

    if ( max_dist < dist[2] ) {
        max_dist = dist[2];
        max_dim = 2;
    }

    return max_dim;

}

void tree_order(  particle_data *particles, int **tree_index, int order_dim, int left, int right ) {

    int mid;
    int long_dim;
    int len = right - left;

    // base case
    if ( len == 1) 
        return;


    // divide
    mid = left + ( len >> 1 );

    long_dim = longest_dim( particles, *tree_index, left, right );
    nth_element( particles, tree_index, long_dim, left, right, mid );

    long_dim = longest_dim( particles, *tree_index, left, mid );
    tree_order( particles, tree_index, long_dim, left, mid );

    long_dim = longest_dim( particles, *tree_index, mid, right );
    tree_order( particles, tree_index, long_dim, mid, right );
    
}

void shuffle_partitions( int **tree_data, int **partition_boundaries, int num_partitions, int** size_partition_sample ) {

    int i, j, r, left, right;

    for ( i = 0; i < num_partitions / 2; ++i ) {
        left =  (*partition_boundaries)[ i * 2 ];
        right = (*partition_boundaries)[ i * 2 + 1];

        for ( j = left; j < left + (*size_partition_sample)[ i ]; ++j ) {
            r = ( rand() % ( right + 1 - j ) ) + j;
            SWAP( (*tree_data)[ j ], (*tree_data)[ r ], int );
        }
    }
}

double fRand( double fMin, double fMax ) {
    double f = ( double )rand() / RAND_MAX;
    return fMin + f * ( fMax - fMin );        
}

void sub_sample_data( particle_data *particles, int *tree_data_index, int *partition_boundaries, int num_partitions, 
    int *size_partition_sample, double **sample_data ) {

    int pid = 0;

    int i, j, r, index;
    for ( i = 0; i < num_partitions / 2; ++i ) {
        for ( j = 0; j < size_partition_sample[ i ]; ++j ) {

            index = tree_data_index[ partition_boundaries[ i * 2 ] + j ];

            //(*sample_data)[ pid * 4 + 0 ] = ( double )particles->x[ index ] + fRand( -1e-5, 1e-5 );
            //(*sample_data)[ pid * 4 + 1 ] = ( double )particles->y[ index ] + fRand( -1e-5, 1e-5 );
            //(*sample_data)[ pid * 4 + 2 ] = ( double )particles->z[ index ] + fRand( -1e-5, 1e-5 );
            (*sample_data)[ pid * 4 + 0 ] = ( double )particles->x[ index ] + ( rand() % 2 == 0 ? -( rand() % 1000 * 1e-8 ) : rand() % 1000 * 1e-8 );
            (*sample_data)[ pid * 4 + 1 ] = ( double )particles->y[ index ] + ( rand() % 2 == 0 ? -( rand() % 1000 * 1e-8 ) : rand() % 1000 * 1e-8 );
            (*sample_data)[ pid * 4 + 2 ] = ( double )particles->z[ index ] + ( rand() % 2 == 0 ? -( rand() % 1000 * 1e-8 ) : rand() % 1000 * 1e-8 );
            ++pid;
        }
    }
}

// Simply count all the particles inside a box centered at x,y,z of length box_len
int count_data_by_vol( particle_data *particles, int num_particles, float center_x, float center_y, float center_z, double box_len, double box_len_z ) {

    int i, j;

    float min_coords[ 3 ];
    float max_coords[ 3 ];

    // the half box length -- has an overload for density edges
    float center_offset_xy = box_len / 2 + box_len * INNER_BOX_OVERLOAD;
    float center_offset_z = box_len_z / 2 + box_len_z * INNER_BOX_OVERLOAD;

    // define the min and max coordinates of the box, i.e., front-left-bottom and back-right-top corners
    min_coords[ 0 ] = center_x - center_offset_xy;
    min_coords[ 1 ] = center_y - center_offset_xy;
    min_coords[ 2 ] = center_z - center_offset_z;
    max_coords[ 0 ] = center_x + center_offset_xy;
    max_coords[ 1 ] = center_y + center_offset_xy;
    max_coords[ 2 ] = center_z + center_offset_z;

    int num_particles_in_vol = 0;

    for ( i = 0; i < num_particles; ++i ) { 
      if (
          ( particles->x[ i ] < max_coords[ 0 ] )&&
          ( particles->y[ i ] < max_coords[ 1 ] )&&
          ( particles->z[ i ] < max_coords[ 2 ] )&&
          ( particles->x[ i ] > min_coords[ 0 ] )&&
          ( particles->y[ i ] > min_coords[ 1 ] )&&
          ( particles->z[ i ] > min_coords[ 2 ] )
      ){
        ++num_particles_in_vol;
      }
    }
    return num_particles_in_vol;
}
                                                                    

// Simply copy all the particles inside a box centered at x,y,z of length box_len
// copied into dtfe format -- flat array of double where each point is 4 consecutive elements, (x,y,z,extra) 
int dtfe_particles_in_vol( particle_data *particles, int num_particles, double **data, float center_x, float center_y, float center_z, double box_len, double box_len_z ) {

    int i, j;

    float min_coords[ 3 ];
    float max_coords[ 3 ];

    // the half box length -- has an overload for density edges
    float center_offset_xy = box_len / 2 + box_len * INNER_BOX_OVERLOAD;
    float center_offset_z = box_len_z / 2 + box_len_z * INNER_BOX_OVERLOAD;

    // define the min and max coordinates of the box, i.e., front-left-bottom and back-right-top corners
    min_coords[ 0 ] = center_x - center_offset_xy;
    min_coords[ 1 ] = center_y - center_offset_xy;
    min_coords[ 2 ] = center_z - center_offset_z;
    max_coords[ 0 ] = center_x + center_offset_xy;
    max_coords[ 1 ] = center_y + center_offset_xy;
    max_coords[ 2 ] = center_z + center_offset_z;

    int num_particles_in_vol = 0;

    for ( i = 0; i < num_particles; ++i ) { 
      if (
          ( particles->x[ i ] < max_coords[ 0 ] )&&
          ( particles->y[ i ] < max_coords[ 1 ] )&&
          ( particles->z[ i ] < max_coords[ 2 ] )&&
          ( particles->x[ i ] > min_coords[ 0 ] )&&
          ( particles->y[ i ] > min_coords[ 1 ] )&&
          ( particles->z[ i ] > min_coords[ 2 ] )
      ){
        (*data)[ num_particles_in_vol * 4 + 0 ] = particles->x[ i ];
        (*data)[ num_particles_in_vol * 4 + 1 ] = particles->y[ i ];
        (*data)[ num_particles_in_vol * 4 + 2 ] = particles->z[ i ];
        (*data)[ num_particles_in_vol * 4 + 3 ] = 0.0f;
        ++num_particles_in_vol;
      }
    }
    return num_particles_in_vol;
}
                                                                    
void shuffle_buff( double **data, int num_particles_in_buff, int num_shuffle ) {
    int i, j, r;
    double temp[ 4 ];
    for ( i = 1; i < num_shuffle; ++i ) {
        r = rand() % ( num_particles_in_buff - i ) + i;
        for ( j = 0; j < 4; ++j )
            temp[ j ] = (*data)[ r * 4 + j ];
        for ( j = 0; j < 4; ++j )
            (*data)[ r * 4 + j ] = (*data)[ ( i - 1 ) * 4 + j ];
        for ( j = 0; j < 4; ++j )
            (*data)[ ( i - 1 ) * 4 + j ] = temp[ j ];
    }
}

int sub_sample_data_by_vol( particle_data *particles, int *tree_data_index, int *partition_boundaries, int num_partitions, 
    int *size_partition_sample, int num_samples, double **sample_data, float center_x, float center_y, float center_z, double box_len ) {

    int i, j, r, index;

    float min_coords[ 3 ];
    float max_coords[ 3 ];

    float center_offset = box_len / 2 + box_len * 0.05;

    min_coords[ 0 ] = center_x - center_offset;
    min_coords[ 1 ] = center_y - center_offset;
    min_coords[ 2 ] = center_z - center_offset;
    max_coords[ 0 ] = center_x + center_offset;
    max_coords[ 1 ] = center_y + center_offset;
    max_coords[ 2 ] = center_z + center_offset;

    // this bitmap is a membership function for sample particles in the volume
    word_t *sample_map = create_bitmap( num_samples );

    int num_particles_in_vol = 0;
    // count the number of sample points in the sub-volume
    for ( i = 0; i < num_partitions / 2; ++i ) {
        for ( j = 0; j < size_partition_sample[ i ]; ++j ) {
            index = tree_data_index[ partition_boundaries[ i * 2 ] + j ];
            if ( particles->x[ index ] < max_coords[ 0 ] && particles->x[ index ] > min_coords[ 0 ] && \
                 particles->y[ index ] < max_coords[ 1 ] && particles->y[ index ] > min_coords[ 1 ] && \
                 particles->z[ index ] < max_coords[ 2 ] && particles->z[ index ] > min_coords[ 2 ] ) {
                 ++num_particles_in_vol;
                 set_bit( sample_map, index );
            }
        }
    }

    // this memory is free'd by the caller
    (*sample_data) = (double*)malloc( num_particles_in_vol * sizeof( double ) * 4 );

    int pid = 0;
    for ( i = 0; i < num_samples; ++i ) {
        if ( get_bit( sample_map, i ) ) {
            (*sample_data)[ pid * 4 + 0 ] = ( double )particles->x[ i ] + ( rand() % 2 == 0 ? -( rand() % 1000 * 1e-8 ) : rand() % 1000 * 1e-8 );
            (*sample_data)[ pid * 4 + 1 ] = ( double )particles->y[ i ] + ( rand() % 2 == 0 ? -( rand() % 1000 * 1e-8 ) : rand() % 1000 * 1e-8 );
            (*sample_data)[ pid * 4 + 2 ] = ( double )particles->z[ i ] + ( rand() % 2 == 0 ? -( rand() % 1000 * 1e-8 ) : rand() % 1000 * 1e-8 );
            ++pid;
        }
    }

    free( sample_map );

    return num_particles_in_vol;
}

// perform a 3-d rotation of the points
void rotate3d( double **sample_data, int num_points, double theta_x, double theta_y, double theta_z ) {

    int i, j, m;
    double x, y, z;
    double Qx[ 3 ][ 3 ];
    double Qy[ 3 ][ 3 ];
    double Qz[ 3 ][ 3 ];
    double Qt1[ 3 ][ 3 ];
    double Qt2[ 3 ][ 3 ];

    Qx[ 0 ][ 0 ] = 1.0f;            Qx[ 0 ][ 1 ] = 0.0f;           Qx[ 0 ][ 2 ] = 0.0f;
    Qx[ 1 ][ 0 ] = 0.0f;            Qx[ 1 ][ 1 ] = cos(theta_x);   Qx[ 1 ][ 2 ] = sin(theta_x);
    Qx[ 2 ][ 0 ] = 0.0f;            Qx[ 2 ][ 1 ] = -sin(theta_x);  Qx[ 2 ][ 2 ] = cos(theta_x);

    Qy[ 0 ][ 0 ] = cos(theta_y);    Qy[ 0 ][ 1 ] = 0.0f;            Qy[ 0 ][ 2 ] = -sin(theta_y); 
    Qy[ 1 ][ 0 ] = 0.0f;            Qy[ 1 ][ 1 ] = 1.0f;            Qy[ 1 ][ 2 ] = 0.0f; 
    Qy[ 2 ][ 0 ] = sin(theta_y);    Qy[ 2 ][ 1 ] = 0.0f;            Qy[ 2 ][ 2 ] = cos(theta_y);

    Qz[ 0 ][ 0 ] = cos(theta_z);    Qz[ 0 ][ 1 ] = sin(theta_z);    Qz[ 0 ][ 2 ] = 0.0f; 
    Qz[ 1 ][ 0 ] = -sin(theta_z);   Qz[ 1 ][ 1 ] = cos(theta_z);    Qz[ 1 ][ 2 ] = 0.0f; 
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

    for ( i = 0; i < num_points; ++i ) {
        x = (*sample_data)[ i * 4 + 0 ];
        y = (*sample_data)[ i * 4 + 1 ];
        z = (*sample_data)[ i * 4 + 2 ];
        (*sample_data)[ i * 4 + 0 ] = x * Qt2[ 0 ][ 0 ] + y * Qt2[ 0 ][ 1 ] + z * Qt2[ 0 ][ 2 ];
        (*sample_data)[ i * 4 + 1 ] = x * Qt2[ 1 ][ 0 ] + y * Qt2[ 1 ][ 1 ] + z * Qt2[ 1 ][ 2 ];
        (*sample_data)[ i * 4 + 2 ] = x * Qt2[ 2 ][ 0 ] + y * Qt2[ 2 ][ 1 ] + z * Qt2[ 2 ][ 2 ];
    }
}
