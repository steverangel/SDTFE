#include <assert.h>
#include "qhull_a.h"
#include "io.h"
#include "bitmap.h"
#include "constants.h"

// This function triangulates the point list ( pts ) and constructs a tetrahedron list ( tetra_lst )
int *triangulate( double *particle_data, int n_particles, size_t *num_facets, char *options ) {

    /* ******************************************************************************
     * The Qhull delanuay triangulation uses (+1) dimensions
     * one plus the dim is for projecting points onto a
     * higher dimensional (one dimenstion higher than the input data) parabaloid
     * */

    int i, j, k;
    (*num_facets) = 0;
    int verts1[4];
    int verts2[4];
    int comp[4];
    int curlong, totlong;

    // called before error handling initialized
    qh_init_A(stdin, stdout, stderr, 0, NULL);

    // set flags and constants from command line
    qh_initflags(options);

    // called after points are defined, dim is 4
    qh_init_B (particle_data, n_particles, 4, 0);

    // project points to paraboloid for Delaunay triangulation
    qh_setdelaunay(4, n_particles, particle_data);

    // construct the convex hull of a set of points
    qh_qhull();
    qh_triangulate(); 
    
    // for qhull macros to iterate over qhull data structures
    facetT  *facet, *neighbor;
    vertexT *vertex, **vertexp;
    ridgeT  *ridge, **ridgep;

    // count the number of facets in the delanuay triangulation
    //  voronoi tetrahedra are facets that are not upper delaunay
    FORALLfacets {
        if ( !facet->upperdelaunay ) {
            ++(*num_facets);
        }
    }

   // ( A, B, C, D, neighbor_ba, neighbor_ac, neighbor_dc, neighbor_bd ) 
    int *tetra_data = (int*)malloc((*num_facets) * sizeof(int) * 8);

    i = 0;
    FORALLfacets {
        assert ( facet->simplicial );
        if (!facet->upperdelaunay) {
            facet->visitid = i++;
        }
    }
    // copy the tetrahedra points and re-orient them if needed
    FORALLfacets {
        if ( !facet->upperdelaunay ) {
            j = 0;
            FOREACHvertex_( facet->vertices ) {
                verts1[ j++ ] = qh_pointid( vertex->point );
            }
            // sanity test
            assert( j == 4 );
            // find the index of the tetrahedron
            double ab[3], ac[3], ad[3], norm[3];
            //#pragma simd
            for ( k = 0; k < 3; ++k ) {
                ab[ k ] = particle_data[ verts1[ B ] * 4 + k ] - particle_data[ verts1[ A ] * 4 + k ];
                ac[ k ] = particle_data[ verts1[ C ] * 4 + k ] - particle_data[ verts1[ A ] * 4 + k ];
                ad[ k ] = particle_data[ verts1[ D ] * 4 + k ] - particle_data[ verts1[ A ] * 4 + k ];
            }

            // cross product
            norm[ 0 ] = ab[ 1 ]*ac[ 2 ] - ab[ 2 ]*ac[ 1 ]; //uy*vz - uz*vy
            norm[ 1 ] = ab[ 2 ]*ac[ 0 ] - ab[ 0 ]*ac[ 2 ]; //uz*vx - ux*vz
            norm[ 2 ] = ab[ 0 ]*ac[ 1 ] - ab[ 1 ]*ac[ 0 ]; //ux*vy - uy*vx

            // dot product
            double temp = 0;
            //#pragma simd
            for ( k = 0; k < 3; ++k ) {
                temp += norm[ k ] * ad[ k ];
            }
            
            if ( temp > 0 ) {    // TRI ABC is clockwise
                //#pragma simd
                for ( k = 0; k < 3; ++k )
                    norm[ k ] = -norm[ k ];
                
                int temp2 = verts1[ B ];
                verts1[ B ] = verts1[ C ];
                verts1[ C ] = temp2;
            }

            i = facet->visitid;

            //#pragma simd
            for ( k = 0; k < 4; ++k )
              tetra_data[ i * 8 + k ] = verts1[ k ];
        }
    }

    // compute the neighbor information
    FORALLfacets {
        if ( !facet->upperdelaunay ) {

            // find the index of the tetrahedron
            i = facet->visitid;

            tetra_data[ i * 8 + 4 + BDC ] = -1;
            tetra_data[ i * 8 + 4 + ACD ] = -1;
            tetra_data[ i * 8 + 4 + ADB ] = -1;
            tetra_data[ i * 8 + 4 + ABC ] = -1;

            verts1[ 0 ] = tetra_data[ i * 8 + A ];
            verts1[ 1 ] = tetra_data[ i * 8 + B ];
            verts1[ 2 ] = tetra_data[ i * 8 + C ];
            verts1[ 3 ] = tetra_data[ i * 8 + D ];

            qh_makeridges( facet );
            FOREACHridge_( facet->ridges ) {
                neighbor = otherfacet_( ridge, facet );
                if ( !neighbor->upperdelaunay ) {
                    j = neighbor->visitid;

                    comp[ 0 ] = 0;
                    comp[ 1 ] = 0;
                    comp[ 2 ] = 0;
                    comp[ 3 ] = 0;
                   
                    verts2[ 0 ] = tetra_data[ j * 8 + A ];
                    verts2[ 1 ] = tetra_data[ j * 8 + B ];
                    verts2[ 2 ] = tetra_data[ j * 8 + C ];
                    verts2[ 3 ] = tetra_data[ j * 8 + D ];

                    for ( k = 0; k < 4; k++ ) {
                        if ( ( verts1[ k ] == verts2[ 0 ] ) || ( verts1[ k ] == verts2[ 1 ] ) || ( verts1[ k ] == verts2[ 2 ] ) || ( verts1[ k ] == verts2[ 3 ] ) ) {
                            comp[ k ] = 1;
                        }
                    }

                    if ( comp[ 0 ] + comp[ 1 ] + comp[ 2 ] + comp[ 3 ] != 3) {
                        qh_freeqhull(!qh_ALL);
                        qh_memfreeshort( &curlong, &totlong );
                        printf("TRIANGULATION ERROR: comp_sum = %d\n", comp[ 0 ] + comp[ 1 ] + comp[ 2 ] + comp[ 3 ] );
                        return NULL;
                    }

                    if ( comp[ 0 ] == 0 ) {
                        assert( tetra_data[ i * 8 + 4 + BDC ] == -1 );
                        tetra_data[ i * 8 + 4 + BDC ] = j;
                    }
                    else if ( comp[ 1 ] == 0 ) {
                        assert( tetra_data[ i * 8 + 4 + ACD ] == -1 );
                        tetra_data[ i * 8 + 4 + ACD ] = j;
                    }
                    else if ( comp[ 2 ] == 0 ) {
                        assert( tetra_data[ i * 8 + 4 + ADB ] == -1 );
                        tetra_data[ i * 8 + 4 + ADB ] = j;
                    }
                    else if (comp[ 3 ] == 0) {
                        assert( tetra_data[ i * 8 + 4 + ABC ] == -1 );
                        tetra_data[ i * 8 + 4 + ABC ] = j;
                    } 
                }
            }
        }
    }

    qh_freeqhull(!qh_ALL);
    qh_memfreeshort( &curlong, &totlong );

    return tetra_data;
}

void pre_compute_vol( double **particle_data, int n_particles, int *tetra_data, int n_tetra, float p_mass ) {
    double vol;
    int i;

    for ( i = 0; i < n_particles; ++i )
       (*particle_data)[ i * 4 + 3 ] = 0;

    for ( i = 0; i < n_tetra; ++i ) {
        vol = ( 
                ( (*particle_data)[  tetra_data[ i * 8 + A ] * 4 + 0 ]  - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 0 ]  ) *
                (
                    ( (*particle_data)[  tetra_data[ i * 8 + B ] * 4 + 1 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 1 ] ) *
                    ( (*particle_data)[  tetra_data[ i * 8 + C ] * 4 + 2 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 2 ] ) - 
                    ( (*particle_data)[  tetra_data[ i * 8 + B ] * 4 + 2 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 2 ] ) *
                    ( (*particle_data)[  tetra_data[ i * 8 + C ] * 4 + 1 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 1 ] ) 
                ) +
                ( (*particle_data)[  tetra_data[ i * 8 + A ] * 4 + 1 ]  - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 1 ]  ) *
                (
                    ( (*particle_data)[  tetra_data[ i * 8 + B ] * 4 + 2 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 2 ] ) *
                    ( (*particle_data)[  tetra_data[ i * 8 + C ] * 4 + 0 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 0 ] ) - 
                    ( (*particle_data)[  tetra_data[ i * 8 + B ] * 4 + 0 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 0 ] ) *
                    ( (*particle_data)[  tetra_data[ i * 8 + C ] * 4 + 2 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 2 ] )
                ) +
                ( (*particle_data)[  tetra_data[ i * 8 + A ] * 4 + 2 ]  - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 2 ]  ) *
                (
                    ( (*particle_data)[  tetra_data[ i * 8 + B ] * 4 + 0 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 0 ] ) *
                    ( (*particle_data)[  tetra_data[ i * 8 + C ] * 4 + 1 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 1 ] ) - 
                    ( (*particle_data)[  tetra_data[ i * 8 + B ] * 4 + 1 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 1 ] ) *
                    ( (*particle_data)[  tetra_data[ i * 8 + C ] * 4 + 0 ] - (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 0 ] )
                )
             ) / 6;

        vol = fabs(vol);

        (*particle_data)[  tetra_data[ i * 8 + A ] * 4 + 3 ] += vol;
        (*particle_data)[  tetra_data[ i * 8 + B ] * 4 + 3 ] += vol;
        (*particle_data)[  tetra_data[ i * 8 + C ] * 4 + 3 ] += vol;
        (*particle_data)[  tetra_data[ i * 8 + D ] * 4 + 3 ] += vol;
    }
 
    for ( i = 0; i < n_particles; ++i )
       (*particle_data)[ i * 4 + 3 ] = p_mass * ( 1 + DIM ) / (*particle_data)[ i * 4 + 3 ]; 
}


static inline double determ3x3(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
   return a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h; 
}


double *gradients( int *tetra_data, int n_tetra, double *particle_data ) {
  double *nabla = (double*)malloc(n_tetra*sizeof(double)*3);
  double D_xb, D_yb, D_zb, D_xc, D_yc, D_zc, D_xd, D_yd, D_zd;
  double D_rb, D_rc, D_rd;
  int i;
  for ( i = 0; i < n_tetra; ++i ) {
    D_xb = particle_data[ tetra_data[ i * 8 + B ] * 4 + 0 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 0 ];
    D_yb = particle_data[ tetra_data[ i * 8 + B ] * 4 + 1 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 1 ];
    D_zb = particle_data[ tetra_data[ i * 8 + B ] * 4 + 2 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 2 ];
    D_xc = particle_data[ tetra_data[ i * 8 + C ] * 4 + 0 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 0 ];
    D_yc = particle_data[ tetra_data[ i * 8 + C ] * 4 + 1 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 1 ];
    D_zc = particle_data[ tetra_data[ i * 8 + C ] * 4 + 2 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 2 ];
    D_xd = particle_data[ tetra_data[ i * 8 + D ] * 4 + 0 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 0 ];
    D_yd = particle_data[ tetra_data[ i * 8 + D ] * 4 + 1 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 1 ];
    D_zd = particle_data[ tetra_data[ i * 8 + D ] * 4 + 2 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 2 ];

    D_rb = particle_data[ tetra_data[ i * 8 + B ] * 4 + 3 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 3 ];
    D_rc = particle_data[ tetra_data[ i * 8 + C ] * 4 + 3 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 3 ];
    D_rd = particle_data[ tetra_data[ i * 8 + D ] * 4 + 3 ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + 3 ];

    double denom = determ3x3(D_xb, D_yb, D_zb, D_xc, D_yc, D_zc, D_xd, D_yd, D_zd);

    nabla[ i * 3 + 0 ] = determ3x3(D_rb, D_yb, D_zb, D_rc, D_yc, D_zc, D_rd, D_yd, D_zd) / denom;
    nabla[ i * 3 + 1 ] = determ3x3(D_xb, D_rb, D_zb, D_xc, D_rc, D_zc, D_xd, D_rd, D_zd) / denom;
    nabla[ i * 3 + 2 ] = determ3x3(D_xb, D_yb, D_rb, D_xc, D_yc, D_rc, D_xd, D_yd, D_rd) / denom;
  }
  return nabla;
}

// This function determines which tetrahedra have a face on the forward hull
int find_forward_hull( int *tetra_data, int n_tetra, double *particle_data, int n_particles, word_t *bitmap ) {
  int num_forward_faces = 0;
  int i, j, k;
  int written = 0;
  int neighbor_id;

  double ab[ 3 ], ac [ 3 ], norm[ 3 ], t1[ 3 ], t2[ 3 ], unit_z[ 3 ];

  // count the number of forward faces of the hull
  // find the faces of the hull that are in the -z direction
  // i.e. whose surface normal dotted with [0 0 1] is negative 
  for ( i = 0; i < n_tetra; ++i ) {
    // tetra row ( A, B, C, D, N1, N2, N3, N4 )
    for ( j = 0; j < 4; ++j ) {
      neighbor_id = tetra_data[ i * 8 + ( 4 + j ) ];
      if ( neighbor_id == -1 ) {
        switch ( j ) {
          case ABC:
            for ( k = 0; k < 3; ++k ) {
              t1[ k ] = particle_data[ tetra_data[ i * 8 + B ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + k ];
              t2[ k ] = particle_data[ tetra_data[ i * 8 + C ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + k ];
            }
            break;
          case ADB:
            for ( k = 0; k < 3; ++k ) {
              t1[ k ] = particle_data[ tetra_data[ i * 8 + D ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + k ];
              t2[ k ] = particle_data[ tetra_data[ i * 8 + B ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + k ];
            }
            break;
          case ACD:
            for ( k = 0; k < 3; ++k ) {
              t1[ k ] = particle_data[ tetra_data[ i * 8 + C ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + k ];
              t2[ k ] = particle_data[ tetra_data[ i * 8 + D ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + A ] * 4 + k ];
            }
            break;
          case BDC:
            for ( k = 0; k < 3; ++k ) {
              t1[ k ] = particle_data[ tetra_data[ i * 8 + D ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + B ] * 4 + k ];
              t2[ k ] = particle_data[ tetra_data[ i * 8 + C ] * 4 + k ] - particle_data[ tetra_data[ i * 8 + B ] * 4 + k ];
            }
            break;
          default:
              fprintf( stderr, "Something really crazy happened!\n" );
              exit( EXIT_FAILURE ); 
        }
        if ( t1[ 0 ]*t2[ 1 ] - t1[ 1 ]*t2[ 0 ] < 0.0) {
            set_bit( bitmap, i * 4 + j );
            ++num_forward_faces;
        }
      }
    }
  }
  return num_forward_faces;
}

int *copy_fh_triangles( int *tetra_data, int n_tetra, double *particle_data, int n_particles, int n_tri, word_t *bitmap ) {
  int *tri_data = (int*)malloc(n_tri*8*sizeof(int));

  int T_id = 0;
  int a, b, c;
  // copy the 3-d data triangulation into the 2-d forward hull projection
  int i, j;
  for ( i = 0; i < n_tetra; ++i ) {
    for ( j = 0; j < 4; ++j ) {
      if ( get_bit( bitmap, i * 4 + j ) ) {    
        switch ( j ) {
          case ABC:
            a = tetra_data[ i * 8 + A ]; b = tetra_data[ i * 8 + B ]; c = tetra_data[ i * 8 + C ];
            break;
          case ADB:
            a = tetra_data[ i * 8 + A ]; b = tetra_data[ i * 8 + D ]; c = tetra_data[ i * 8 + B ];
            break;
          case ACD:
            a = tetra_data[ i * 8 + A ]; b = tetra_data[ i * 8 + C ]; c = tetra_data[ i * 8 + D ];
            break;
          case BDC:
            a = tetra_data[ i * 8 + B ]; b = tetra_data[ i * 8 + D ]; c = tetra_data[ i * 8 + C ];
            break;
          default:
            break;
        }
        tri_data[ T_id * 8 + A ] = a;
        tri_data[ T_id * 8 + B ] = b;
        tri_data[ T_id * 8 + C ] = c;
        // set the neighbors to -1
        tri_data[ T_id * 8 + 3 ] = -1;
        tri_data[ T_id * 8 + 4 ] = -1;
        tri_data[ T_id * 8 + 5 ] = -1;
        tri_data[ T_id * 8 + 6 ] = i;
        T_id++;
      }
    }
  }
  //*******************************************************************************
  // find the neighbors in the set of faces (exhaustive, but ralatively small in number)
  for ( i = 0; i < T_id; ++i ) {
    for ( j = 0; j < T_id; ++j ) {
      if ( i != j ) {
        int comp[ 3 ];
        comp[ 0 ] = 0;
        comp[ 1 ] = 0; 
        comp[ 2 ] = 0;

        if ( ( tri_data[ i * 8 + A ] == tri_data[ j * 8 + A ] ) ||
           (   tri_data[ i * 8 + A ] == tri_data[ j * 8 + B ] ) ||
           (   tri_data[ i * 8 + A ] == tri_data[ j * 8 + C ] ) )
        comp[ 0 ] = 1;
        if ( ( tri_data[ i * 8 + B ] == tri_data[ j * 8 + A ] ) ||
           (   tri_data[ i * 8 + B ] == tri_data[ j * 8 + B ] ) ||
           (   tri_data[ i * 8 + B ] == tri_data[ j * 8 + C ] ) )
        comp[ 1 ] = 1;
        if ( ( tri_data[ i * 8 + C ] == tri_data[ j * 8 + A ] ) ||
           (   tri_data[ i * 8 + C ] == tri_data[ j * 8 + B ] ) ||
           (   tri_data[ i * 8 + C ] == tri_data[ j * 8 + C ] ) )
        comp[ 2 ] = 1;

        if ( ( comp[ 0 ] == 1 ) && ( comp[ 1 ] == 1 ) )
          tri_data[ i * 8 + ( 3 + BA ) ] = j;

        if ( ( comp[ 0 ] == 1 ) && ( comp[ 2 ] == 1 ) )
          tri_data[ i * 8 + ( 3 + AC ) ] = j;

        if ( ( comp[ 1 ] == 1 ) && ( comp[ 2 ] == 1 ) )
          tri_data[ i * 8 + ( 3 + CB ) ] = j;

        // sanity test    
        assert( ( comp[ 0 ] + comp[ 1 ] + comp[ 2 ] ) < 3 );
      }
    }
  }
  // end
  return tri_data;
}
