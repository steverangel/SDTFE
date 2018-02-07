#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <tiffio.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include "io.h"
#include "bitmap.h"
#include "constants.h"

// Creates a tiff image in RGB using the Jet colormap. 
// Uses log space values of the density
void write_image(char *file_name, double *rho, int grid_dim) {

    int i;
    for (i = 0; i < grid_dim * grid_dim; ++i) {
      if ( rho[i] < 0.0 ) {
        printf("Negative density!\n");
        return;
      }
      if (rho[i] != rho[i]) {
        printf("Nan denaity!\n");
        return;
      }
    }

    double max_val = 0.0;
    double min_val = DBL_MAX;
    for (i=0; i<grid_dim*grid_dim; ++i) {
        max_val = (rho[i] > max_val) ? rho[i] : max_val;
        min_val = (rho[i] > 0.0 && rho[i] < min_val) ? rho[i] : min_val;
    }

    int sampleperpixel = 1;
    unsigned char *image = (unsigned char*)malloc(grid_dim*grid_dim*sampleperpixel*sizeof(unsigned char));
    double max_val_log = log(max_val+1.0);
    double min_val_log = log(min_val+1.0);
    double nf = max_val_log-min_val_log;
    for (i=0; i<grid_dim*grid_dim; ++i) {
      if (rho[i]>0.0) {
        image[i] = (unsigned char)(128.0*(log(rho[i]+1.0)-min_val_log)/nf);
      }
      else {
        image[i]=0;
      }
    }

    TIFF *image_out= TIFFOpen( file_name, "w" );
    TIFFSetField( image_out, TIFFTAG_IMAGEWIDTH, grid_dim ); 
    TIFFSetField( image_out, TIFFTAG_IMAGELENGTH, grid_dim );
    TIFFSetField( image_out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel ); 
    TIFFSetField( image_out, TIFFTAG_BITSPERSAMPLE, 8 );
    TIFFSetField( image_out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );
    TIFFSetField( image_out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );
    TIFFSetField( image_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
    tsize_t linebytes = sampleperpixel*grid_dim;

    // The Jet colormap
    unsigned char r_tbl[256]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,131,135,139,143,147,151,155,159,163,167,171,175,179,183,187,191,195,199,203,207,211,215,219,223,227,231,235,239,243,247,251,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,251,247,243,239,235,231,227,223,219,215,211,207,203,199,195,191,187,183,179,175,171,167,163,159,155,151,147,143,139,135,131,128};

    unsigned char g_tbl[256]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,131,135,139,143,147,151,155,159,163,167,171,175,179,183,187,191,195,199,203,207,211,215,219,223,227,231,235,239,243,247,251,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,251,247,243,239,235,231,227,223,219,215,211,207,203,199,195,191,187,183,179,175,171,167,163,159,155,151,147,143,139,135,131,128,124,120,116,112,108,104,100,96,92,88,84,80,76,72,68,64,60,56,52,48,44,40,36,32,28,24,20,16,12,8,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    unsigned char b_tbl[256]={131,135,139,143,147,151,155,159,163,167,171,175,179,183,187,191,195,199,203,207,211,215,219,223,227,231,235,239,243,247,251,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,251,247,243,239,235,231,227,223,219,215,211,207,203,199,195,191,187,183,179,175,171,167,163,159,155,151,147,143,139,135,131,128,124,120,116,112,108,104,100,96,92,88,84,80,76,72,68,64,60,56,52,48,44,40,36,32,28,24,20,16,12,8,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    TIFFSetField( image_out, TIFFTAG_COLORMAP, r_tbl, g_tbl, b_tbl);

    unsigned char *buf = NULL; 
    if ( TIFFScanlineSize( image_out ) > linebytes )
        buf = ( unsigned char* )_TIFFmalloc( linebytes );
    else
        buf = ( unsigned char* )_TIFFmalloc( TIFFScanlineSize( image_out ) );

    TIFFSetField( image_out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize( image_out, grid_dim*sampleperpixel ) );

    uint32_t row;
    for ( row = 0; row < grid_dim; row++)
    {
        memcpy( buf, &image[ row * linebytes ], linebytes);   
        if ( TIFFWriteScanline( image_out, buf, row, 0 ) < 0 )
            break;
    }
    ( void ) TIFFClose( image_out );
    if ( buf )
        _TIFFfree( buf );

    free( image );
}


// Writes the raw data of the field
int write_rho( char* fname, double **rho, int rho_size ) {
    FILE *outfile;
    outfile = fopen( fname, "wb" );
    int bytes_written = fwrite( *rho, sizeof( double ), rho_size, outfile );
    fclose( outfile );

    return bytes_written;
}

