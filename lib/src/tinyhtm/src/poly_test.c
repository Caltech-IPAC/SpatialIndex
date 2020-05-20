
//  $Id: poly_test.c,v 20130613.99094427 2013/06/13 16:44:27 modell ipac $

/** \file
    \brief      Check Polygon (e.g. for convexity)

    \authors    Mark O'Dell
    \copyright  IPAC/Caltech
  */

//  expects RA-Dec pairs for vertices
//  followed by a blank line when done
//  returns transformed coordinates (sc>v3) 
//  and conclusion(s) about polygonm -- 
//  so far only re convexity

//  $ poly_test        # e.g.
//  10 10
//  20 10
//  20 20
//  10 20
//  
//  
//   # = 04
//  
//    00     10.000000000     10.000000000
//    01     20.000000000     10.000000000
//    02     20.000000000     20.000000000
//    03     10.000000000     20.000000000
//  
//    00      0.969846310      0.171010072      0.173648178
//    01      0.925416578      0.336824089      0.173648178
//    02      0.883022222      0.321393805      0.342020143
//    03      0.925416578      0.163175911      0.342020143
//  
//   Polygon IS  Convex

#if 1    //  for syntax folding

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <readline/readline.h>
#include <readline/history.h>

#include "tinyhtm/geometry.h"
#include "tinyhtm/htm.h"

#ifndef STRINGIFY
#define STRING(whatever)         # whatever
#define STRINGIFY(whatever) STRING(whatever)
#endif

#ifndef FALSE
#define FALSE 0
#define TRUE  (!FALSE)
#endif

#define MAX_LINE  1024
#define NP        16
#define ALWAYS    1

#define ERROUT(nm)     printf( "\n* ERROR -- malloc(" STRINGIFY(nm) ")\n\n" )
#define ERRMSG(ty,nm)  if ((struct ty *)0==nm) ERROUT(nm)

#define A_LT(ty,nm)         struct ty * nm 
#define A_RT(nu,ty,nm)     (struct ty *)malloc( nu * sizeof(struct ty))

#define NPALLOC(a_n,a_type,a_name) A_LT(a_type,a_name) = A_RT(a_n,a_type,a_name); ERRMSG(a_type,a_name)
#define  PALLOC(a_type,a_name)     NPALLOC(1, a_type, a_name )

#ifdef __cplusplus
extern "C" {
#endif

#endif

enum htm_errcode ec = HTM_OK;

char * mReadline( char * in ) {
    char * txt;
    char * out = (char *)malloc( MAX_LINE );
    
    //  prompt
    printf( "%s", in );
    
    if ((txt = fgets( out, MAX_LINE, stdin))) {
         while(*txt == ' ' || *txt == '\t') {
            //  trim preceding whitespace 
             txt++;
         }
         return txt;    
     }
     else {
         return 0;
     }
}

int is_convex( struct htm_v3 * verts, int n_vert ) {
    return( htm_v3_convex( verts, n_vert, &ec ));
}

void report_pts( struct htm_sc * p, int np ) {
    int ii;
    for ( ii=0; ii<np; ii++ ) {
        printf( "\n  %02d  %15.9lf  %15.9lf", ii, p[ii].lon, p[ii].lat );
        printf( "\n" );
    }
    printf( "\n" );
}

void report_verts( struct htm_v3 * v, int nv ) {
    int ii;
    for ( ii=0; ii<nv; ii++ ) {
        printf( "\n  %02d  %15.9lf  %15.9lf  %15.9lf", ii, v[ii].x, v[ii].y, v[ii].z );
    }
    printf( "\n" );
}

int main( int argc, char *argv[], char  *envp[] ) {

    char * rcvd;

    int    nr;
    int    nn  = 0;

    double ra  = 0.0;
    double dec = 0.0;

    char   given[ MAX_LINE ];

    int convexity = FALSE;

    struct htm_s2cpoly * poly;

    struct htm_sc         pt;
    struct htm_v3         vert;

    NPALLOC( NP, htm_sc, pts   );

    NPALLOC( NP, htm_v3, verts );

    while ( ALWAYS ) {

//      rcvd = readline( "" );
        rcvd = mReadline( "" );

        if (0 == rcvd) break;

        nr = sscanf( rcvd, "%lf %lf", &ra, &dec );
        free( rcvd );
        if (2 != nr) break;

        ec = htm_sc_init( &pt, ra, dec );
        if (HTM_OK != ec) {
            printf( 
                "Did not find acceptable spherical coordinates (%lf,%lf) -- %s\n\n",
                ra, dec,
                htm_errmsg(ec) );
            exit( -1 );
        }

        pts[   nn ] = pt;

        ec = htm_sc_tov3( &vert, &pt );
        if (HTM_OK != ec) {
            printf( 
                "Failed to convert spherical coordinates (%lf,%lf) to unit vector (%lf,%lf,%lf) -- %s\n\n",
                pt.lon,  pt.lat,
                vert.x,  vert.y, vert.z, 
                htm_errmsg(ec) );
            exit( -2 );
        }

        verts[ nn ] = vert;

        nn++;

        if (NP <= nn) break;
    }

    printf( "\n" );

    printf( "\n # = %02d", nn );

    printf( "\n" );

    fflush( stdout );

    report_pts( pts, nn );

    report_verts( verts, nn );

    poly = htm_s2cpoly_hull( verts, (size_t)nn, &ec );
//  free(verts);
    if (ec != HTM_OK) {
        printf( "\nCould not compute convex hull: %s", htm_errmsg(ec) );
        exit( -3 );
    }

    report_verts( poly->ve, nn );

//      printf( "\n Polygon (%d) is ", nn );
//      if ( ! is_convex( poly, nn )) {
//          printf( "NOT " );
//      }
//      printf( "Convex\n" );

    //  printf( "\n" );

    printf( "\n Polygon.verts (%d) is ", nn );
    if ( ! is_convex( poly->ve, nn )) {
        printf( "NOT " );
    }
    printf( "Convex\n" );

    convexity = is_convex( verts, nn );

    printf( "\n Polygon " );
    //  
    if ( convexity ) {
        printf( "IS  " );
    }
    else {
        printf( "is NOT " );
    }
    printf( "Convex\n\n" );

    fflush( stdout );
    return( 0 );
}

#if 1    //  for syntax folding

#ifdef __cplusplus
}
#endif

#endif

//  vi: set tabstop=4 shiftwidth=4 expandtab :

