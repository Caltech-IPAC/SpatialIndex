#ifndef SPTQUERY_H
#define SPTQUERY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>

#include <tinyhtm.h>

#define HTM 0
#define HPX 1

#define DECIMAL 0
#define BASE4   1

#define CONE    0
#define POLYGON 1

#define MAXSTR 4096

#define DTR  0.0174532925199432957692369076849


 struct sptSkiplist
 {
    int     index;
    int64_t count;
 };

struct sptRanges
{
   int64_t min;
   int64_t max;

   int removeFlag;

   int next;
   int prev;
};


struct hpxRanges
{
   int64_t min;
   int64_t max;
};


typedef struct vec
{
   double lon, lat;
   double x, y, z;
}
Vec;


struct sptConstraints
{
   int  status;
   char errorMsg       [1024];
   char indexConstraint[32768];
   char geomConstraint [32768];
};


struct sptConstraints sptConeSearch(char *indname, int indexMode, int indexEncoding, int level, 
                                    char *xcol, char *ycol, char *zcol, double ra, double dec, double radius);

struct sptConstraints sptPolygonSearch(char *indname, int indexMode, int indexEncoding, int level, 
                                       char *xcol, char *ycol, char *zcol, int npoly, double *ra, double *dec);

int    sptSortCmpFunc   (const void *a, const void *b);

int    hpxConeSearch    (int omax, double ra, double dec, double radius, struct hpxRanges **hpxrng);
int    hpxPolygonSearch (int omax, int npoly, double *ra, double *dec,   struct hpxRanges **hpxrng);
double hpxMaxPixRad     (int order);
void   hpxBoundingCircle(int np, Vec *point, Vec *center, double *cosrad);
void   hpxGetCircle     (int np, Vec *point, int q, Vec *center, double *cosrad);
void   hpxPix2Loc       (int order, int64_t pix, double *z, double *phi);
int    hpxCompressBits  (int64_t v);

void   vCopy     (Vec *v, Vec *c);
void   vCalcRADec(Vec *v);
void   vCalcXYZ  (Vec *v);
void   vMidpoint (Vec *a, Vec *b, Vec *c);
void   vPixCenter(Vec *a, Vec *b, Vec *c, Vec *d, Vec *v);
int    vCross    (Vec *a, Vec *b, Vec *c);
double vDot      (Vec *a, Vec *b);
double vNormalize(Vec *a);
void   vReverse  (Vec *v);
void   vPrint    (Vec *v, char *lbl);

int sptDebug;

#endif
