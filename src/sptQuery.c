/*
  Copyright (c) 2020, Caltech IPAC.
  This code is released with a BSD 3-clause license. License information is at
    https://github.com/Caltech-IPAC/nexsciTAP/blob/master/LICENSE
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sptQuery.h>

void printError(char *msg);


int main(int argc, char **argv)
{
   int     i, indxMode, srchMode, level, npoly;

   double  rac, decc;
   double  rad;
   double *ra, *dec;

   char    indname       [1024];

   char    indxModeStr   [32];
   char    srchModeStr   [32];
   char    msg           [1024];

   char   *end;


   struct sptConstraints constraints;


   // Process generic command-line arguments

   if(argc > 2 && strcmp(argv[1], "-d") == 0)
   {
      sptDebug = 1;

      ++argv;
      --argc;
   }

   if(argc < 7)
      printError("Usage: sptQuery [-d] HTM|HPX level -c ra dec radius | sptQuery [-d] HTM|HPX level -p ra1 dec1 ... raN decN (N >= 3)");

   strcpy(indxModeStr, argv[1]);

   if(strcmp(indxModeStr, "HTM") == 0)
   {
      indxMode = HTM;
      strcpy(indname, "htmind");
   }

   else if(strcmp(indxModeStr, "HPX") == 0)
   {
      indxMode = HPX;
      strcpy(indname, "hpxind");
   }
   else
   {
      sprintf(msg, "Invalid index type '%s' (must be HTM or HPX)", indxModeStr);
      printError(msg);
   }

   level = atoi(argv[2]);

   strcpy(srchModeStr, argv[3]);

   srchMode = CONE;
   if(strcmp(srchModeStr, "-p") == 0)
      srchMode = POLYGON;


   if(srchMode == CONE)   // CONE search arguments
   {
      rac  = atof(argv[4]);
      decc = atof(argv[5]);
      rad  = atof(argv[6]);

      constraints = sptConeSearch(indname, indxMode, DECIMAL, level, "x", "y", "z", rac, decc, rad);

      if(constraints.status)
         printError(constraints.errorMsg);
   }

   else   // POLYGON search arguments
   {
      npoly = (argc - 4)/2;

      if(npoly < 3 || (npoly * 2 + 4) != argc)
         printError("Usage: sptQuery [-d] HTM|HPX level -c ra dec radius | sptQuery [-d] HTM|HPX level -p ra1 dec1 ... raN decN (N >= 3)");

      ra  = (double *)malloc(npoly * sizeof(double));
      dec = (double *)malloc(npoly * sizeof(double));

      if(sptDebug)
      {
         printf("\nnpoly = %d\n\n", npoly);
         fflush(stdout);
      }


      for(i=0; i<npoly; ++i)
      {
         ra[i] = strtod(argv[2*i+4], &end);

         if(end < argv[2*i+4]+strlen(argv[2*i+4]))
         {
            sprintf(msg, "Invalid RA (deg) string: [%s]", argv[2*i+4]);
            printError(msg);
         }

         while(ra[i] <    0.) ra[i] += 360.;
         while(ra[i] >= 360.) ra[i] -= 360.;


         dec[i] = strtod(argv[2*i+5], &end);

         if(end < argv[2*i+5]+strlen(argv[2*i+5]))
         {
            sprintf(msg, "Invalid Dec (deg) string: [%s]", argv[2*i+5]);
            printError(msg);
         }

         if(dec[i] < -90 || dec[i] > 90.)
         {
            sprintf(msg, "Invalid Dec (deg) value: [%s]", argv[2*i+5]);
            printError(msg);
         }
      }

      constraints = sptPolygonSearch(indname, indxMode, DECIMAL, level, "x", "y", "z", npoly, ra, dec);

      if(constraints.status)
         printError(constraints.errorMsg);
   }


   printf("\nINDEX CONSTRAINT> %s\n\n", constraints.indexConstraint);
   printf("GEOM  CONSTRAINT> %s\n\n", constraints.geomConstraint);

   exit(0);
}


void printError(char *msg)
{
   printf("[struct stat=\"ERROR\", msg=\"%s\"]\n", msg);
   fflush(stdout);
   exit(1);
}
